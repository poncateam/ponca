/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
\file examples/Grenaille/ssgls.cu
\brief Screen space GLS using c++/CUDA
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <FreeImagePlus.h>
#include <vector>

#include "Eigen/Core"
#include "Patate/grenaille.h"

//! [mypoint]
class ScreenSpacePoint
{
public:
    enum {Dim = 3};
    typedef float Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
    typedef Eigen::Matrix<Scalar, 2,   1>   ScreenVectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

    MULTIARCH inline ScreenSpacePoint(const VectorType       &_pos    = VectorType::Zero(),
                                      const VectorType       &_normal = VectorType::Zero(),
                                      const ScreenVectorType &_spos   = ScreenVectorType::Zero())
        : m_pos(_pos), m_normal(_normal), m_spos(_spos){}

    MULTIARCH inline const VectorType& pos()	const { return m_pos; }
    MULTIARCH inline const VectorType& normal()	const { return m_normal; }
    MULTIARCH inline const ScreenVectorType& spos() const { return m_spos; }

    MULTIARCH inline VectorType& pos()	 { return m_pos; }
    MULTIARCH inline VectorType& normal()	 { return m_normal; }
    MULTIARCH inline ScreenVectorType& spos() { return m_spos; }

private:
    ScreenVectorType m_spos;
    VectorType	m_pos, m_normal;
};
//! [mypoint]

typedef ScreenSpacePoint::Scalar Scalar;
typedef ScreenSpacePoint::VectorType VectorType;
typedef ScreenSpacePoint::ScreenVectorType ScreenVectorType;

//! [w_def]
class ProjectedWeightFunc: public Grenaille::DistWeightFunc<ScreenSpacePoint,Grenaille::SmoothWeightKernel<Scalar> >
{
public:
    typedef ScreenSpacePoint::Scalar Scalar;
    typedef ScreenSpacePoint::VectorType VectorType;

    MULTIARCH inline ProjectedWeightFunc(const Scalar& _t = Scalar(1.), const Scalar _dz = 0.f)
        : Grenaille::DistWeightFunc<ScreenSpacePoint,Grenaille::SmoothWeightKernel<Scalar> >(_t),
          m_dz(_dz) {}

    MULTIARCH inline Scalar w(const VectorType& _relativePos, const ScreenSpacePoint&  _attributes) const
    {
        Scalar d  = _attributes.spos().norm();
        const float dz = abs(_relativePos[2]);
        if (d > m_t || (m_dz != Scalar(0) && dz > m_dz))
        {
            return Scalar(0.);
        }
        return m_wk.f(d/m_t);
    }
private:
    float m_dz;
};
//! [w_def]

//! [fit_def]
typedef Grenaille::Basket< ScreenSpacePoint,
                           ProjectedWeightFunc,
                           Grenaille::OrientedSphereFit,
                           Grenaille::GLSParam> ScreenSpaceFit;
//! [fit_def]

//! [data_acces]
__device__ int getId(const int _x,
                     const int _y,
                     const int _width,
                     const int _height,
                     const int _component,
                     const int _nbComponent)
{
    return (_component) + _nbComponent*(_x + _y * _width);
}

__device__ VectorType getVector(const int _x,
                                const int _y,
                                const int _width,
                                const int _height,
                                const float* _buffer)
{
    VectorType r;
    r << Scalar(_buffer[getId(_x,_y,_width,_height,0,3)]),
        Scalar(_buffer[getId(_x,_y,_width,_height,1,3)]),
        Scalar(_buffer[getId(_x,_y,_width,_height,2,3)]);
    return r;
}
//! [data_acces]

//! [kernel]
__global__ void doGLS_kernel(	int _imgw, int _imgh, int _scale,
								float _maxDepthDiff, float* _positions, float* _normals,
								float* _result)
{
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int bw = blockDim.x;
	int bh = blockDim.y;
	int x = blockIdx.x * bw + tx;
	int y = blockIdx.y * bh + ty;

	int idx = y * _imgw + x;

	if((x >= _imgw || y >= _imgh))
	{
		return;
	}
	else if(getVector(x, y, _imgw, _imgh, _normals).squaredNorm() == 0.f)
	{
        _result[idx] = 0.f;
		return;
	}

	VectorType one = VectorType::Ones();
	const float scale2 = float(_scale * _scale);

	ScreenSpaceFit fit;
	fit.init(getVector(x, y, _imgw, _imgh, _positions) * 2.f - one);
	fit.setWeightFunc(ProjectedWeightFunc(_scale, _maxDepthDiff));

	_result[idx] = 0.f;

	// collect neighborhood
	for(int dy = -_scale; dy != _scale + 1; dy++)
	{
		for(int dx = -_scale; dx != _scale + 1; dx++)
		{
			float dist2 = dy*dy + dx*dx;
			// Check if we are in the circular screen-space neighborhood
			if (dist2 < scale2)
			{
				int nx, ny; // neighbor ids

				nx = x + dx;
				ny = y + dy;

				// Check image boundaries
				if(nx >= 0 && ny >= 0 && nx < _imgw && ny < _imgh)
				{
					ScreenSpacePoint::VectorType n = getVector(nx, ny, _imgw, _imgh, _normals);
					// add nei only when the normal is properly defined
					if(n.squaredNorm() != 0.f)
					{
						// RGB to XYZ remapping
						n =  2.f * n - one;
						n.normalize();

						ScreenSpacePoint::ScreenVectorType xyCoord;
						xyCoord[0] = dx;
						xyCoord[1] = dy;

						ScreenSpacePoint::VectorType p = getVector(nx, ny, _imgw, _imgh, _positions) * 2.f - one;
						// GLS computation
						fit.addNeighbor(ScreenSpacePoint(p, n, xyCoord));
					}
				}
			}
		}
	}

	// closed form minimization
	fit.finalize();
	_result[idx] = fit.kappa();
}
//! [kernel]

/**
* \brief RGB basic color representation
*/
typedef struct
{
    double r,g,b;
}Color;

/**
* \brief Return Color corresponding to the _value param. Simulating a "seismic" like color map
*/
__host__ Color getColor(float _value, float _valueMin, float _valueMax)
{
    Color c = {1.0, 1.0, 1.0};
    double dv;

    // Unknown values in our kernel
    if(_value == 0.)
    {
        return c;
    }

    // Threshold
    if (_value < _valueMin)
    {
        _value = _valueMin;
    }

    if (_value > _valueMax)
    {
        _value = _valueMax;
    }

    // Interval
    dv = _valueMax - _valueMin;

    // Seismic color map like
    if(_value < (_valueMin + 0.5 * dv))
    {
        c.r = 2 * (_value - _valueMin) / dv;
        c.g = 2 * (_value - _valueMin) / dv;
        c.b = 1;
    }
    else
    {
        c.b = 2 - 2 * (_value - _valueMin) / dv;
        c.g = 2 - 2 * (_value - _valueMin) / dv;
        c.r = 1;
    }

    return c;
}

/**
* \brief Load input images with freeimageplus lib
*/
__host__ bool loadImages(fipImage& _positions, fipImage& _normals, const char* _positionsFilename, const char* _normalsFilename)
{
    if(!_positions.load(_positionsFilename))
    {
        fprintf(stderr, "Cannot load _positions\n");
        return 0;
    }

    if(!_normals.load(_normalsFilename))
    {
        fprintf(stderr, "Cannot load _normal map\n");
        return 0;
    }

    _positions.convertTo24Bits();
    _normals.convertTo24Bits();

    return 1;
}

/**
* \brief Init input datas to be used on host
*/
__host__ bool initInputDatas(const fipImage& _positions, const fipImage& _normals, float** _positionsInfos, float** _normalsInfos,
                             unsigned int& _width, unsigned int& _height)
{
    BYTE* positionsPixels = 0;
    positionsPixels = _positions.accessPixels();
    if(!positionsPixels)
    {
        fprintf(stderr, "Cannot get _positions datas\n");
        return 0;
    }

    BYTE* normalsPixels = 0;
    normalsPixels = _normals.accessPixels();
    if(!normalsPixels)
    {
        fprintf(stderr, "Cannot get normals datas\n");
        return 0;
    }

    _width = _positions.getWidth();
    _height = _positions.getHeight();

    (*_positionsInfos) = new float[_width*_height*3];
    (*_normalsInfos) = new float[_width*_height*3];
    if(!*_positionsInfos || !*_normalsInfos)
    {
        fprintf(stderr, "Cannot alloc memory in initInputDatas\n");
        return 0;
    }

    for(int i = 0; i < _width * _height; ++i)
    {
        (*_positionsInfos)[i * 3 + 0] = positionsPixels[i * 3 + 0] / 255.f * 2.f - 1.f;
        (*_positionsInfos)[i * 3 + 1] = positionsPixels[i * 3 + 1] / 255.f * 2.f - 1.f;
        (*_positionsInfos)[i * 3 + 2] = positionsPixels[i * 3 + 2] / 255.f * 2.f - 1.f;

        (*_normalsInfos)[i * 3 + 0] = normalsPixels[i * 3 + 0] / 255.f;
        (*_normalsInfos)[i * 3 + 1] = normalsPixels[i * 3 + 1] / 255.f;
        (*_normalsInfos)[i * 3 + 2] = normalsPixels[i * 3 + 2] / 255.f;
    }

    positionsPixels = 0;
    normalsPixels = 0;

    return 1;
}

/**
* \brief Save _results into png image
*/
__host__ bool saveResult(float* _results, const unsigned int& _width, const unsigned int& _height,
                         const char* _positionsFilename, const char* _resultFilename)
{
    float kappaMin = *std::min_element(_results, _results + _width*_height);
    float kappaMax = *std::max_element(_results, _results + _width*_height);
    std::cout << "Kappa min : " << kappaMin << std::endl;
    std::cout << "Kappa max : " << kappaMax << std::endl;

    fipImage result;
    if(!result.load(_positionsFilename))
    {
        fprintf(stderr, "Cannot load positions\n");
        return 0;
    }

    result.convertTo24Bits();

    BYTE* resultInfos = 0;
    resultInfos = result.accessPixels();
    if(!resultInfos)
    {
        fprintf(stderr, "Cannot get result datas\n");
        return 0;
    }

    for(int i = 0; i < _width * _height; ++i)
    {
        //check nan
        if(_results[i] != _results[i])
        {
            _results[i] = 0.f;
        }

        Color c = getColor(_results[i], -10., 10.);

        resultInfos[i * 3 + 0] = c.r * 255.;
        resultInfos[i * 3 + 1] = c.g * 255.;
        resultInfos[i * 3 + 2] = c.b * 255.;
        //resultInfos[i * 4 + 3] = 255.;
    }

    if(!result.save(_resultFilename, 0))
    {
        fprintf(stderr, "Cannot save image\n");
    }

    resultInfos = 0;
    result.clear();

    return 1;
}

__host__ int adjust(int n, int blockSize)
{
   if (n < blockSize) { return n; }
   return (n / blockSize + (n % blockSize == 0 ? 0 : 1)) * blockSize;
}

int main()
{
    std::string positionsFilename = "./data/ssgls_sample_wc.png";
    std::string normalsFilename = "./data/ssgls_sample_normal.png";
    std::string resultFilename = "./data/ssgls_sample_results.png";

    fipImage positions, normals;

    if(!loadImages(positions, normals, positionsFilename.c_str(), normalsFilename.c_str()))
    {
        return 0;
    }

    float fScale = 10.f;
    float fMaxDepthDiff = 0.f;
    unsigned int width = 0;
    unsigned int height = 0;
    float* positionsInfos = 0;
    float* normalsInfos = 0;

    if(!initInputDatas(positions, normals, &positionsInfos, &normalsInfos, width, height))
    {
        return 0;
    }

    std::cout << "Image size : " << width << "*" << height << std::endl;

    /*********** Init Output ************/
    float *results = new float[width*height];
    for(int i = 0; i < width * height; ++i)
    {
        results[i] = 0.f;
    }

    /************* Init device mem *************/
    size_t sizeResults = width * height * sizeof(float);
    size_t sizeImg = width * height * 3 * sizeof(float);

    float* positionsInfos_device;
    float* normalsInfos_device;
    float* results_device;

    cudaMalloc(&positionsInfos_device, sizeImg);
    cudaMemcpy(positionsInfos_device, positionsInfos, sizeImg, cudaMemcpyHostToDevice);

    cudaMalloc(&normalsInfos_device, sizeImg);
    cudaMemcpy(normalsInfos_device, normalsInfos, sizeImg, cudaMemcpyHostToDevice);

    cudaMalloc(&results_device, sizeResults);
    cudaMemcpy(results_device, results, sizeResults, cudaMemcpyHostToDevice);

    cudaError_t err = cudaGetLastError();
    /************* Memory conf *************/

    // calculate grid size
    dim3 block(16, 16, 1);
    dim3 grid(adjust(width, block.x) / block.x, adjust(height, block.y) / block.y, 1);

    /************* Kernel Call *************/

    std::cout << "ssCurvature running..." << std::endl;

    doGLS_kernel<<<grid, block>>>(width, height, fScale, fMaxDepthDiff, positionsInfos_device, normalsInfos_device, results_device);

    cudaThreadSynchronize();	// Wait for the GPU launched work to complete
    err = cudaGetLastError();

    std::cout << "ssCurvature completed..." << std::endl;

    /************* Get Results *************/
    cudaMemcpy(results, results_device, sizeResults, cudaMemcpyDeviceToHost);

    err = cudaGetLastError();

    std::cout << "Finalizing..." << std::endl;

    /********** Cuda Free ************/
    cudaFree(positionsInfos_device);
    cudaFree(normalsInfos_device);
    cudaFree(results_device);

    err = cudaGetLastError();

    /********** Saving _result ************/
    if(!saveResult(results, width, height, positionsFilename.c_str(), resultFilename.c_str()))
    {
        return 0;
    }

    /********** Free Memory *********/
    positions.clear();
    normals.clear();

    delete [] positionsInfos;
    delete [] normalsInfos;
    delete [] results;

    cudaDeviceReset();
    err = cudaGetLastError();

    std::cout << "Finished !" << std::endl;

    return 0;
}

