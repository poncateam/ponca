/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*!
\file examples/Ponca/ssgls.cu
\brief Screen space GLS using c++/CUDA
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <chrono>

#include <png.h>

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int

#include <Ponca/src/Fitting/basket.h>
#include <Ponca/src/Fitting/gls.h>
#include <Ponca/src/Fitting/orientedSphereFit.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>



/**************************************************************************************************/
/* IO (source: http://zarb.org/~gc/html/libpng.html )                                             */
/**************************************************************************************************/

class PNGImage
{
public:
  inline bool load(const char *file_name);
  inline bool loaded () const { return ! row_pointers.empty(); }
  inline bool save(const char *file_name);

  inline png_uint_32 width()  const { return m_width; };
  inline png_uint_32 height() const { return m_height; };

  inline const std::vector<png_bytep>& buffer() const { return row_pointers; }
  inline std::vector<png_bytep>& buffer() { return row_pointers; }
  inline png_byte colorType() const { return png_get_color_type(png_ptr, info_ptr);}

  ~PNGImage() { for (auto e: row_pointers) delete e; row_pointers.clear(); }
private:
  png_uint_32 m_width, m_height;
  png_byte color_type;
  png_byte bit_depth;

  png_structp png_ptr;
  png_infop info_ptr;
  int number_of_passes;
  std::vector<png_bytep> row_pointers;

  using vecSizeT = typename std::vector<png_bytep>::size_type;
};

bool
PNGImage::load(const char* file_name)
{
    unsigned char header[8];    // 8 is the maximum size that can be checked

    /* open file and test for it being a png */
    FILE *fp = fopen(file_name, "rb");
    if (!fp)
    {
        std::cerr << "[read_png_file] File " \
                  <<  file_name
                  << " could not be opened for reading"
                  << std::endl;
        return false;
    }

    fread(header, 1, 8, fp);
    if (png_sig_cmp(header, 0, 8))
    {
        std::cerr << "[read_png_file] File " \
                  <<  file_name
                  << " is not recognized as a PNG file"
                  << std::endl;
        return false;
    }


    /* initialize stuff */
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

    if (!png_ptr)
    {
        std::cerr << "[read_png_file] png_create_read_struct failed"
                  << std::endl;
        return false;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        std::cerr << "[read_png_file] png_create_info_struct failed"
                  << std::endl;
        return false;
    }

    if (setjmp(png_jmpbuf(png_ptr)))
    {
        std::cerr << "[read_png_file] Error during init_iod"
                  << std::endl;
        return false;
    }

    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);

    png_read_info(png_ptr, info_ptr);

    m_width = png_get_image_width(png_ptr, info_ptr);
    m_height = png_get_image_height(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr, info_ptr);
    bit_depth = png_get_bit_depth(png_ptr, info_ptr);

    number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);


    /* read file */
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        std::cerr << "[read_png_file] Error during read_image"
                  << std::endl;
        return false;
    }

    row_pointers.resize( m_height );
    for (vecSizeT y=0; y< vecSizeT(m_height); y++)
      row_pointers[y] = (png_byte*) (malloc(png_get_rowbytes(png_ptr,info_ptr)));

    png_read_image(png_ptr, row_pointers.data());

    fclose(fp);

    return true;
}

bool
PNGImage::save(const char* file_name) {
  /* create file */
  FILE *fp = fopen(file_name, "wb");
  if (!fp)
  {
      std::cerr << "[write_png_file] File " \
                <<  file_name
                << " could not be opened for reading"
                << std::endl;
      return false;
  }


  /* initialize stuff */
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

  if (!png_ptr)
  {
      std::cerr << "[write_png_file] png_create_write_struct failed"
                << std::endl;
      return false;
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  {
      std::cerr << "[write_png_file] png_create_info_struct failed"
                << std::endl;
      return false;
  }

  if (setjmp(png_jmpbuf(png_ptr)))
  {
      std::cerr << "[write_png_file] Error during init_io"
                << std::endl;
      return false;
  }

  png_init_io(png_ptr, fp);


  /* write header */
  if (setjmp(png_jmpbuf(png_ptr)))
  {
      std::cerr << "[write_png_file] Error during writing header"
                << std::endl;
      return false;
  }

  png_set_IHDR(png_ptr, info_ptr, m_width, m_height,
               bit_depth, color_type, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_write_info(png_ptr, info_ptr);


  /* write bytes */
  if (setjmp(png_jmpbuf(png_ptr)))
  {
      std::cerr << "[write_png_file] Error during writing bytes"
                << std::endl;
      return false;
  }

  png_write_image(png_ptr, row_pointers.data());


  /* end write */
  if (setjmp(png_jmpbuf(png_ptr)))
  {
      std::cerr << "[write_png_file] Error during end of write"
                << std::endl;
      return false;
  }


  png_write_end(png_ptr, nullptr);

  fclose(fp);
  return true;
}

/**************************************************************************************************/
/* Ponca initialization                                                                           */
/**************************************************************************************************/
//! [mypoint]
class ScreenSpacePoint
{
public:
    enum {Dim = 3};
    typedef float Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
    typedef Eigen::Matrix<Scalar, 2,   1>   ScreenVectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

    PONCA_MULTIARCH inline ScreenSpacePoint(const VectorType       &_pos    = VectorType::Zero(),
                                      const VectorType       &_normal = VectorType::Zero(),
                                      const ScreenVectorType &_spos   = ScreenVectorType::Zero())
        : m_pos(_pos), m_normal(_normal), m_spos(_spos){}

    PONCA_MULTIARCH inline const VectorType& pos()	const { return m_pos; }
    PONCA_MULTIARCH inline const VectorType& normal()	const { return m_normal; }
    PONCA_MULTIARCH inline const ScreenVectorType& spos() const { return m_spos; }

    PONCA_MULTIARCH inline VectorType& pos()	 { return m_pos; }
    PONCA_MULTIARCH inline VectorType& normal()	 { return m_normal; }
    PONCA_MULTIARCH inline ScreenVectorType& spos() { return m_spos; }

private:
    VectorType	m_pos, m_normal;
    ScreenVectorType m_spos;
};
//! [mypoint]

typedef ScreenSpacePoint::Scalar Scalar;
typedef ScreenSpacePoint::VectorType VectorType;
typedef ScreenSpacePoint::ScreenVectorType ScreenVectorType;

//! [w_def]
class ProjectedWeightFunc: public Ponca::DistWeightFunc<ScreenSpacePoint,Ponca::SmoothWeightKernel<Scalar> >
{
public:
    typedef ScreenSpacePoint::Scalar Scalar;
    typedef ScreenSpacePoint::VectorType VectorType;
    typedef Ponca::DistWeightFunc<ScreenSpacePoint,Ponca::SmoothWeightKernel<Scalar> >::WeightReturnType WeightReturnType;


    PONCA_MULTIARCH inline ProjectedWeightFunc(const Scalar& _t = Scalar(1.), const Scalar _dz = 0.f)
        : Ponca::DistWeightFunc<ScreenSpacePoint,Ponca::SmoothWeightKernel<Scalar> >(_t),
          m_dz(_dz) {}

    PONCA_MULTIARCH inline WeightReturnType w(const VectorType& _relativePos, const ScreenSpacePoint&  _attributes) const
    {
        PONCA_MULTIARCH_STD_MATH(abs);
        Scalar d  = _attributes.spos().norm();
        const float dz = abs(_relativePos[2]);
        if (d > m_t || (m_dz != Scalar(0) && dz > m_dz))
        {
            return {Scalar(0.), _relativePos};
        }
        return {m_wk.f(d/m_t), _relativePos};
    }
private:
    float m_dz;
};
//! [w_def]

//! [fit_def]
typedef Ponca::Basket< ScreenSpacePoint,
                           ProjectedWeightFunc,
                           Ponca::OrientedSphereFit,
                           Ponca::GLSParam> ScreenSpaceFit;
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
__global__ void doGLS_kernel( int _imgw, int _imgh, int _scale,
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

//    VectorType vvvvv = getVector(x, y, _imgw, _imgh, _positions);
//    VectorType nnnnn = getVector(x, y, _imgw, _imgh, _normals);
//    _result[idx] = vvvvv(2);
//    return;

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
* \brief Init input datas to be used on host
*/
__host__ bool initInputDatas(const PNGImage& positions, const PNGImage& normals,
                             std::vector<float>& positionsInfos,
                             std::vector<float>& normalsInfos,
                             unsigned int& width, unsigned int& height)
{

    if (positions.colorType() != PNG_COLOR_TYPE_RGB) {
      std::cerr << "[process_file] color_type of input file must be PNG_COLOR_TYPE_RGB ("
                << PNG_COLOR_TYPE_RGB
                << ") (is "
                << positions.colorType()
                << ")"
                << std::endl;
      return false;
    }


    width = positions.width();
    height = positions.height();

    positionsInfos.resize(width*height*3);
    normalsInfos.resize(width*height*3);

    auto pbuf = positions.buffer();
    auto nbuf = normals.buffer();

    for (int j = 0; j < height; ++j) {
        png_bytep pcol = pbuf[j];
        png_bytep ncol = nbuf[j];

        float* pout = positionsInfos.data()+j*width*3;
        float* nout = normalsInfos.data()+j*width*3;

        auto scaleValues = [](const png_byte& in){ return in / 255.f * 2.f - 1.f; };
        std::transform(pcol, pcol+width*3, pout, scaleValues );
        std::transform(ncol, ncol+width*3, nout, scaleValues );
    }

    return true;
}

/**
* \brief Save _results into png image
*/
__host__ bool saveResult(float* _results,
                         const char* _positionsFilename, const char* _resultFilename)
{

    PNGImage result;
    if(!result.load(_positionsFilename))
    {
        fprintf(stderr, "Cannot load positions\n");
        return false;
    }

    int width = result.width();
    int height = result.height();

    auto pbuf = result.buffer().data();

    for (int j = 0; j < height; ++j) {
        float* pin = _results+j*width;
        png_bytep col = pbuf[j];
        for (int i = 0; i < width; ++i) {
            //check nan
            if(std::isnan(pin[i]))
            {
                pin[i] = 0.f;
            }
            Color c = getColor(pin[i], -10., 10.);

            col[i * 3 + 0] = c.r * 255.;
            col[i * 3 + 1] = c.g * 255.;
            col[i * 3 + 2] = c.b * 255.;
        }
    }

    if(!result.save(_resultFilename))
    {
        fprintf(stderr, "Cannot save image\n");
    }

    return true;
}

__host__ int adjust(int n, int blockSize)
{
   if (n < blockSize) { return n; }
   return (n / blockSize + (n % blockSize == 0 ? 0 : 1)) * blockSize;
}

int main()
{
    const char *positionsFilename = "./data/ssgls_sample_wc.png";
    const char *normalsFilename = "./data/ssgls_sample_normal.png";
    const char *resultFilename = "./ssgls_results.png";

    PNGImage positions, normals;

    if(!positions.load(positionsFilename) || ! normals.load(normalsFilename))
    {
        return 0;
    }

    float fScale = 10.f;
    float fMaxDepthDiff = 0.00f;
    unsigned int width = 0;
    unsigned int height = 0;
    std::vector<float> positionsInfos, normalsInfos;

    if(!initInputDatas(positions, normals, positionsInfos, normalsInfos, width, height))
    {
        return 0;
    }

    std::cout << "Image size : " << width << "*" << height << std::endl;

    /*********** Init Output ************/
    float *results = new float[width*height];
    std::fill( results, results + width*height, 0.f );

    /************* Init device mem *************/
    size_t sizeResults = width * height * sizeof(float);
    size_t sizeImg = width * height * 3 * sizeof(float);

    float* positionsInfos_device;
    float* normalsInfos_device;
    float* results_device;

    cudaMalloc(&positionsInfos_device, sizeImg);
    cudaMemcpy(positionsInfos_device, positionsInfos.data(), sizeImg, cudaMemcpyHostToDevice);

    cudaMalloc(&normalsInfos_device, sizeImg);
    cudaMemcpy(normalsInfos_device, normalsInfos.data(), sizeImg, cudaMemcpyHostToDevice);

    cudaMalloc(&results_device, sizeResults);
    cudaMemcpy(results_device, results, sizeResults, cudaMemcpyHostToDevice);

    cudaError_t err = cudaGetLastError();
    /************* Memory conf *************/

    // calculate grid size
    dim3 block(32, 32, 1);
    dim3 grid(adjust(width, block.x) / block.x, adjust(height, block.y) / block.y, 1);

    /************* Kernel Call *************/

    std::cout << "ssCurvature running..." << std::endl;

    // dry run: first call is always slower
    doGLS_kernel<<<grid, block>>>(width, height, fScale, fMaxDepthDiff, positionsInfos_device, normalsInfos_device, results_device);

    int nbrun = 100;
    auto start = std::chrono::system_clock::now();
    for( int i = 0; i != nbrun; ++i) {
      doGLS_kernel<<<grid, block>>>(width, height, fScale, fMaxDepthDiff, positionsInfos_device, normalsInfos_device, results_device);
      cudaDeviceSynchronize();	// Wait for the GPU launched work to complete
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = (end-start)/double(nbrun);

    err = cudaGetLastError();

    std::cout << "ssCurvature completed in " << diff.count() << " s" << std::endl;

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
    if(!saveResult(results, positionsFilename, resultFilename))
    {
        return 0;
    }

    /********** Free Memory *********/
    delete [] results;

    cudaDeviceReset();
    err = cudaGetLastError();

    std::cout << "Finished !" << std::endl;

    return 0;
}

