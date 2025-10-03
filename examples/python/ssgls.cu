/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 
 Compile using 
 nvcc --ptx -I ../../Patate/ -I /path/to/eigen-nvcc/ ssgls.cu
*/


#include <cuda.h>
#include <algorithm>

#include "Eigen/Core"
#include "Patate/grenaille.h"


//! [mypoint]
class MyPoint
{
public:
    enum {Dim = 3};
    typedef float Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;
    typedef Eigen::Matrix<Scalar, 2, 1>   ScreenVectorType;

    MULTIARCH inline MyPoint(   const VectorType& _pos        = VectorType::Zero(),
                                const VectorType& _normal     = VectorType::Zero(),
                                const ScreenVectorType& _spos = ScreenVectorType::Zero(),
                                const Scalar _dz = 0.f)
        : m_pos(_pos), m_normal(_normal), m_spos(_spos), m_dz(_dz){}

    MULTIARCH inline const VectorType& pos()	const { return m_pos; }
    MULTIARCH inline const VectorType& normal()	const { return m_normal; }
    MULTIARCH inline const ScreenVectorType& spos() const { return m_spos; }
    MULTIARCH inline const float & dz()	const { return m_dz; }


    MULTIARCH inline VectorType& pos()	 { return m_pos; }
    MULTIARCH inline VectorType& normal()	 { return m_normal; }
    MULTIARCH inline ScreenVectorType& spos() { return m_spos; }
    MULTIARCH inline float& dz()	 { return m_dz; }


private:
    ScreenVectorType m_spos;
    VectorType	m_pos, m_normal;
    float m_dz; // depth threshold
};
//! [mypoint]

typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;
typedef MyPoint::ScreenVectorType ScreenVectorType;

//! [w_def]
class ProjectWeightFunc: public Grenaille::DistWeightFunc<MyPoint, Grenaille::SmoothWeightKernel<Scalar> >
{
public:
    typedef MyPoint::Scalar Scalar;
    typedef MyPoint::VectorType VectorType;

    /*
    Default constructor (needed by Grenaille). Note that the screenspace
    evaluation position is specified as parameter
    */
    MULTIARCH inline ProjectWeightFunc( const Scalar& _t                = 1.f,
                                        const ScreenVectorType& _refPos = ScreenVectorType::Zero(),
                                        const Scalar& _dz               = 0.f)
        : Grenaille::DistWeightFunc<MyPoint, Grenaille::SmoothWeightKernel<Scalar> >(_t), m_refPos(_refPos), m_dz(_dz) {}

    MULTIARCH inline Scalar w(const VectorType& _q, const MyPoint&  _attributes) const
    {
        Scalar d  = (_attributes.spos()-m_refPos).norm();
        const Scalar dz = _attributes.dz();

        if (d > m_t || dz > m_dz)
            return Scalar(0.);

        return m_wk.f(d/m_t);
    }
private:
    ScreenVectorType m_refPos;
    float m_dz;
};
//! [w_def]

//! [fit_def]
typedef Grenaille::Basket<MyPoint,ProjectWeightFunc,Grenaille::OrientedSphereFit, Grenaille::GLSParam> Gls;
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
extern "C" { 
__global__ void doGLS_kernel(int* _params, //[w, h, scale]
                             float* _positions,
                             float* _normals,
                             float* _result)
{
    
    uint x = (blockIdx.x * blockDim.x) + threadIdx.x;
    uint y = (blockIdx.y * blockDim.y) + threadIdx.y;
    
    const int &width     = _params[0];
    const int &height    = _params[1];

    if (x < width && y < height)
    {
        const int &scale = _params[2];

        ScreenVectorType refPos;
        refPos << x, y;

        int dx, dy; // neighbor offset ids
        int nx, ny; // neighbor ids

        Gls gls;
        gls.setNeighborFilter({scale, refPos});
        gls.init( getVector(x,y,width,height,_positions) );

        if (getVector(x,y,width,height,_normals).squaredNorm() == 0.f )
        {
            _result[getId(x,y,width,height,0,1)] = -1.0;
        }
        else
        {
            //_result[getId(x,y,width,height,0,1)] = getVector(x,y,width,height,_normals)(0);
            VectorType p, n;

            // collect neighborhood
            VectorType one = VectorType::Zero();

            for(dy = -scale; dy != scale; dy++)
            {
                for(dx = -scale; dx != scale; dx++)
                {
                    nx = x+dx;
                    ny = y+dy;


                    // Check image boundaries
                    if (nx >= 0 && ny >= 0 && nx < width && ny < height)
                    {
                        n = getVector(nx,ny,width,height,_normals);

                        // add nei only when the _normal is properly defined
                        // need to use an explicit floating point comparison with pycuda
                        if (n.squaredNorm() != 0.f )
                        {

                            // RGB to XYZ remapping
                            n =  2.f * n - one;
                            n.normalize();

                            // GLS computation
                            gls.addNeighbor(MyPoint(getVector(nx,ny,width,height,_positions),
                                n,
                                ScreenVectorType(nx,ny)));
                        }
                    }
                }
            }
            // closed form minimization
            gls.finalize();
            _result[getId(x,y,width,height,0,1)] = gls.kappa();
        }
    }
}

}
//! [kernel]

