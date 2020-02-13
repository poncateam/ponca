/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
\file examples/Grenaille/basic_cpu.h
\brief Basic use of Grenaille

\author: Nicolas Mellado, Gautier Ciaudo
*/
#include <cmath>
#include <algorithm>
#include <iostream>

#include "ponca.h"
#include "Eigen/Eigen"

#include <vector>

using namespace std;
using namespace Ponca;

#define DIMENSION 3

/*
   \brief Variant of the MyPoint class allowing to work with external raw data.

   Using this approach, ones can use the patate library with already existing
   data-structures and without any data-duplication.

   In this example, we use this class to map an interlaced raw array containing
   both point normals and coordinates.
 */
class MyPoint
{
public:
    enum {Dim = DIMENSION};
    typedef double Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

    MULTIARCH inline MyPoint(Scalar* _interlacedArray, int _pId)
        : m_pos   (Eigen::Map< const VectorType >(_interlacedArray + Dim*2*_pId  )),
        m_normal(Eigen::Map< const VectorType >(_interlacedArray + Dim*2*_pId+Dim))
    {}

    MULTIARCH inline const Eigen::Map< const VectorType >& pos()    const { return m_pos; }
    MULTIARCH inline const Eigen::Map< const VectorType >& normal() const { return m_normal; }

private:
    Eigen::Map< const VectorType > m_pos, m_normal;
};

typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;

// Define related structure
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc;
typedef Basket<MyPoint,WeightFunc,OrientedSphereFit,   GLSParam> Fit;


template<typename Fit>
void test_fit(Fit& _fit,
              Scalar* _interlacedArray,
              int _n,
              const VectorType& _p)
{
    Scalar tmax = 100.0;

    // Set a weighting function instance
    _fit.setWeightFunc(WeightFunc(tmax));

    // Set the evaluation position
    _fit.init(_p);

    // Iterate over samples and _fit the primitive
    // A MyPoint instance is generated on the fly to bind the raw arrays to the
    // library representation. No copy is done at this step.
    for(int i = 0; i!= _n; i++)
    {
        _fit.addNeighbor(MyPoint(_interlacedArray, i));
    }

    //finalize fitting
    _fit.finalize();

    //Test if the fitting ended without errors
    if(_fit.isStable())
    {
        cout << "Center: [" << _fit.center().transpose() << "] ;  radius: " << _fit.radius() << endl;

        cout << "Pratt normalization"
            << (_fit.applyPrattNorm() ? " is now done." : " has already been applied.") << endl;

        // Play with fitting output
        cout << "Value of the scalar field at the initial point: "
            << _p.transpose()
            << " is equal to " << _fit.potential(_p)
            << endl;

        cout << "It's gradient at this place is equal to: "
            << _fit.primitiveGradient(_p).transpose()
            << endl;

        cout << "Fitted Sphere: " << endl
            << "\t Tau  : "      << _fit.tau()             << endl
            << "\t Eta  : "      << _fit.eta().transpose() << endl
            << "\t Kappa: "      << _fit.kappa()           << endl;

        cout << "The initial point " << _p.transpose()              << endl
            << "Is projected at   " << _fit.project(_p).transpose() << endl;
    }
}

// Build an interlaced array containing _n position and normal vectors
Scalar* buildInterlacedArray(int _n)
{
    Scalar* interlacedArray = new Scalar[uint(2*DIMENSION*_n)];

    for(int k=0; k<_n; ++k)
    {
        // For the simplicity of this example, we use Eigen Vectors to compute
        // both coordinates and normals, and then copy the raw values to an
        // interlaced array, discarding the Eigen representation.
        Eigen::Matrix<Scalar, DIMENSION, 1> nvec = Eigen::Matrix<Scalar, DIMENSION, 1>::Random().normalized();
        Eigen::Matrix<Scalar, DIMENSION, 1> pvec = nvec * Eigen::internal::random<Scalar>(0.9,1.1);

        // Grab coordinates and store them as raw buffer
        memcpy(interlacedArray+2*DIMENSION*k,           pvec.data(), DIMENSION*sizeof(Scalar));
        memcpy(interlacedArray+2*DIMENSION*k+DIMENSION, nvec.data(), DIMENSION*sizeof(Scalar));

    }

    return interlacedArray;
}

int main()
{
    // Build arrays containing normals and positions, simulating data coming from
    // outside the library.
    int n = 1000;
    Scalar *interlacedArray = buildInterlacedArray(n);

    // set evaluation point and scale at the first coordinate
    VectorType p (interlacedArray);

    // Here we now perform the fit, starting from a raw interlaced buffer, without
    // any data duplication
    Fit fit;
    test_fit(fit, interlacedArray, n, p);
}
