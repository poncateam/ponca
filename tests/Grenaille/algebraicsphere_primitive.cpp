/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


/*!
    \file test/Grenaille/okabe_primitive.cpp
    \brief Test validity Plane Primitive
 */

#include "../common/testing.h"
#include "../common/testUtils.h"

#include <vector>


using namespace std;
using namespace Grenaille;

// This class defines the input data format
template<typename _Scalar, int _Dim>
class MyPoint
{
public:
    enum {Dim = _Dim};
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1, Eigen::DontAlign>   VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim, Eigen::DontAlign> MatrixType;

    MULTIARCH inline MyPoint(   const VectorType &pos    = VectorType::Zero(),
                                const VectorType& normal = VectorType::Zero())
        : m_pos(pos), m_normal(normal) {}

    MULTIARCH inline const VectorType& pos()    const { return m_pos; }
    MULTIARCH inline const VectorType& normal() const { return m_normal; }

    MULTIARCH inline VectorType& pos()    { return m_pos; }
    MULTIARCH inline VectorType& normal() { return m_normal; }

    static inline MyPoint Random(Scalar radius)
    {
        VectorType scale;
        scale.setRandom();
        VectorType n = VectorType::Random().normalized();
        VectorType p =   scale.asDiagonal() * n * radius   // create an ellipse
                         * Eigen::internal::random<Scalar>(Scalar(0.99),Scalar(1.01)); // add high frequency noise

        n = (scale.asDiagonal().inverse() * n).normalized();

        return MyPoint (p, n);
    };


private:
    VectorType m_pos, m_normal;
};

template<typename Point, typename Fit, typename WeightFunc>
void testFunction()
{
        // Define related structure
    typedef typename Point::Scalar Scalar;
    typedef typename Point::VectorType VectorType;

     // dummy precision is to hard for this test
    Scalar epsilon = 1e-4; //Eigen::NumTraits<Scalar>::dummy_precision();

    // generate sample data
    int n = Eigen::internal::random<int>(10,100);
    Scalar radius = Eigen::internal::random<Scalar>(1,10);
    Scalar tmax = Scalar(10.) * std::sqrt(Scalar(4. * M_PI) * radius * radius/n);
    vector<Point> vecs (n);

    for(int k=0; k<n; ++k)
    {
        vecs[k] = Point::Random(radius);
    }


#pragma omp parallel for
    for(unsigned int k=0; k<vecs.size(); ++k)
    {
        Fit fit;
        fit.setWeightFunc(WeightFunc(tmax));

        fit.init(vecs[k].pos());
        fit.compute(vecs.cbegin(), vecs.cend());

        if(fit.isStable())
        {
            Fit f2 = fit;
            // do to iterations to test successive calls
            VectorType candidate = VectorType::Random();
            for (unsigned int j = 0; j != 2; ++j){
                f2.changeBasis(fit.basisCenter() + Scalar(0.01)*VectorType::Random());

                VERIFY( Eigen::internal::isMuchSmallerThan(
                            (fit.project(candidate) - f2.project(candidate)).norm(),
                            Scalar(1.),
                            epsilon) );
            }
        }
    }
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef MyPoint<Scalar, Dim> Point;

    // We test only primitive functions and not the fitting procedure
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightFunc;
    typedef Basket<Point, WeightFunc, OrientedSphereFit> Sphere;

    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, Sphere, WeightFunc>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test Algebraic Sphere Primitive functions in 3 dimensions..." << endl;
    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
    cout << "Ok..." << endl;
}
