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

template<typename DataPoint, typename Fit, typename WeightFunc> 
void testFunction()
{
    // Define related structure
    typedef typename DataPoint::Scalar Scalar;
    typedef typename DataPoint::VectorType VectorType;
    typedef typename DataPoint::HVectorType HVectorType;
    
    Fit f;
    f.m_p = HVectorType::Random();
    f.m_p.template head<DataPoint::Dim>().normalize();
    
    VectorType query = VectorType::Random();
    
    Scalar epsilon = testEpsilon<Scalar>();
    
    // Check if we get the same point when projecting points x and x+gradient.
    VectorType proj1 = f.project( query );
    VectorType proj2 = f.project( query + f.primitiveGradient( query ) );
    
    VERIFY( ( proj1 - proj2 ).norm() <= epsilon);

    // Check that the potential value is equal to the distance between the query
    // and its projection    
    VERIFY( std::abs( std::abs(f.potential(query)) - (query - proj1).norm()) <= epsilon);

    // Check that the potential value of a projected point is equal to 0
    VERIFY( std::abs( f.potential(proj1) ) <= epsilon);
}

template<typename Scalar, int Dim>
void callSubTests()
{
    typedef PointPosistionNormal<Scalar, Dim> Point;

    // We test only primitive functions and not the fitting procedure
    typedef DistWeightFunc<Point, SmoothWeightKernel<Scalar> > WeightFunc;
    typedef Basket<Point, WeightFunc, Plane> Plane;        
    
    for(int i = 0; i < g_repeat; ++i)
    {
        CALL_SUBTEST(( testFunction<Point, Plane, WeightFunc>() ));
    }
}

int main(int argc, char** argv)
{
    if(!init_testing(argc, argv))
    {
        return EXIT_FAILURE;
    }

    cout << "Test Plane Primitive functions in 3 dimensions..." << endl;
    callSubTests<float, 3>();
    callSubTests<double, 3>();
    callSubTests<long double, 3>();
    cout << "Ok..." << endl;

    cout << "Test Plane Primitive functions in 4 dimensions..." << endl;
    callSubTests<float, 4>();
    callSubTests<double, 4>();
    callSubTests<long double, 4>();
    cout << "Ok..." << endl;

    cout << "Test Plane Primitive functions in 5 dimensions..." << endl;
    callSubTests<float, 5>();
    callSubTests<double, 5>();
    callSubTests<long double, 5>();
    cout << "Ok..." << endl;

    cout << "Test Plane Primitive functions in 42 dimensions..." << endl;
    callSubTests<float, 42>();
    callSubTests<double, 42>();
    callSubTests<long double, 42>();
    cout << "Ok..." << endl;
}
