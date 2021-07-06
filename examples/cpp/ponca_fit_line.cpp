/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include<stdio.h>
#include <random>

#include <Ponca/Fitting>

#include "Eigen/Eigen"
#include <Eigen/Core> 
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
using namespace Ponca;

#define DIMENSION 3
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This class defines the input data format
class MyPoint
{
public:
  enum {Dim = 3};
  typedef double Scalar;
  typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
  typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;
 
  PONCA_MULTIARCH inline MyPoint(const VectorType& _pos    = VectorType::Zero(),
                                 const VectorType& _normal = VectorType::Zero())
    : m_pos(_pos), m_normal(_normal) {}
 
  PONCA_MULTIARCH inline const VectorType& pos()    const { return m_pos; }
  PONCA_MULTIARCH inline const VectorType& normal() const { return m_normal; }
 
  PONCA_MULTIARCH inline VectorType& pos()    { return m_pos; }
  PONCA_MULTIARCH inline VectorType& normal() { return m_normal; }
 
  static inline MyPoint Random()
  {
    VectorType n = VectorType::Random().normalized();
    VectorType p = n * Eigen::internal::random<Scalar>(0.9,1.1);
    return MyPoint (p, (n + VectorType::Random()*0.1).normalized());
  }
 
private:
  VectorType m_pos, m_normal;
};
typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;

typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc;
typedef Basket<MyPoint,WeightFunc, LeastSquareLine> fit;


int main(int argc, char **argv) {

   
    int n = 10000;
    vector<MyPoint> points(n);
    std::generate(points.begin(), points.end(), []() {return MyPoint::Random(); });
    const VectorType& p = points.at(0).pos();

    /*===========================================================*/
    std::cout << "====================\nLeastSquareLineFit:\n";
    
    Scalar tmax = 0.1;
    fit _fit;

    /* Set a weighting function instance   */
    _fit.setWeightFunc(WeightFunc(tmax));

    /* Set the evaluation position */
    _fit.init(p);

    for( int idx = 0; idx < points.size(); idx++ )
    {
     
        _fit.addNeighbor( points.at(idx) );
    }

    _fit.finalize();

    if( _fit.isStable() )
    {
         cout << "\nA point on the fitted 3D line: \n"
            << _fit.point()    
            << endl;

        cout << "\nThe direction of the fitted 3D line: \n"
            << _fit.direction()
            << endl;

      

        
    } 

    return 0;   
}
