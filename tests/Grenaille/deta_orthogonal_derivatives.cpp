/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
 \file tets/Grenaille/deta_orthogonal_derivatives.cpp
 \brief Test validity of the spatial differentiation of eta

 \author: Nicolas Mellado
 */
#include <cmath>
#include <algorithm>
#include <iostream>

#include "Patate/grenaille.h"
#include "Eigen/Eigen"

#include <vector>

using namespace std;
using namespace Grenaille;

#define SCALAR float

SCALAR radius = 10;
SCALAR tmax   = 2;


// This class defines the input data format
class MyPoint{
public:
  enum {Dim = 3};
  typedef SCALAR Scalar;
  typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
  typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;

  MULTIARCH inline MyPoint(const VectorType &pos    = VectorType::Zero(), 
		 const VectorType& normal = VectorType::Zero())
    : _pos(pos), _normal(normal) {}
    
  MULTIARCH inline const VectorType& pos()    const { return _pos; }  
  MULTIARCH inline const VectorType& normal() const { return _normal; }

  MULTIARCH inline VectorType& pos()    { return _pos; }  
  MULTIARCH inline VectorType& normal() { return _normal; }

  static inline MyPoint Random() {
    VectorType scale; scale << 3, 1, 2;
    VectorType n = VectorType::Random().normalized();
    VectorType p =   scale.asDiagonal() * n* radius   // create an ellipse
                   * Eigen::internal::random<Scalar>(0.99,1.01); // add high frequency noise
    n = (scale.asDiagonal().inverse() * n).normalized();
    return MyPoint (p, n);
  };


private:
  VectorType _pos, _normal;
};


// Define related structure
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<MyPoint::Scalar> > WeightFunc; 

// Add unused arbitrary extensions
typedef Basket<MyPoint,WeightFunc,OrientedSphereFit, GLSParam, OrientedSphereScaleSpaceDer, GLSDer, GLSGeomVar> Fit;

template<typename Fit>
bool
test_orthoEta(Fit& fit, vector<MyPoint>& vecs){
  
  MyPoint::Scalar epsilon = std::numeric_limits<MyPoint::Scalar>::epsilon();
    
  for(int k=0; k<vecs.size(); ++k){  
    fit.setWeightFunc(WeightFunc(tmax));  
    
    fit.init(vecs[k].pos());
    for(vector<MyPoint>::iterator it = vecs.begin(); it != vecs.end(); it++)
      fit.addNeighbor(*it);      
    fit.finalize();
  
    typename Fit::VectorType eta  = fit.eta();
    typename Fit::MatrixType deta = fit.deta().template middleCols<MyPoint::Dim>(fit.isScaleDer() ? 1: 0);

    // FIXME Use eigen colwise test
    if (eta.dot(deta.col(0)) > epsilon ||
        eta.dot(deta.col(1)) > epsilon ||
        eta.dot(deta.col(2)) > epsilon )
      return false;    
  }
  return true;
}

int main() {
  cout << "Deta Orthogonal Derivatives Test .... " << flush;
  
  int n = 10000;
  vector<MyPoint> vecs (n);

  for(int k=0; k<n; ++k)
    vecs[k] = MyPoint::Random();
    
  Fit fit;
  if (test_orthoEta(fit, vecs)){
    cout << "OK" << endl;
    return EXIT_SUCCESS;
  }
  else{
    cout << "FAILED" << endl;
    return EXIT_FAILURE;
  }
}
