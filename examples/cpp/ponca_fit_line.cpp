/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <numeric> // transform_reduce
#include <algorithm>
#include <iostream>
#include <vector>

#include <Ponca/Fitting>

using namespace Eigen;
using namespace std;
using namespace Ponca;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 
// This class defines the input data format
class MyPoint
{
public:
  enum {Dim = 3};
  typedef double Scalar;
  typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
  typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;
 
  PONCA_MULTIARCH inline MyPoint(const VectorType& _pos    = VectorType::Zero())
    : m_pos(_pos) {}
 
  PONCA_MULTIARCH inline const VectorType& pos()    const { return m_pos; }
 
  PONCA_MULTIARCH inline VectorType& pos()    { return m_pos; }
 
private:
  VectorType m_pos;
};
typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;

typedef DistWeightFunc<MyPoint,ConstantWeightKernel<Scalar> > WeightFunc;
typedef Basket<MyPoint,WeightFunc, CovarianceLineFit> Fit;


int main(int argc, char **argv) {

    int n = 10000;

    // Generate a random set of n points along a line
    vector<MyPoint> points(n);
    VectorType dir = VectorType::Random().normalized();
    std::generate(points.begin(), points.end(), [&dir]() {
        return (dir * Eigen::internal::random<Scalar>(0.1,2)).eval();
    });
    const VectorType& p = points.at(0).pos();

    std::cout << "====================\nLeastSquareLineFit:\n";
    cout << "Direction of the generated line: " << dir.transpose() << endl;

    // Fit line on data
    Fit _fit;
    _fit.setWeightFunc(WeightFunc(p));
    _fit.compute(points.cbegin(), points.cend());

    // Check Fit output
    if( _fit.isStable() ) {
        cout << "Direction of the fitted 3D line: " << _fit.direction().transpose() << endl;
        cout << "Origin of the fitted line: " << _fit.origin().transpose() << endl;
        cout << "Projection of the basis center on the fitted line: " << _fit.project(p).transpose() << endl;

        /// Compute sum of the distances between samples and line: should be 0
        Scalar dist = std::transform_reduce(points.cbegin(), points.cend(),
                                                   Scalar(0),
                                                   std::plus<>(),
                                                   [&_fit](const MyPoint &q) {
                                                       return (q.pos() - _fit.project(q.pos())).norm();
                                                   });
        cout << "Mean error between samples and fitted line: " << dist / Scalar(n) << endl;
        return EXIT_SUCCESS;
    }
    cerr << "Fit is not stable: " << _fit.getCurrentState() << endl;
    return EXIT_FAILURE;
}
