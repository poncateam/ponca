#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include<stdio.h>
#include "happly.h"
#include <random>
#include <string>
#include <iterator>
//Ponca
#include <Ponca/Fitting>
#include "Eigen/Eigen"
#include <Ponca/src/SpatialPartitioning/KdTree/kdTree.h>


// Polyscope
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"

#include <Eigen/Core> 
#include <Eigen/Dense>
using namespace Eigen;

using namespace std;
using namespace Ponca;


#define DIMENSION 3
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 

class MyPoint
{
public:
    enum {Dim = DIMENSION};
    typedef double Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1>   VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;
 
    PONCA_MULTIARCH inline MyPoint(const std::array<Scalar, 3>&poss,
                                   const std::array<Scalar, 3>&  norm)
        : m_pos    (Eigen::Map< const VectorType >(poss.begin())),
          m_normal (Eigen::Map< const VectorType >(norm.begin()))
    {}
 
    PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& pos()    const { return m_pos; }
    PONCA_MULTIARCH inline const Eigen::Map< const VectorType >& normal() const { return m_normal; }
 
private:
    Eigen::Map< const VectorType > m_pos, m_normal;
}; 


typedef MyPoint::Scalar Scalar;
typedef MyPoint::VectorType VectorType;

void loadPointCloud(std::string filename,
                    std::vector<std::array<double, 3>>& vertexPositionsOut) {

  happly::PLYData plyIn(filename);

  /* Get mesh-style data from the object */
  vertexPositionsOut = plyIn.getVertexPositions();


}
typedef DistWeightFunc<MyPoint,SmoothWeightKernel<Scalar> > WeightFunc;
typedef Basket<MyPoint,WeightFunc, LeastSquareLine> fit;





int main(int argc, char **argv) {

   
    polyscope::init();

    string filename = "line.ply";

    std::ifstream testStream(filename);
    if (!testStream) {
        return 0;
    }

    std::vector< std::array<double, 3> > positions;
    /* Load positions from file */
    loadPointCloud(filename, positions);
    testStream.close();
    vector<MyPoint> points;
    /*===========================================================*/
    std::cout << "====================\nLeastSquareLineFit:\n";
    
    for(const auto &p : positions)
    {
       
        points.push_back({p , {0,0,0}});
    }

    const VectorType& p = points.at(0).pos();
     
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
         cout << "A point on the fitted 3D line: "
            << _fit.point()    
            << endl;

        cout << "The direction of the fitted 3D line: "
            << _fit.direction()
            << endl;

        /* Adding edges to view on polyscope */
        std::vector<VectorType> nodes ({_fit.point(), _fit.direction()});
        std::vector<std::array<size_t, 2>> edges ({{0,1}});

        /* Add the curve network */
        polyscope::registerCurveNetwork("my network", nodes, edges);
        
    } 
    


    /* visualize! */
    polyscope::registerPointCloud("positions", positions);

    /* Show the gui */
    polyscope::show(); 
    return 0;   
}
