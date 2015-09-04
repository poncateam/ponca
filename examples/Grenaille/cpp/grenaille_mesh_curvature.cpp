/*
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


/*!
\file examples/Grenaille/cpp/grenaille_mesh_curvature.cpp
\brief Basic use of Grenaille
*/
#include <limits>
#include <queue>
#include <iostream>
#include <fstream>

#include "Patate/common/surface_mesh/surfaceMesh.h"
#include "Patate/common/surface_mesh/objReader.h"
#include "Patate/common/surface_mesh/objWriter.h"
#include "Patate/grenaille.h"

#include <vector>

using namespace std;
using namespace Grenaille;


typedef double Scalar;
typedef Eigen::Matrix<Scalar, 3, 1> Vector;

class Mesh : public PatateCommon::SurfaceMesh {
public:
    typedef ::Scalar Scalar;
    typedef ::Vector Vector;

private:
    typedef PatateCommon::SurfaceMesh::VertexProperty<Scalar> ScalarProp;
    typedef PatateCommon::SurfaceMesh::VertexProperty<Vector> VectorProp;

public:
    Mesh() {
        m_positions    = addVertexProperty<Vector>("v:position",     Vector::Zero());
        m_normals      = addVertexProperty<Vector>("v:normals",      Vector::Zero());
        m_minDirection = addVertexProperty<Vector>("v:minDirection", Vector::Zero());
        m_maxDirection = addVertexProperty<Vector>("v:maxDirection", Vector::Zero());

        m_minCurvature = addVertexProperty<Scalar>("v:minCurvature", 0);
        m_maxCurvature = addVertexProperty<Scalar>("v:maxCurvature", 0);
        m_geomVar      = addVertexProperty<Scalar>("v:geomVar",      0);
    }

    template <typename Derived>
    inline Vertex addVertex(const Eigen::DenseBase<Derived>& pos) {
        Vertex vx = PatateCommon::SurfaceMesh::addVertex();
        m_positions[vx] = pos;
        return vx;
    }

    inline const Vector& position(Vertex vx) const
    { return m_positions[vx]; }

    inline Vector& position(Vertex vx)
    { return m_positions[vx]; }

    inline const Vector& normal(Vertex vx) const
    { return m_normals[vx]; }

    inline Vector& normal(Vertex vx)
    { return m_normals[vx]; }

    inline const Vector& minDirection(Vertex vx) const
    { return m_minDirection[vx]; }

    inline Vector& minDirection(Vertex vx)
    { return m_minDirection[vx]; }

    inline const Vector& maxDirection(Vertex vx) const
    { return m_maxDirection[vx]; }

    inline Vector& maxDirection(Vertex vx)
    { return m_maxDirection[vx]; }

    inline Scalar minCurvature(Vertex vx) const
    { return m_minCurvature[vx]; }

    inline Scalar& minCurvature(Vertex vx)
    { return m_minCurvature[vx]; }

    inline Scalar maxCurvature(Vertex vx) const
    { return m_maxCurvature[vx]; }

    inline Scalar& maxCurvature(Vertex vx)
    { return m_maxCurvature[vx]; }

    inline Scalar geomVar(Vertex vx) const
    { return m_geomVar[vx]; }

    inline Scalar& geomVar(Vertex vx)
    { return m_geomVar[vx]; }


    void computeNormals() {
        for(VertexIterator vit = verticesBegin();
            vit != verticesEnd(); ++vit) {
            normal(*vit) = Vector::Zero();
        }
        for(FaceIterator fit = facesBegin();
            fit != facesEnd(); ++fit) {

            VertexAroundFaceCirculator vit = vertices(*fit);
            const Vector& p0 = position(*(++vit));
            const Vector& p1 = position(*(++vit));
            const Vector& p2 = position(*(++vit));

            Vector n = (p1 - p0).cross(p2 - p0);

            VertexAroundFaceCirculator vEnd = vit;
            do {
                normal(*vit) += n;
                ++vit;
            } while(vit != vEnd);
        }
        for(VertexIterator vit = verticesBegin();
            vit != verticesEnd(); ++vit) {
            normal(*vit).normalize();
        }
    }


private:
    VectorProp m_positions;
    VectorProp m_normals;
    VectorProp m_minDirection;
    VectorProp m_maxDirection;
    ScalarProp m_minCurvature;
    ScalarProp m_maxCurvature;
    ScalarProp m_geomVar;
};

typedef Mesh::Vertex Vertex;

// This class defines the input data format
class GLSPoint
{
public:
    enum {Dim = 3};
    typedef ::Scalar                         Scalar;
    typedef Eigen::Matrix<Scalar, Dim, 1>    VectorType;
    typedef Eigen::Matrix<Scalar, Dim, Dim>  MatrixType;

    MULTIARCH inline GLSPoint(Mesh& mesh, Vertex vx)
        : m_mesh(&mesh), m_vx(vx) {}

    MULTIARCH inline const VectorType& pos()    const { return m_mesh->position(m_vx); }
    MULTIARCH inline VectorType& pos()    { return m_mesh->position(m_vx); }

    MULTIARCH inline const VectorType& normal() const { return m_mesh->normal(m_vx); }
    MULTIARCH inline VectorType& normal() { return m_mesh->normal(m_vx); }

private:
    Mesh*   m_mesh;
    Vertex  m_vx;
};


typedef Grenaille::DistWeightFunc<GLSPoint, Grenaille::SmoothWeightKernel<Scalar> >
               WeightFunc;
typedef Basket<GLSPoint, WeightFunc, OrientedSphereFit, GLSParam,
               OrientedSphereScaleSpaceDer, GLSDer, GLSCurvatureHelper, GLSGeomVar> Fit;


Scalar fit(Mesh& mesh, Scalar radius) {
    // compute curvature
    Vertex v0, v1;
    Mesh::VertexAroundVertexCirculator vc, vc_end;

    std::queue<Vertex> nrings;

    Scalar r = Scalar(2) * radius;
    Scalar gCurvMax = Scalar(0);

    for(Mesh::VertexIterator vit = mesh.verticesBegin();
        vit != mesh.verticesEnd(); ++vit)
    {
        v0 = *vit;

        // initialize fit
        Fit fit;
        fit.init(mesh.position(v0));
        nrings.push(v0);

        fit.setWeightFunc(WeightFunc(r));

        int nbNeighboors = 0;

        // process stack
        if(mesh.halfedge(v0).isValid()) {
            Mesh::VertexProperty<bool> processed =
                    mesh.addVertexProperty<bool>("v:processed",false);
            processed[v0] = true;
            while(!nrings.empty()) {
                v1 = nrings.front();
                nrings.pop();
                // add to the fit
                fit.addNeighbor(GLSPoint(mesh, v1));
                nbNeighboors++;

                // add its 1-ring neighboors to the stack if they haven't been processed yet...
                vc = mesh.vertices(v1);
                vc_end = vc;
                do {
                    if(!processed[*vc]) {
                        processed[*vc] = true;
                        // ...and if they are close enough
                        if((mesh.position(*vc) - mesh.position(v0)).norm() < r)
                            nrings.push(*vc);
                    }
                } while(++vc != vc_end);
            }
            mesh.removeVertexProperty<bool>(processed);
        }

        fit.finalize();

        if(fit.isReady()) {
            if(!fit.isStable())
                std::cerr << "Unstable fit: " << nbNeighboors << " neighbhoors\n";

            mesh.maxDirection(v0)  = fit.GLSk1Direction();
            mesh.minDirection(v0)  = fit.GLSk2Direction();
            mesh.maxCurvature(v0)  = fit.GLSk1();
            mesh.minCurvature(v0)  = fit.GLSk2();
            mesh.geomVar(v0)       = fit.geomVar();

            Scalar gCurv = fabs(fit.GLSGaussianCurvature());
            if(gCurv > gCurvMax)
                gCurvMax = gCurv;
        } else {
            std::cerr << "Fit impossible: " << nbNeighboors << " neighbhoors\n";
            mesh.maxDirection(v0)  = Vector::Zero();
            mesh.minDirection(v0)  = Vector::Zero();
            mesh.maxCurvature(v0)  = std::numeric_limits<Scalar>::quiet_NaN();
            mesh.minCurvature(v0)  = std::numeric_limits<Scalar>::quiet_NaN();
            mesh.geomVar(v0)       = std::numeric_limits<Scalar>::quiet_NaN();
        }
    }

    return gCurvMax;
}


Vector getColor(Scalar _value, Scalar _valueMin, Scalar _valueMax) {
    Vector c = Vector(1.0, 1.0, 1.0);
    Scalar dv;
    if(_value == 0.)
        return c;
    if(isnan(_value))
        return Vector(0.0, 1.0, 0.0);
    if (_value < _valueMin)
        _value = _valueMin;
    if (_value > _valueMax)
        _value = _valueMax;
    dv = _valueMax - _valueMin;
    if(_value < (_valueMin + 0.5 * dv)){
        c(0) = 2 * (_value - _valueMin) / dv;
        c(1) = 2 * (_value - _valueMin) / dv;
        c(2) = 1;
    }else{
        c(2) = 2 - 2 * (_value - _valueMin) / dv;
        c(1) = 2 - 2 * (_value - _valueMin) / dv;
        c(0) = 1;
    }
    return c;
}


int main(int argc, char** argv)
{
    if(argc != 3) {
        std::cerr << "Usage: " << argv[0] << " OBJ_IN OBJ_OUT\n";
        return 1;
    }

    const char* inFilename  = argv[1];
    const char* outFilename = argv[2];
    double radius = 0.089;

    Mesh mesh;
    {
        std::ifstream in(inFilename);
        PatateCommon::OBJReader<Mesh>().read(in, mesh);
    }

    mesh.triangulate();
    mesh.computeNormals();

    Vector meshBarycenter = Vector::Zero();
    Scalar meshRadius = 0.f;

    Mesh::VertexIterator vit, vend = mesh.verticesEnd();
    for(vit = mesh.verticesBegin(); vit != vend; ++vit)
    {
        meshBarycenter += mesh.position(*vit);
    }
    meshBarycenter /= mesh.nVertices();

    for (vit = mesh.verticesBegin(); vit != vend; ++vit)
    {
        Scalar dist = (mesh.position(*vit) - meshBarycenter).norm();
        if(dist > meshRadius)
            meshRadius = dist;
    }

    Scalar gCurveMax = fit(mesh, meshRadius * radius);
    std::cout << "gCurveMax: " << gCurveMax << "\n";

    {
        std::ofstream out(outFilename);
//        PatateCommon::OBJWriter<Mesh>().write(out, mesh);

        mesh.garbageCollection();
        for(Mesh::VertexIterator vit = mesh.verticesBegin();
            vit != mesh.verticesEnd(); ++vit)
        {
            Scalar curv = mesh.minCurvature(*vit) * mesh.maxCurvature(*vit);
            Vector color = getColor(curv, -gCurveMax, gCurveMax);
            out << "v " << mesh.position(*vit).transpose()
                << " " << color.transpose() << "\n";
        }
        for(Mesh::VertexIterator vit = mesh.verticesBegin();
            vit != mesh.verticesEnd(); ++vit)
        {
            out << "vn " << mesh.normal(*vit).transpose() << "\n";
        }
        for(Mesh::FaceIterator fit = mesh.facesBegin();
            fit != mesh.facesEnd(); ++fit)
        {
            out << "f";
            Mesh::VertexAroundFaceCirculator vit  = mesh.vertices(*fit);
            Mesh::VertexAroundFaceCirculator vEnd = vit;
            do {
                out << " " << (*vit).idx() + 1
                    << "//" << (*vit).idx() + 1;
                ++vit;
            } while(vit != vEnd);
            out << "\n";
        }
    }
}
