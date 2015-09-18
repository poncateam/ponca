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

#ifdef __linux
#include "time.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

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

    MULTIARCH inline GLSPoint()
        : m_pos(), m_norm() {}
    MULTIARCH inline GLSPoint(const Mesh& mesh, Vertex vx)
        : m_pos(mesh.position(vx)), m_norm(mesh.normal(vx)) {}

    MULTIARCH inline const VectorType& pos()    const { return m_pos; }
    MULTIARCH inline const VectorType& normal() const { return m_norm; }

private:
    Vector m_pos;
    Vector m_norm;
};


typedef Grenaille::DistWeightFunc<GLSPoint, Grenaille::SmoothWeightKernel<Scalar> >
               WeightFunc;
typedef Basket<GLSPoint, WeightFunc, OrientedSphereFit, GLSParam,
               OrientedSphereScaleSpaceDer, GLSDer, GLSCurvatureHelper, GLSGeomVar> Fit;


class Grid {
public:
    typedef Eigen::Array3i Cell;
    typedef Eigen::AlignedBox3i CellBlock;
    typedef std::vector<unsigned> IndexVector;
    typedef std::vector<GLSPoint, Eigen::aligned_allocator<GLSPoint> > PointVector;

public:
    Grid(const Mesh& mesh, Scalar cellSize)
        : _cellSize(cellSize),
          _gridSize(),
          _gridBase() {

        // Compute bounding box.
        typedef Eigen::AlignedBox<Scalar, 3> BoundingBox;
        BoundingBox bb;
        for(Mesh::VertexIterator vit = mesh.verticesBegin();
            vit != mesh.verticesEnd(); ++vit) {
            bb.extend(mesh.position(*vit));
        }

        _gridBase = bb.min();
        _gridSize = gridCoord(bb.max()) + Cell::Constant(1);
        std::cout << "GridSize: " << _gridSize.transpose()
                  << " (" << nCells() << ")\n";

        // Compute cells size (number of points in each cell)
        _cellsIndex.resize(nCells() + 1, 0);
        for(Mesh::VertexIterator vit = mesh.verticesBegin();
            vit != mesh.verticesEnd(); ++vit) {
            unsigned i = cellIndex(mesh.position(*vit));
            ++_cellsIndex[i + 1];
        }

        // Prefix sum to get indices
        for(unsigned i = 1; i <= nCells(); ++i) {
            _cellsIndex[i] += _cellsIndex[i - 1];
        }

        // Fill the point array
        IndexVector pointCount(nCells(), 0);
        _points.resize(_cellsIndex.back());
        for(Mesh::VertexIterator vit = mesh.verticesBegin();
            vit != mesh.verticesEnd(); ++vit) {
            unsigned i = cellIndex(mesh.position(*vit));
            *(_pointsBegin(i) + pointCount[i]) = GLSPoint(mesh, *vit);
            ++pointCount[i];
        }

        for(unsigned i = 0; i < nCells(); ++i) {
            assert(pointCount[i] == _cellsIndex[i+1] - _cellsIndex[i]);
        }
    }

    inline unsigned nCells() const {
        return _gridSize.prod();
    }
    inline Cell gridCoord(const Vector& p) const {
        return Eigen::Array3i(((p - _gridBase) / _cellSize).cast<int>());
    }
    inline int cellIndex(const Cell& c) const {
        return c(0) + _gridSize(0) * (c(1) + _gridSize(1) * c(2));
    }
    inline int cellIndex(const Vector& p) const {
        Cell c = gridCoord(p);
        assert(c(0) >= 0 && c(1) >= 0 && c(2) >= 0
            && c(0) < _gridSize(0) && c(1) < _gridSize(1) && c(2) < _gridSize(2));
        return c(0) + _gridSize(0) * (c(1) + _gridSize(1) * c(2));
    }

    inline const GLSPoint* pointsBegin(unsigned ci) const {
        assert(ci < nCells());
        return &_points[0] + _cellsIndex[ci];
    }
    inline const GLSPoint* pointsEnd(unsigned ci) const {
        assert(ci < nCells());
        return &_points[0] + _cellsIndex[ci + 1];
    }

    CellBlock cellsAround(const Vector& p, Scalar radius) {
        return CellBlock(gridCoord(p - Vector::Constant(radius))
                                 .max(Cell(0, 0, 0)).matrix(),
                         (gridCoord(p + Vector::Constant(radius)) + Cell::Constant(1))
                                 .min(_gridSize).matrix());
    }

private:
    inline GLSPoint* _pointsBegin(unsigned ci) {
        assert(ci < nCells());
        return &_points[0] + _cellsIndex[ci];
    }

private:
    Scalar      _cellSize;
    Cell        _gridSize;
    Vector      _gridBase; // min corner
    IndexVector _cellsIndex;
    PointVector _points;
};


Scalar fit(Mesh& mesh, Scalar radius) {
    // compute curvature
    const Scalar r = Scalar(2) * radius;
    const Scalar r2 = r * r;
    const unsigned nVertices = mesh.verticesSize();

    Scalar gCurvMax = Scalar(0);

    Grid grid(mesh, r);

#ifdef _OPENMP
#pragma omp parallel for reduction(max:gCurvMax)
#endif
    for(unsigned i = 0; i < nVertices; ++i)
    {
        Vertex v0 = Vertex(i);
        if(mesh.isDeleted(v0)) {
            continue;
        }

        // initialize fit
        Fit fit;
        Vector p0 = mesh.position(v0);
        fit.init(p0);
        fit.setWeightFunc(WeightFunc(r));

        int nbNeighboors = 0;

        Eigen::AlignedBox3i cells = grid.cellsAround(p0, r);
        Grid::Cell c = cells.min().array();
        for(c(2) = cells.min()(2); c(2) < cells.max()(2); ++c(2)) {
            for(c(1) = cells.min()(1); c(1) < cells.max()(1); ++c(1)) {
                for(c(0) = cells.min()(0); c(0) < cells.max()(0); ++c(0)) {
                    const GLSPoint* p = grid.pointsBegin(grid.cellIndex(c));
                    const GLSPoint* end = grid.pointsEnd(grid.cellIndex(c));
                    for(; p != end; ++p) {
                        if((p->pos() - p0).squaredNorm() < r2) {
                            fit.addNeighbor(*p);
                            ++nbNeighboors;
                        }
                    }
                }
            }
        }

        fit.finalize();

        if(fit.isReady()) {
            if(!fit.isStable()) {
#ifdef _OPENMP
#pragma omp critical(cerr)
#endif
                {
                    std::cerr << "Unstable fit: " << nbNeighboors << " neighbhoors\n";
                }
            }

            mesh.maxDirection(v0)  = fit.GLSk1Direction();
            mesh.minDirection(v0)  = fit.GLSk2Direction();
            mesh.maxCurvature(v0)  = fit.GLSk1();
            mesh.minCurvature(v0)  = fit.GLSk2();
            mesh.geomVar(v0)       = fit.geomVar();

            Scalar gCurv = fabs(fit.GLSGaussianCurvature());
            if(gCurv > gCurvMax)
                gCurvMax = gCurv;
        } else {
#ifdef _OPENMP
#pragma omp critical(cerr)
#endif
            {
                std::cerr << "Fit impossible: " << nbNeighboors << " neighbhoors\n";
            }
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


Vector getColor2(Scalar value, Scalar p) {
    bool sign = value < 0;
    value = 2 * std::min(std::max(std::pow(std::abs(value), p), Scalar(0)), Scalar(1));
    return Vector(std::abs(1 - value),
                  std::min(2 - value, Scalar(1)),
                  Scalar(sign));
}


bool noopErrorCallback(const std::string& /*msg*/, unsigned /*line*/, void* /*ptr*/) {
    return false;
}


void usage(const char* progName) {
    std::cerr << "Usage: " << progName << " [-r RADIUS] [-m MAX_CURV] [-R] OBJ_IN OBJ_OUT\n";
    exit(1);
}


int main(int argc, char** argv)
{
    double radius   = 0;
    bool   relative = false;
    const char* inFilename  = 0;
    const char* outFilename = 0;
    double maxCurv = 0;

    for(int i = 1; i < argc; ++i) {
        if(argv[i][0] == '-') {
            if(argv[i] == std::string("-r")
            || argv[i] == std::string("--radius")) {
                if(++i == argc) {
                    std::cerr << "Error: " << argv[i-1] << " need an argument.\n";
                    usage(argv[0]);
                }
                if(radius != 0) std::cerr << "Warning: radius specified several times.\n";
                radius = std::atof(argv[i]);
                if(radius == 0) {
                    std::cerr << "Error: invalid radius.\n";
                    exit(1);
                }
            } else if(argv[i] == std::string("-m")
                   || argv[i] == std::string("--max-curvature")) {
                    if(++i == argc) {
                        std::cerr << "Error: " << argv[i-1] << " need an argument.\n";
                        usage(argv[0]);
                    }
                    if(maxCurv != 0) std::cerr << "Warning: radius specified several times.\n";
                    maxCurv = std::atof(argv[i]);
                    if(maxCurv == 0) {
                        std::cerr << "Error: invalid radius.\n";
                        exit(1);
                    }
            } else if(argv[i] == std::string("-R")
                   || argv[i] == std::string("--relative")) {
                relative = true;
            } else {
                std::cerr << "Unknown option " << argv[i] << ".\n";
                usage(argv[0]);
            }
        } else {
            if(!inFilename) {
                inFilename = argv[i];
            } else if(!outFilename) {
                outFilename = argv[i];
            } else {
                std::cerr << "Error: too much parameters.\n";
                usage(argv[0]);
            }
        }
    }

    if(!inFilename) {
        std::cerr << "Missing input file.\n";
        usage(argv[0]);
    }
    if(!outFilename) {
        std::cerr << "Missing output file.\n";
        usage(argv[0]);
    }

    Mesh mesh;
    {
        std::ifstream in(inFilename);
        PatateCommon::OBJReader<Mesh> reader;
        reader.setErrorCallback(PatateCommon::defaultErrorCallback,
                                noopErrorCallback, NULL);
        reader.read(in, mesh);
    }

    mesh.triangulate();
    mesh.computeNormals();

    if(relative) {
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

        radius *= meshRadius;
    }

    std::cout << "Radius: " << radius << "\n";

#ifdef __linux
    timespec startTime;
    int res = clock_gettime(CLOCK_MONOTONIC_RAW, &startTime);
    if(res != 0) { startTime.tv_sec = 0; startTime.tv_nsec = 0; }
#endif

    Scalar gCurvMax = fit(mesh, radius);

#ifdef __linux
    timespec endTime;
    res = clock_gettime(CLOCK_MONOTONIC_RAW, &endTime);
    if(res != 0) { endTime.tv_sec = 0; endTime.tv_nsec = 0; }
    std::cout << "Curvature computation time: " <<
                 double(double(endTime.tv_sec) - double(startTime.tv_sec))
                 + double(double(endTime.tv_nsec) - double(startTime.tv_nsec)) * 1.e-9 << " s\n";
#endif

    std::cout << "Maximum absolute computed curvature: " << gCurvMax << "\n";

    if(maxCurv == 0) {
        maxCurv = gCurvMax;
    }

    int ret = 0;
    {
        std::ofstream out(outFilename);
//        PatateCommon::OBJWriter<Mesh>().write(out, mesh);

        mesh.garbageCollection();
        for(Mesh::VertexIterator vit = mesh.verticesBegin();
            vit != mesh.verticesEnd(); ++vit)
        {
            Scalar curv = mesh.minCurvature(*vit) * mesh.maxCurvature(*vit);
//            Vector color = getColor(curv, -maxCurv, maxCurv);
            Vector color = getColor2(curv / maxCurv, .5);
            out << "v " << mesh.position(*vit).transpose()
                << " " << color.transpose() << "\n";
        }
        for(Mesh::VertexIterator vit = mesh.verticesBegin();
            vit != mesh.verticesEnd(); ++vit)
        {
            if(std::isnan(mesh.maxCurvature(*vit)) || std::isnan(mesh.minCurvature(*vit)))
                ret = 1;
            out << "vt " << mesh.maxCurvature(*vit) << " " << mesh.minCurvature(*vit) << "\n";
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
                    << "/" << (*vit).idx() + 1
                    << "/" << (*vit).idx() + 1;
                ++vit;
            } while(vit != vEnd);
            out << "\n";
        }
    }
    return ret;
}
