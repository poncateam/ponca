#ifndef FITTINGMANAGER_H
#define FITTINGMANAGER_H

#include "Patate/grenaille.h"
#include "Patate/common/gl_utils/glmesh.h"

#include <QObject>

class FittingManager : public QObject
{
    Q_OBJECT
public:
    explicit FittingManager(QObject *parent = 0);

    //! Type of supported fit types, used to generate Baskets.
    //! \warning Must be aligned with UI
    enum FIT_TYPE {
        PLANE_MEAN,
        PLANE_COV,
        MONGE_PATCH,
        SPHERE_ORIENTED,
        SPHERE,
        UNSUPPORTED
    };

    typedef PatateCommon::GLTri3DMesh Mesh;
    typedef typename Mesh::Scalar Scalar;

    inline void setMesh(Mesh* mesh) { _mesh = mesh; }

    //! \brief Getter on the projected mesh.
    Mesh *getNeighborhoodMeshApprox();

signals:
    //! \brief Signal emitted when the basket has been applied to the input data.
    void fitPerformed();
    void scaleChanged();
    void evaluationPointChanged();

public slots:
    void setBasketType(FIT_TYPE type);
    void setEvaluationPoint(const PatateCommon::GLTri3DMesh::Vector &pos);
    void setScale(float scale);
    void setScale(double scale);
    void fitPrimitive();

private:
    FIT_TYPE _fitType;
    PatateCommon::GLTri3DMesh *_mesh;
    PatateCommon::GLTri3DMesh *_neiApproximation;
    Scalar _scale;
    typename Mesh::Vector _evalPos;

    bool _stateUpdateNei;
};

namespace fittingmanagerspace {
typedef PatateCommon::GLTri3DMesh::GrenaillePoint MyPoint;

template <FittingManager::FIT_TYPE type>
struct BasketMaker {
    typedef Ponca::DistWeightFunc<MyPoint,Ponca::SmoothWeightKernel<MyPoint::Scalar> > WeightFunc;
};

template <>
struct BasketMaker<FittingManager::PLANE_MEAN>{
    typedef Ponca::DistWeightFunc<MyPoint,Ponca::SmoothWeightKernel<MyPoint::Scalar> > WeightFunc;
    typedef Ponca::Basket<MyPoint,WeightFunc, Ponca::MeanPlaneFit> Basket;
};

template <>
struct BasketMaker<FittingManager::PLANE_COV>{
    typedef Ponca::DistWeightFunc<MyPoint,Ponca::SmoothWeightKernel<MyPoint::Scalar> > WeightFunc;
    typedef Ponca::Basket<MyPoint,WeightFunc, Ponca::CovariancePlaneFit> Basket;
};

template <>
struct BasketMaker<FittingManager::MONGE_PATCH>{
    typedef Ponca::DistWeightFunc<MyPoint,Ponca::SmoothWeightKernel<MyPoint::Scalar> > WeightFunc;
    typedef Ponca::Basket<MyPoint,WeightFunc, Ponca::CovariancePlaneFit, Ponca::MongePatch> Basket;
};

template <>
struct BasketMaker<FittingManager::SPHERE_ORIENTED>{
    typedef Ponca::DistWeightFunc<MyPoint,Ponca::SmoothWeightKernel<MyPoint::Scalar> > WeightFunc;
    typedef Ponca::Basket<MyPoint,WeightFunc, Ponca::OrientedSphereFit,
                                                  Ponca::GLSParam> Basket;
};

template <>
struct BasketMaker<FittingManager::SPHERE>{
    typedef Ponca::DistWeightFunc<MyPoint,Ponca::SmoothWeightKernel<MyPoint::Scalar> > WeightFunc;
    typedef Ponca::Basket<MyPoint,WeightFunc, Ponca::SphereFit,
                                                  Ponca::GLSParam> Basket;
};

} // namespace fittingmanagerspace

#endif // FITTINGMANAGER_H