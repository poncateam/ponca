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
        PLANE_COV,
        SPHERE_ORIENTED,
        UNSUPPORTED
    };

    inline void setMesh(PatateCommon::GLTri3DMesh* mesh) { _mesh = mesh; }

signals:

public slots:
    void setBasketType(FIT_TYPE type);

private:
    FIT_TYPE _fitType;
    PatateCommon::GLTri3DMesh *_mesh;
};

namespace fittingmanagerspace {
typedef PatateCommon::GLTri3DMesh::GrenaillePoint MyPoint;

template <FittingManager::FIT_TYPE type>
struct BasketMaker { };

template <>
struct BasketMaker<FittingManager::PLANE_COV>{
    typedef Grenaille::DistWeightFunc<MyPoint,Grenaille::SmoothWeightKernel<MyPoint::Scalar> > WeightFunc;
    typedef Grenaille::Basket<MyPoint,WeightFunc, Grenaille::CompactPlane,
                                                  Grenaille::CovariancePlaneFit> Basket;
};

template <>
struct BasketMaker<FittingManager::SPHERE_ORIENTED>{
    typedef Grenaille::DistWeightFunc<MyPoint,Grenaille::SmoothWeightKernel<MyPoint::Scalar> > WeightFunc;
    typedef Grenaille::Basket<MyPoint,WeightFunc, Grenaille::AlgebraicSphere,
                                                  Grenaille::OrientedSphereFit> Basket;
};

} // namespace fittingmanagerspace

#endif // FITTINGMANAGER_H
