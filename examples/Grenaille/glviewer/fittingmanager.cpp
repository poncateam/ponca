#include "fittingmanager.h"

FittingManager::FittingManager(QObject *parent) :
    QObject(parent),
    _fitType(FittingManager::PLANE_COV),
    _mesh(NULL),
    _scale(0.02),
    _evalPos(Mesh::Vector::Zero())
{
}
void
FittingManager::setBasketType(FIT_TYPE type){
    _fitType = type;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    fitPrimitive();
}

void
FittingManager::setEvaluationPoint(const PatateCommon::GLTri3DMesh::Vector &pos){
    _evalPos = pos;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    fitPrimitive();
}

void
FittingManager::setScale(float scale){
    _scale = scale;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    fitPrimitive();
}

void
FittingManager::setScale(double scale){
    _scale = scale;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    fitPrimitive();
}

void
FittingManager::fitPrimitive(){
#define COMPUTE(T)                                                             \
    typedef typename fittingmanagerspace::BasketMaker<T>::Basket Fit;          \
    typedef typename fittingmanagerspace::BasketMaker<T>::WeightFunc WFunc;    \
    Fit f;                                                                     \
    f.setWeightFunc(WFunc(_scale));                                            \
    f.init(_evalPos);                                                          \
    f.compute(_mesh->begin(), _mesh->end());


    switch(_fitType){
    case PLANE_COV:
    {
        COMPUTE(PLANE_COV);
        std::cout << "Normal vector:     " << f.primitiveGradient(_evalPos).transpose() << std::endl;
        std::cout << "Surface variation: " << f.surfaceVariation() << std::endl;
        break;
    }
    case SPHERE_ORIENTED:
    {
        COMPUTE(SPHERE_ORIENTED);
        std::cout << "Tau:   " << f.tau() << std::endl;
        std::cout << "Eta:   " << f.eta().transpose() << std::endl;
        std::cout << "Kappa: " << f.kappa() << std::endl;
        break;
    }
    case UNSUPPORTED:
    default:
        return;
    }

    emit fitPerformed();
}
