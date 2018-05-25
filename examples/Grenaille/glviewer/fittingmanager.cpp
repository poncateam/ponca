#include "fittingmanager.h"

FittingManager::FittingManager(QObject *parent) :
    QObject(parent),
    _fitType(FittingManager::PLANE_COV),
    _mesh(nullptr),
    _neiApproximation(new Mesh()),
    _scale(0.02),
    _evalPos(Mesh::Vector::Zero()),
    _stateUpdateNei(true)
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
    _stateUpdateNei = true;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    emit evaluationPointChanged();

    fitPrimitive();
}

void
FittingManager::setScale(float scale){
    _scale = scale;
    _stateUpdateNei = true;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    emit scaleChanged();

    fitPrimitive();
}

void
FittingManager::setScale(double scale){
    _scale = scale;
    _stateUpdateNei = true;
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    emit scaleChanged();

    fitPrimitive();
}

void
FittingManager::fitPrimitive(){

    // fit the primitive and project the neighborhood on the primitive.
    // In the current implementation (dealing only with meshes), we want to
    // keep the faces. Since we currently have no way to get the faces from the
    // points, we loop a second time over the mesh to do the reprojection.
    // This is really ugly, but this was the easy way to get good rendering for
    // meshes, without too much developement efforts.
    //
    // Workflow is as follow:
    //   - extract nei mesh
    //   - fit primitive
    //   - project nei vertices on primitive
    //

    if (_stateUpdateNei){
        Mesh nei;
        // accessing using UNSUPPORTED returns default params
        typedef typename fittingmanagerspace::BasketMaker<UNSUPPORTED>::WeightFunc WFunc;
        typedef typename Mesh::Vector VectorType;
        WFunc wfunc (_scale);

        int vertexId = 0;

        // loop over the faces, and push all the faces that have a vertex in the nei
        // like for descriptor computation, would deserve acceleration structure
        for (Mesh::faceIterator it = _mesh->faceBegin();
             it != _mesh->faceEnd();
             ++it){
            Mesh::faceIterator::value_type list = *it;

            bool add = false;
            for (int i = 0; i != 3 && !add; ++i){
                VectorType q = list[i].pos() - _evalPos;
                Scalar w = wfunc.w(q, list[i]);
                add = w > Scalar(0.);
            }


            if (add){
                nei.addVertex(list[0].pos());
                nei.addVertex(list[1].pos());
                nei.addVertex(list[2].pos());

                nei.addFace(vertexId, vertexId+1, vertexId+2);

                vertexId += 3;
            }
        }
        _neiApproximation->setFrom(nei);
    }



#define COMPUTE(T)                                                             \
    typedef typename fittingmanagerspace::BasketMaker<T>::Basket Fit;          \
    typedef typename fittingmanagerspace::BasketMaker<T>::WeightFunc WFunc;    \
    Fit f;                                                                     \
    f.setWeightFunc(WFunc(_scale));                                            \
    f.init(_evalPos);                                                          \
    f.compute(_mesh->begin(), _mesh->end());                                   \
    for (Mesh::posIterator it = _neiApproximation->vertexBegin();              \
                           it != _neiApproximation->vertexEnd(); ++it)         \
        *it = f.project(*it);


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
    case SPHERE:
    {
        COMPUTE(SPHERE);
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

/*!
 * \brief Generate a new mesh that is an approximation of the current nei.
 *
 * Internally, the weighting function is created, and used to select the vertices
 * having an influence on the fit
 */
FittingManager::Mesh *
FittingManager::getNeighborhoodMeshApprox() {
    return _neiApproximation;
}
