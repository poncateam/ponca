
#include <Eigen/SVD>
#include <Eigen/Geometry>

#include "mongePatch.h"

template < class DataPoint, class _NFilter, typename T>
void
MongePatchQuadraticFitImpl<DataPoint, _NFilter, T>::init()
{
    Base::init();

    m_A = SampleMatrix(6,6);
    m_A.setZero();
    m_b.setZero();
    m_planeIsReady = false;
}

template < class DataPoint, class _NFilter, typename T>
bool
MongePatchQuadraticFitImpl<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                      const VectorType &localQ,
                                                      const DataPoint &attributes)
{
    if(! m_planeIsReady)
    {
        return Base::addLocalNeighbor(w, localQ, attributes);
    }
    else // base plane is ready, we can now fit the patch
    {
        // express neighbor in local coordinate frame
        const VectorType local = Base::worldToTangentPlane(attributes.pos());
        const Scalar& h = Base::getHFromLocalCoordinates(local);
        const Scalar& u = Base::getUFromLocalCoordinates(local);
        const Scalar& v = Base::getVFromLocalCoordinates(local);

        Eigen::Matrix<Scalar, 6, 1 > p;
        p << u*u, v*v, u*v, u, v, 1;
        m_A    += w*p*p.transpose();
        m_b    += w*h*p;

        return true;
    }
    return false;
}

template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
MongePatchQuadraticFitImpl<DataPoint, _NFilter, T>::finalize ()
{
    // end of the fitting process, check plane is ready
    if (! m_planeIsReady) {
        FIT_RESULT res = Base::finalize();

        if(res == STABLE) {  // plane is ready
            m_planeIsReady = true;

            return Base::m_eCurrentState = NEED_OTHER_PASS;
        }
        return res;
    }
    // end of the monge patch fitting process
    else {
        // we use BDCSVD as the matrix size is 36
        // http://eigen.tuxfamily.org/dox/classEigen_1_1BDCSVD.html
        Base::quadraticHeightField().setQuadric
        (m_A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV | Eigen::NoQRPreconditioner).solve(m_b) );

        return Base::m_eCurrentState = STABLE;
    }
}



template < class DataPoint, class _NFilter, typename T>
void
MongePatchRestrictedQuadraticFitImpl<DataPoint, _NFilter, T>::init()
{
    Base::init();

    m_A = SampleMatrix(3,3);
    m_A.setZero();
    m_b.setZero();
    m_planeIsReady = false;
}

template < class DataPoint, class _NFilter, typename T>
bool
MongePatchRestrictedQuadraticFitImpl<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                                     const VectorType &localQ,
                                                                     const DataPoint &attributes)
{
    if(! m_planeIsReady)
    {
        return Base::addLocalNeighbor(w, localQ, attributes);
    }
    else // base plane is ready, we can now fit the patch
    {
        // express neighbor in local coordinate frame
        const VectorType local = Base::worldToTangentPlane(attributes.pos());
        const Scalar& h = Base::getHFromLocalCoordinates(local);
        const Scalar& u = Base::getUFromLocalCoordinates(local);
        const Scalar& v = Base::getVFromLocalCoordinates(local);

        Eigen::Matrix<Scalar, 3, 1 > p;
        p << u*u, v*v, u*v;
        m_A    += w*p*p.transpose();
        m_b    += w*h*p;

        return true;
    }
    return false;
}

template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
MongePatchRestrictedQuadraticFitImpl<DataPoint, _NFilter, T>::finalize ()
{
    // end of the fitting process, check plane is ready
    if (! m_planeIsReady) {
        FIT_RESULT res = Base::finalize();

        if(res == STABLE) {  // plane is ready
            m_planeIsReady = true;

            return Base::m_eCurrentState = NEED_OTHER_PASS;
        }
        return res;
    }
        // end of the monge patch fitting process
    else {
        // we use SVD as the matrix size is 3x3
        Base::quadraticHeightField().setQuadric
        (m_A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV | Eigen::NoQRPreconditioner).solve(m_b));

        return Base::m_eCurrentState = STABLE;
    }
}

//template < class DataPoint, class _NFilter, typename T>
//typename MongePatch<DataPoint, _NFilter, T>::Scalar
//MongePatch<DataPoint, _NFilter, T>::kMean() const {
//  PONCA_MULTIARCH_STD_MATH(pow);
//  static const Scalar one (1);
//  static const Scalar two (2);
//  static const Scalar threeOverTwo (Scalar(3)/Scalar(2));
//  return ((one + pow(h_v(),two) ) * h_uu() * two*h_u()*h_v()*h_uv() + (one+pow(h_u(),two))*h_vv()) /
//      (two * pow(one +pow(h_u(),two) + pow(h_v(),two),threeOverTwo));
//}
//
//template < class DataPoint, class _NFilter, typename T>
//typename MongePatch<DataPoint, _NFilter, T>::Scalar
//MongePatch<DataPoint, _NFilter, T>::GaussianCurvature() const {
//    PONCA_MULTIARCH_STD_MATH(pow);
//    static const Scalar one (1);
//    static const Scalar two (2);
//    return (h_uu()*h_vv() - pow(h_uv(),two)) /
//        pow((one + pow(h_u(),two) + pow(h_v(),two) ), two);
//}
