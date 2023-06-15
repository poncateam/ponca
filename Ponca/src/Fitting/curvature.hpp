#include PONCA_MULTIARCH_INCLUDE_STD(cmath)   //abs

namespace Ponca {
    template<class DataPoint, class _WFunctor, int DiffType, typename T>
    void
    CurvatureEstimatorBase<DataPoint, _WFunctor, DiffType, T>::init(const VectorType &_evalPos) {
        Base::init(_evalPos);
        m_kmin = 0;
        m_kmax = 0;
        m_vmin = VectorType::Zero();
        m_vmax = VectorType::Zero();
        m_isValid = false;
    }


    template<class DataPoint, class _WFunctor, int DiffType, typename T>
    void
    CurvatureEstimatorBase<DataPoint, _WFunctor, DiffType, T>::setCurvatureValues(
            Scalar kmin, Scalar kmax, const VectorType &vmin, const VectorType &vmax) {
        PONCA_MULTIARCH_STD_MATH(abs)
        PONCA_DEBUG_ASSERT(kmin <= kmax);
        PONCA_DEBUG_ASSERT(abs(vmin.dot(vmax)) < 0.001);
        m_kmin = kmin;
        m_kmax = kmax;
        m_vmin = vmin;
        m_vmax = vmax;
        m_isValid = true;
    }
}