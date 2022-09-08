#include PONCA_MULTIARCH_INCLUDE_STD(cmath)   //abs

namespace Ponca {
    template<class DataPoint, class _WFunctor, int DiffType, typename T>
    void
    CurvatureEstimatorBase<DataPoint, _WFunctor, DiffType, T>::init(const VectorType &_evalPos) {
        Base::init(_evalPos);
        m_k1 = 0;
        m_k2 = 0;
        m_v1 = VectorType::Zero();
        m_v2 = VectorType::Zero();
        m_isValid = false;
    }


    template<class DataPoint, class _WFunctor, int DiffType, typename T>
    void
    CurvatureEstimatorBase<DataPoint, _WFunctor, DiffType, T>::setCurvatureValues(
            Scalar k1, Scalar k2, const VectorType &v1, const VectorType &v2) {
        PONCA_MULTIARCH_STD_MATH(abs)

        if (abs(k1) < abs(k2)) {
            m_k1 = k2;
            m_k2 = k1;
            m_v1 = v2;
            m_v2 = v1;
        } else {
            m_k1 = k1;
            m_k2 = k2;
            m_v1 = v1;
            m_v2 = v2;
        }

        m_isValid = true;
    }
}