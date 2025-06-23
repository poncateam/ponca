
namespace Ponca {
    template<class DataPoint, class _WFunctor, int DiffType, typename T>
    void
    CurvatureEstimatorBase<DataPoint, _WFunctor, DiffType, T>::init() {
        Base::init();
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
        if(kmin <= kmax) {
            m_kmin = kmin;
            m_kmax = kmax;
            m_vmin = vmin;
            m_vmax = vmax;
        } else {
            m_kmin = kmax;
            m_kmax = kmin;
            m_vmin = vmax;
            m_vmax = vmin;
        }
        m_isValid = true;
    }
}