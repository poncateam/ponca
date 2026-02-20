
namespace Ponca {
    namespace internal
    {
        ///////// CurvatureEstimatorBase
        template<class DataPoint, class _NFilter, typename T>
        void
#if defined(__GNUC__) and not defined(__clang__)
        // attribute 'no-tree-vectorize' is required for the init function.
        // It resolves a crash when compiling with GCC 11.4.0
        // with the default optimization level (-O3) for release builds.
        __attribute__((optimize("no-tree-vectorize")))
#endif
        CurvatureEstimatorBase<DataPoint, _NFilter, T>::init() {
            Base::init();
            m_kmin = 0;
            m_kmax = 0;
            m_vmin = VectorType::Zero();
            m_vmax = VectorType::Zero();
            m_isValid = false;
        }


        template<class DataPoint, class _NFilter, typename T>
        void
        CurvatureEstimatorBase<DataPoint, _NFilter, T>::setCurvatureValues(
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
}