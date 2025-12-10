
namespace Ponca {
    ///////// CurvatureEstimatorBase
    template<class DataPoint, class _NFilter, typename T>
    void
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


    ///////// FundamentalFormCurvatureEstimator
    template<class DataPoint, class _NFilter, typename T>
    typename FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::Matrix2
    FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::firstFundamentalForm() const {
        Matrix2 first;
        firstFundamentalForm(first);
        return first;
    }

    template<class DataPoint, class _NFilter, typename T>
    template<typename Matrix2Derived>
    void FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::firstFundamentalForm(Matrix2Derived& first) const {
        Base::firstFundamentalFormComponents(first(0,0), first(1,0), first(1,1));
        first(0,1) = first(1,0); // diagonal
    }

    template<class DataPoint, class _NFilter, typename T>
    typename FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::Matrix2
    FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::secondFundamentalForm() const {
        Matrix2 second;
        secondFundamentalForm(second);
        return second;
    }

    template<class DataPoint, class _NFilter, typename T>
    template<typename Matrix2Derived>
    void FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::secondFundamentalForm(Matrix2Derived &second) const {
        Base::secondFundamentalFormComponents(second(0,0), second(1,0), second(1,1));
        second(0,1) = second(1,0); // diagonal
    }


    template<class DataPoint, class _NFilter, typename T>
    typename FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::Matrix2
    FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::weingartenMap() const {
        Matrix2 w;
        weingartenMap(w);
        return w;
    }

    template<class DataPoint, class _NFilter, typename T>
    template<typename Matrix2Derived>
    void FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::weingartenMap(Matrix2Derived &w) const {
        w = firstFundamentalForm().inverse() *  secondFundamentalForm();
    }

    template<class DataPoint, class _NFilter, typename T>
    typename FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::Scalar
    FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::kMean() const {
        Scalar E, F, G, L, M, N;
        Base::firstFundamentalFormComponents(E, F, G);
        Base::secondFundamentalFormComponents(L, M, N);
        return (G*L-Scalar(2)*F*M+E*N)/(Scalar(2)*(E*G-F*F));
    }

    template<class DataPoint, class _NFilter, typename T>
    typename FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::Scalar
    FundamentalFormWeingartenEstimator<DataPoint, _NFilter, T>::GaussianCurvature() const {
        Scalar E, F, G, L, M, N;
        Base::firstFundamentalFormComponents(E, F, G);
        Base::secondFundamentalFormComponents(L, M, N);
        return (L*N-M*M)/(E*G-F*F);
    }

    ///////// WeingartenCurvatureEstimator
    template<class DataPoint, class _NFilter, typename T>
    typename WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>::Scalar
    WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>::kMean() const {
        return Scalar(0.5)*Base::weingartenMap().trace();
    }

    template<class DataPoint, class _NFilter, typename T>
    typename WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>::Scalar
    WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>::GaussianCurvature() const {
        return Base::weingartenMap().determinant();
    }

    template<class DataPoint, class _NFilter, typename T>
    void
    WeingartenCurvatureEstimatorBase<DataPoint, _NFilter, T>::computeCurvature() const {
        if (! m_computedCurvature ){
            m_computedCurvature = true;

            Matrix2 w;
            Base::weingartenMap(w);

            // w is self adjoint by construction
            Eigen::SelfAdjointEigenSolver<Matrix2> solver(w);

            m_kmin = solver.eigenvalues().x();
            m_kmax = solver.eigenvalues().y();
            m_vmin(0) = Scalar(0); // set height
            m_vmax(0) = Scalar(0); // set height
            m_vmin.template bottomRows<2>() = solver.eigenvectors().col(0);
            m_vmax.template bottomRows<2>() = solver.eigenvectors().col(1);

            m_vmin = Base::tangentPlaneToWorld(m_vmin, false);
            m_vmax = Base::tangentPlaneToWorld(m_vmax, false);
        }
    }
}