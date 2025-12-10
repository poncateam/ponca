
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

    ///////// NormalDerivativeWeingartenEstimator
    template<class DataPoint, class _NFilter, int DiffType, typename T>
    typename NormalDerivativeWeingartenEstimator<DataPoint, _NFilter, DiffType, T>::Matrix2
    NormalDerivativeWeingartenEstimator<DataPoint, _NFilter, DiffType, T>::weingartenMap() const {
        Matrix2 w;
        weingartenMap(w);
        return w;
    }

    template<class DataPoint, class _NFilter, int DiffType, typename T>
    template<typename Matrix2Derived>
    void NormalDerivativeWeingartenEstimator<DataPoint, _NFilter, DiffType, T>::weingartenMap(Matrix2Derived &W) const {
        PONCA_MULTIARCH_STD_MATH(abs);
        PONCA_MULTIARCH_STD_MATH(sqrt);

        using Index = typename VectorType::Index;
        using Matrix32 = Eigen::Matrix<Scalar,3,2>;

        // Get the object space Weingarten map dN
        MatrixType dN = Base::dNormal().template middleCols<DataPoint::Dim>(Base::isScaleDer() ? 1: 0);

        // Compute tangent-space basis
        // \FIXME This code provides tangent plane estimation, but does not make it available to other classes
        //        Could be interesting to make it available through PROVIDES_TANGENT_PLANE_BASIS
        Matrix32 B;
        {
            Index i0=Index(-1), i1=Index(-1), i2=Index(-1);
            // Use the spatial derivative of the normal. This option leads to NaN
            // values if dN is null (like in the case of a perfect plane)

            // Get the object space Weingarten map dN
            MatrixType dN = Base::dNormal().template middleCols<DataPoint::Dim>(Base::isScaleDer() ? 1 : 0);

            // Compute tangent-space basis from dN
            //   1 - pick the column with maximal norm as the first tangent vector,
            Scalar sqNorm = dN.colwise().squaredNorm().maxCoeff(&i0);
            B.col(0) = dN.col(i0) / sqrt(sqNorm);
            //   2 - orthogonalize the other column vectors, and pick the most reliable one
            i1 = (i0 + 1) % 3;
            i2 = (i0 + 2) % 3;
            VectorType v1 = dN.col(i1) - B.col(0).dot(dN.col(i1)) * B.col(0);
            VectorType v2 = dN.col(i2) - B.col(0).dot(dN.col(i2)) * B.col(0);
            Scalar v1norm2 = v1.squaredNorm();
            Scalar v2norm2 = v2.squaredNorm();
            if (v1norm2 > v2norm2) B.col(1) = v1 / sqrt(v1norm2);
            else B.col(1) = v2 / sqrt(v2norm2);
        }

        // Compute the 2x2 matrix representing the shape operator by transforming dN to the basis B.
        // Recall that dN is a bilinear form, it thus transforms as follows:
        W = B.transpose() * dN * B;

        // Recall that at this stage, the shape operator represented by W describes the normal curvature K_n(u) in the direction u \in R^2 as follows:
        //   K_n(u) = u^T W u
        // The principal curvatures are fully defined by the values and the directions of the extrema of K_n.
        //
        // If the normal field N(x) comes from the gradient of a scalar field, then N(x) is curl-free, and dN and W are symmetric matrices.
        // In this case, the extrema of the previous quadratic form are directly obtained by the eigenvalue decomposition of W.
        // However, if N(x) is only an approximation of the normal field of a surface, then N(x) is not necessarily curl-free, and in this case W is not symmetric.
        // In this case, we first have to find an equivalent symmetric matrix W' such that:
        //   K_n(u) = u^T W' u,
        // for any u \in R^2.
        // It is easy to see that such a W' is simply obtained as:
        //   W' = (W + W^T)/2
        W(0,1) = W(1,0) = (W(0,1) + W(1,0))/Scalar(2);
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