template < class DataPoint, class _WFunctor, int DiffType, typename T>
FIT_RESULT
NormalDerivativesCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::finalize()
{
    if (Base::finalize() == STABLE) {
        if (Base::curvatureEstimatorBase().isValid()) Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
        Base::m_eCurrentState = computeCurvature(false);
    }

    return Base::m_eCurrentState;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
FIT_RESULT NormalDerivativesCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::computeCurvature(bool useNormal)
{
    PONCA_MULTIARCH_STD_MATH(abs);

    // Get the object space Weingarten map dN
    MatrixType dN = Base::dNormal().template middleCols<DataPoint::Dim>(Base::isScaleDer() ? 1: 0);

    // Compute tangent-space basis
    Mat32 B = tangentPlane(useNormal);

    // Compute the 2x2 matrix representing the shape operator by transforming dN to the basis B.
    // Recall that dN is a bilinear form, it thus transforms as follows:
    Mat22 S = B.transpose() * dN * B;

    // Recall that at this stage, the shape operator represented by S describes the normal curvature K_n(u) in the direction u \in R^2 as follows:
    //   K_n(u) = u^T S u
    // The principal curvatures are fully defined by the values and the directions of the extrema of K_n.
    //
    // If the normal field N(x) comes from the gradient of a scalar field, then N(x) is curl-free, and dN and S are symmetric matrices.
    // In this case, the extrema of the previous quadratic form are directly obtained by the the eigenvalue decomposition of S.
    // However, if N(x) is only an approximation of the normal field of a surface, then N(x) is not necessarily curl-free, and in this case S is not symmetric.
    // In this case, we first have to find an equivalent symmetric matrix S' such that:
    //   K_n(u) = u^T S' u,
    // for any u \in R^2.
    // It is easy to see that such a S' is simply obtained as:
    //   S' = (S + S^T)/2

    S(0,1) = S(1,0) = (S(0,1) + S(1,0))/Scalar(2);
    Eigen::SelfAdjointEigenSolver<Mat22> eig2;
    eig2.computeDirect(S);

    if (eig2.info() != Eigen::Success) return UNDEFINED;

    Base::setCurvatureValues(eig2.eigenvalues()(0), eig2.eigenvalues()(1),
                             B * eig2.eigenvectors().col(0), B * eig2.eigenvectors().col(1));

    return Base::m_eCurrentState;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
typename NormalDerivativesCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::Mat32
NormalDerivativesCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::tangentPlane(bool useNormal) const
{
    typedef typename VectorType::Index Index;

    PONCA_MULTIARCH_STD_MATH(sqrt);
    PONCA_MULTIARCH_STD_MATH(abs);

    Mat32 B;
    Index i0=Index(-1), i1=Index(-1), i2=Index(-1);

    // Two choices to compute a basis of the tangent plane
    if(useNormal)
    {
        // Use a vector orthogonal to the surface (the gradient) and compute an
        // orthonormal basis from it
        VectorType n = Base::primitiveGradient();
        n.array().abs().minCoeff(&i0); // i0: dimension where n extends the least
        i1 = (i0+1)%3;
        i2 = (i0+2)%3;

        B.col(0)[i0] = 0;
        B.col(0)[i1] = n[i2];
        B.col(0)[i2] = -n[i1];

        B.col(0).normalize();
        B.col(1) = B.col(0).cross(n);
    }
    else
    {
        // Use the spatial derivative of the normal. This option leads to NaN
        // values if dN is null (like in the case of a perfect plane)

        // Get the object space Weingarten map dN
        MatrixType dN = Base::dNormal().template middleCols<DataPoint::Dim>(Base::isScaleDer() ? 1: 0);

        // Compute tangent-space basis from dN
        //   1 - pick the column with maximal norm as the first tangent vector,
        Scalar sqNorm = dN.colwise().squaredNorm().maxCoeff(&i0);
        B.col(0) = dN.col(i0) / sqrt(sqNorm);
        //   2 - orthogonalize the other column vectors, and pick the most reliable one
        i1 = (i0+1)%3;
        i2 = (i0+2)%3;
        VectorType v1 = dN.col(i1) - B.col(0).dot(dN.col(i1)) * B.col(0);
        VectorType v2 = dN.col(i2) - B.col(0).dot(dN.col(i2)) * B.col(0);
        Scalar v1norm2 = v1.squaredNorm();
        Scalar v2norm2 = v2.squaredNorm();
        if(v1norm2 > v2norm2) B.col(1) = v1 / sqrt(v1norm2);
        else                  B.col(1) = v2 / sqrt(v2norm2);
    }

    return B;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
void
NormalCovarianceCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::init()
{
    Base::init();
    m_cov = MatrixType::Zero();
    m_cog = VectorType::Zero();
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
bool
NormalCovarianceCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::addLocalNeighbor
        (Scalar w, const VectorType &localQ, const DataPoint &attributes, ScalarArray &dw)
{
    if( Base::addLocalNeighbor(w, localQ, attributes, dw) )
    {
        m_cov += w * attributes.normal() * attributes.normal().transpose();
        m_cog += attributes.normal();
        return true;
    }
    return false;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
FIT_RESULT
NormalCovarianceCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::finalize ()
{
    typedef typename VectorType::Index Index;

    PONCA_MULTIARCH_STD_MATH(abs);

    FIT_RESULT res = Base::finalize();

    if(this->isReady())
    {
        // center of gravity (mean)
        m_cog /= Base::getWeightSum();

        // Center the covariance on the centroid
        m_cov = m_cov/Base::getWeightSum() - m_cog * m_cog.transpose();

        m_solver.computeDirect(m_cov);

        Scalar kmin = m_solver.eigenvalues()(1);
        Scalar kmax = m_solver.eigenvalues()(2);

        VectorType vmin = m_solver.eigenvectors().col(1);
        VectorType vmax = m_solver.eigenvectors().col(2);

        // fallback
        // TODO(thib) which epsilon value should be chosen ?
        Scalar epsilon = Scalar(1e-3);
        if(kmin<epsilon && kmax<epsilon)
        {
            kmin = Scalar(0);
            kmax = Scalar(0);

            // set principal directions from normals center of gravity
            VectorType n = m_cog.normalized();
            Index i0 = -1, i1 = -1, i2 = -1;
            n.rowwise().squaredNorm().minCoeff(&i0);
            i1 = (i0+1)%3;
            i2 = (i0+2)%3;
            vmin[i0] = 0;
            vmin[i1] = n[i2];
            vmin[i2] = -n[i1];
            vmax[i0] = n[i1]*n[i1] + n[i2]*n[i2];
            vmax[i1] = -n[i1]*n[i0];
            vmax[i2] = -n[i2]*n[i0];
        }

        Base::setCurvatureValues(kmin, kmax, vmin, vmax);
    }
    return res;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
void
ProjectedNormalCovarianceCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::init()
{
    Base::init();

    m_cog = Vector2::Zero();
    m_cov = Mat22::Zero();
    m_pass = FIRST_PASS;
    m_tframe = Mat32::Zero();
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
bool
ProjectedNormalCovarianceCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::addLocalNeighbor
        (Scalar w, const VectorType &localQ, const DataPoint &attributes, ScalarArray &dw)
{
    if(m_pass == FIRST_PASS)
    {
        return Base::addLocalNeighbor(w, localQ, attributes, dw);
    }
    else if(m_pass == SECOND_PASS)
    {
        // project normal on plane
        VectorType n = attributes.normal();
        Vector2 proj = m_tframe.transpose() * n;

        m_cov += proj * proj.transpose();
        m_cog += w * proj;
        return true;
    }
    return false;
}

template < class DataPoint, class _WFunctor, int DiffType, typename T>
FIT_RESULT
ProjectedNormalCovarianceCurvatureEstimator<DataPoint, _WFunctor, DiffType, T>::finalize ()
{
    if(m_pass == FIRST_PASS)
    {
        FIT_RESULT res = Base::finalize();

        if(res != UNDEFINED)
        {
            // get normal of fitted plane
            VectorType n = this->primitiveGradient(VectorType());

            // compute orthonormal frame of the tangent plane
            Index i0 = -1, i1 = -1, i2 = -1;
            n.rowwise().squaredNorm().minCoeff(&i0);
            i1 = (i0+1)%3;
            i2 = (i0+2)%3;
            m_tframe.col(0)[i0] = 0;
            m_tframe.col(0)[i1] = n[i2];
            m_tframe.col(0)[i2] = -n[i1];
            m_tframe.col(1)[i0] = n[i1]*n[i1] + n[i2]*n[i2];
            m_tframe.col(1)[i1] = -n[i1]*n[i0];
            m_tframe.col(1)[i2] = -n[i2]*n[i0];

            // go to second pass
            m_pass = SECOND_PASS;
            return NEED_OTHER_PASS;
        }
        else
        {
            return UNDEFINED;
        }
    }
    else if(m_pass == SECOND_PASS)
    {
        // center of gravity (mean)
        m_cog /= Base::getWeightSum();

        // Center the covariance on the centroid
        m_cov = m_cov/Base::getWeightSum() - m_cog * m_cog.transpose();

        m_solver.computeDirect(m_cov);

        Base::m_kmin = m_solver.eigenvalues()(0);
        Base::m_kmax = m_solver.eigenvalues()(1);

        // transform from local plane coordinates to world coordinates
        Base::m_v1 = m_tframe * m_solver.eigenvectors().col(0);
        Base::m_v2 = m_tframe * m_solver.eigenvectors().col(1);

        //TODO(thib) which epsilon value should be chosen ?
//        Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();
        Scalar epsilon = Scalar(1e-3);
        if(Base::m_kmin<epsilon && Base::m_kmax<epsilon)
        {
            Base::m_kmin = Scalar(0);
            Base::m_kmax = Scalar(0);

            // set principal directions from fitted plane
            Base::m_v2 = m_tframe.col(0);
            Base::m_v1 = m_tframe.col(1);
        }

        return STABLE;
    }
    return UNDEFINED;
}
