/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
 GLSDer <DataPoint, _WFunctor, T>::dtau() const
{
    MULTIARCH_STD_MATH(sqrt);

    Scalar prattNorm2 = Base::prattNorm2();
    Scalar prattNorm  = sqrt(prattNorm2);
    Scalar cfactor    = Scalar(.5) / prattNorm;

    return (Base::m_dUc * prattNorm - Base::m_uc * cfactor * Base::dprattNorm2()) / prattNorm2;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::VectorArray
GLSDer <DataPoint, _WFunctor, T>::deta() const
{
    Scalar ulN  = Base::m_ul.norm();
    Scalar ulN3 = ulN * ulN * ulN;
    return Base::m_dUl * (Scalar(1.) / ulN) - Base::m_ul * Base::m_ul.transpose() * Base::m_dUl * (Scalar(1.) / ulN3);
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dkappa() const
{
    MULTIARCH_STD_MATH(sqrt);

    Scalar prattNorm2 = Base::prattNorm2();
    Scalar prattNorm  = sqrt(prattNorm2);
    Scalar cfactor    = Scalar(.5) / prattNorm;

    return Scalar(2.) * (Base::m_dUq * prattNorm - Base::m_uq * cfactor * Base::dprattNorm2()) / prattNorm2;
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dtau_normalized() const
{
    return dtau();
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::VectorArray
GLSDer <DataPoint, _WFunctor, T>::deta_normalized() const
{
    return Base::m_t * deta();
}


template < class DataPoint, class _WFunctor, typename T>
typename GLSDer <DataPoint, _WFunctor, T>::ScalarArray
GLSDer <DataPoint, _WFunctor, T>::dkappa_normalized() const
{
    return dkappa() * Base::m_t * Base::m_t;
}




template < class DataPoint, class _WFunctor, typename T>
typename GLSGeomVar <DataPoint, _WFunctor, T>::Scalar
GLSGeomVar <DataPoint, _WFunctor, T>::geomVar(  Scalar wtau, 
                                                Scalar weta,
                                                Scalar wkappa ) const
{
    Scalar dtau   = Base::dtau_normalized().col(0)(0);
    Scalar deta   = Base::deta_normalized().col(0).norm();
    Scalar dkappa = Base::dkappa_normalized().col(0)(0);

    return wtau*dtau*dtau + weta*deta*deta + wkappa*dkappa*dkappa;
}

template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
GLSCurvatureHelper<DataPoint, _WFunctor, T>::finalize()
{
    MULTIARCH_STD_MATH(sqrt);

    FIT_RESULT bResult = Base::finalize();

    if(bResult != UNDEFINED)
    {
        // Extract the spatial variations of eta
        MatrixType jacobian = Base::deta().template middleCols<DataPoint::Dim>(Base::isScaleDer() ? 1: 0);

        // Use a simple solver with 2x2 and 3x3 closed forms compatible with eigen-nvcc
        // This solver requires a triangular matrix
        Eigen::SelfAdjointEigenSolver<MatrixType> eig;
#ifdef __CUDACC__
        eig.computeDirect(jacobian.transpose()*jacobian);
#else
        eig.compute(jacobian.transpose()*jacobian);
#endif

        if (eig.info() != Eigen::Success){
            return UNDEFINED;
        }

        // Need to detect which vector is the most aligned to the smoothed
        // normal direction, and pick the two others as principal directions
        int idK1 = 2,
            idK2 = 1;
        {
            // compute the abs dot product between eta and the eigen vectors
            VectorType dotVec = (eig.eigenvectors().transpose() * Base::eta()).array().abs();
            int maxCoeff = dotVec(0) > dotVec(1) ?
                             (dotVec(0) > dotVec(2) ? 0 : 2) :
                             (dotVec(1) > dotVec(2) ? 1 : 2);
            
            // swap ids if needed
            if (maxCoeff != 0){
                idK2 = 0;
                if (maxCoeff == 2) idK1 = 1;
            }
        }

        // Extract eigenvectors and eigen values
        // Need sqrt because we solve eigendecomposition of JT * J.
        m_k1 = sqrt(eig.eigenvalues()(idK1)); 
        m_k2 = sqrt(eig.eigenvalues()(idK2)); 

        m_v1 = eig.eigenvectors().col(idK1);
        m_v2 = eig.eigenvectors().col(idK2);

        // Now check the sign of the mean curvature to detect if we need to change
        // the sign of the principal curvature values:
        // the eigen decomposition return k1 and k2 wrt k1*k1 > k2*k2

        // Compare with the values of the mean curvature computed without k1 and k2
        Scalar H2 = Scalar(2) * Base::kappa(); // we could also use the trace of the
        // jacobian matrix to get mean curvature

        // Change sign and flip values if needed
        // Thanks to Noam Kremen snoamk@tx.technion.ac.il for this algorithm
        if( H2 == Scalar(0))
        {
            m_k2 = -m_k2;
        }
        else if ( H2 > Scalar(0) )
        {
            if( H2 < m_k1 )
            {
                m_k2 = -m_k2;
            }
        }
        else // 2H < 0. In this case, we have k1<0, and k1 < k2
        {
            if( H2 < - m_k1 )
            {
                m_k2 = -m_k2;
            }
            m_k1 = -m_k1;

            // need to invert k1 and k2, and get the corresponding vectors
            Scalar tmp = m_k1;
            m_k1 = m_k2;
            m_k2 = tmp;        
            m_v1 = eig.eigenvectors().col(idK2);
            m_v2 = eig.eigenvectors().col(idK1);        
        }
    }

    return bResult;
}
