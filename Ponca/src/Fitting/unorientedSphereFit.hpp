/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

namespace Ponca
{

template < class DataPoint, class _NFilter, typename T>
void
UnorientedSphereFitImpl<DataPoint, _NFilter, T>::init()
{
    Base::init();
    m_matA.setZero();
    m_matQ.setZero();
    m_sumDotPP = Scalar(0.0);
}

template<class DataPoint, class _NFilter, typename T>
bool
UnorientedSphereFitImpl<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                        const VectorType &localQ,
                                                        const DataPoint &attributes) {
    if( Base::addLocalNeighbor(w, localQ, attributes) )
    {
        VectorB basis;
        basis << attributes.normal(), attributes.normal().dot(localQ);

        m_matA     += w * basis * basis.transpose();
        m_sumDotPP += w * localQ.squaredNorm();

        return true;
    }

    return false;
}

template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
UnorientedSphereFitImpl<DataPoint, _NFilter, T>::finalize ()
{
    PONCA_MULTIARCH_STD_MATH(sqrt);
    constexpr int Dim = DataPoint::Dim;

    // Compute status
    if(Base::finalize() != STABLE)
        return Base::m_eCurrentState;
    if(Base::getNumNeighbors() < DataPoint::Dim)
        return Base::m_eCurrentState = UNDEFINED;
    if (Base::algebraicSphere().isValid())
        Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
    else
        Base::m_eCurrentState = Base::getNumNeighbors() < 2*DataPoint::Dim ? UNSTABLE : STABLE;

    // 1. finalize sphere fitting
    Scalar invSumW = Scalar(1.) / Base::getWeightSum();

    m_matQ.template topLeftCorner<Dim,Dim>().setIdentity();
    m_matQ.col(Dim).template head<Dim>() = Base::m_sumP * invSumW;
    m_matQ.row(Dim).template head<Dim>() = Base::m_sumP * invSumW;
    m_matQ(Dim,Dim) = m_sumDotPP * invSumW;

    MatrixBB M = m_matQ.inverse() * m_matA;
    m_solver.compute(M);
    VectorB eivals = m_solver.eigenvalues().real();
    int maxId = 0;
    eivals.maxCoeff(&maxId);

    VectorB eivec = m_solver.eigenvectors().col(maxId).real();

    // integrate
    Base::m_ul = eivec.template head<Dim>();
    Base::m_uq = Scalar(0.5) * eivec(Dim);
    Base::m_uc = -invSumW * (Base::m_ul.dot(Base::m_sumP) + m_sumDotPP * Base::m_uq);

    Base::m_isNormalized = false;

    return Base::m_eCurrentState;
}



template < class DataPoint, class _NFilter, int DiffType, typename T>
void
UnorientedSphereDerImpl<DataPoint, _NFilter, DiffType, T>::init()
{
    Base::init();
    for(int dim = 0; dim < Base::NbDerivatives; ++dim)
        m_dmatA[dim] = MatrixBB::Zero();
    m_dSumDotPP = ScalarArray::Zero();
    m_dUc = ScalarArray::Zero();
    m_dUl = VectorArray::Zero();
    m_dUq = ScalarArray::Zero();
}

template < class DataPoint, class _NFilter, int DiffType, typename T>
bool
UnorientedSphereDerImpl<DataPoint, _NFilter, DiffType, T>::addLocalNeighbor(
    Scalar w,
    const VectorType &localQ,
    const DataPoint &attributes,
    ScalarArray& dw)
{
    if( Base::addLocalNeighbor(w, localQ, attributes, dw) )
    {
        VectorB basis;
        basis << attributes.normal(), attributes.normal().dot(localQ);
        const MatrixBB prod = basis * basis.transpose();

        m_dSumDotPP += dw * localQ.squaredNorm();

        for(int dim = 0; dim < Base::NbDerivatives; ++dim)
        {
            m_dmatA[dim] += dw[dim] * prod;
        }

        return true;
    }
    return false;
}

template < class DataPoint, class _NFilter, int DiffType, typename T>
FIT_RESULT
UnorientedSphereDerImpl<DataPoint, _NFilter, DiffType, T>::finalize()
{
    constexpr int Dim = DataPoint::Dim;
    constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();

    Base::finalize();
    if(this->isReady())
    {
        int i = 0;
        Base::m_solver.eigenvalues().real().maxCoeff(&i);
        const Scalar eigenval_i = Base::m_solver.eigenvalues().real()[i];
        const VectorB eigenvec_i = Base::m_solver.eigenvectors().real().col(i);

        const Scalar invSumW = Scalar(1) / Base::getWeightSum();

        // the derivative of (A vi = li Q vi), where A and Q are symmetric real, is
        //
        //     dA vi + A dvi = dli Q vi + li dQ vi + li Q dvi
        //
        // left-multiply by vj (i!=j, associated to eigenvalue lj)
        //
        //     vj . dA vi + vj . A dvi = dli vj . Q vi + li vj . dQ vi + li vj . Q dvi
        //                  ^^^^^^           ^^^^^^^^^
        //                  = lj vj . Q         = 0
        //
        // which simplified to
        //
        //     (Q vj) . dvi = 1/(li - lj) vj . dQ vi
        //
        // therefore
        //
        //     dvi = sum_j (1/(li - lj) vj . dQ vi) Q vj
        //         = Q sum_j (1/(li - lj) vj . ((dA - li dQ) vi)) vj
        //
        for(int dim = 0; dim < Base::NbDerivatives; ++dim)
        {
            MatrixBB dQ;
            dQ.template topLeftCorner<Dim,Dim>().setZero();
            dQ.col(Dim).template head<Dim>() = Base::m_dSumP.col(dim);
            dQ.row(Dim).template head<Dim>() = Base::m_dSumP.col(dim);
            dQ(Dim,Dim) = m_dSumDotPP[dim];
            dQ -= Base::m_dSumW[dim] * Base::m_matQ;
            dQ *= invSumW;

            const VectorB vec = (m_dmatA[dim] - eigenval_i * dQ) * eigenvec_i;

            VectorB deigvec = VectorB::Zero();
            for(int j = 0; j < Dim+1; ++j)
            {
                if(j != i)
                {
                    const Scalar eigenval_j = Base::m_solver.eigenvalues().real()[j];
                    const Scalar eigengap = eigenval_i - eigenval_j; // positive since eigenval_i is the maximal eigenvalue

                    if(eigengap > epsilon)
                    {
                        const VectorB eigenvec_j = Base::m_solver.eigenvectors().real().col(j);
                        deigvec += Scalar(1)/eigengap * eigenvec_j.dot(vec) * eigenvec_j;
                    }
                }
            }

            deigvec = Base::m_matQ * deigvec;

            m_dUq[dim] = Scalar(.5) * deigvec[dim];
            m_dUl.col(dim) = deigvec.template head<Dim>();
        }

        // same as in OrientedSphereDerImpl::finalize()
        m_dUc = -invSumW*( Base::m_sumP.transpose() * m_dUl
            + Base::m_sumDotPP * m_dUq
            + Base::m_ul.transpose() * Base::m_dSumP
            + Base::m_uq * m_dSumDotPP
            + Base::m_dSumW * Base::m_uc);

        return FIT_RESULT::STABLE;
    }
    return Base::m_eCurrentState;
}

template < class DataPoint, class _NFilter, int DiffType, typename T>
typename UnorientedSphereDerImpl<DataPoint, _NFilter, DiffType, T>::ScalarArray
UnorientedSphereDerImpl<DataPoint, _NFilter, DiffType, T>::dPotential() const
{
    // same as OrientedSphereDerImpl::dPotential()
    ScalarArray dfield = m_dUc;
    if(Base::isSpaceDer())
        dfield.template tail<DataPoint::Dim>() += Base::m_ul;
    return dfield;
}

template < class DataPoint, class _NFilter, int DiffType, typename T>
typename UnorientedSphereDerImpl<DataPoint, _NFilter, DiffType, T>::VectorArray
UnorientedSphereDerImpl<DataPoint, _NFilter, DiffType, T>::dNormal() const
{
    // same as OrientedSphereDerImpl::dNormal()
    VectorArray dgrad = m_dUl;
    if(Base::isSpaceDer())
        dgrad.template rightCols<DataPoint::Dim>().diagonal().array() += Scalar(2)*Base::m_uq;
    Scalar norm  = Base::m_ul.norm();
    Scalar norm3 = norm*norm*norm;
    return dgrad / norm - Base::m_ul * (Base::m_ul.transpose() * dgrad) / norm3;
}

} //namespace Ponca
