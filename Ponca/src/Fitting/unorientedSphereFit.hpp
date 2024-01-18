/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

namespace Ponca
{

template < class DataPoint, class _WFunctor, typename T>
void
UnorientedSphereFitImpl<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);
    m_matA.setZero();
    m_matQ.setZero();
    m_sumDotPP = Scalar(0.0);
}

template<class DataPoint, class _WFunctor, typename T>
bool
UnorientedSphereFitImpl<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
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

template < class DataPoint, class _WFunctor, typename T>
FIT_RESULT
UnorientedSphereFitImpl<DataPoint, _WFunctor, T>::finalize ()
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

} //namespace Ponca