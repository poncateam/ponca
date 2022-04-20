/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2015 Gael Guennebaud <gael.guennebaud@inria.fr>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template<class DataPoint, class _WFunctor, typename T>
void
MeanPosition<DataPoint, _WFunctor, T>::init(const VectorType &_evalPos) {
    // Setup primitive
    Base::resetPrimitive();
    Base::basisCenter() = _evalPos;

    // Setup fitting internal values
    m_sumW = Scalar(0.0);
    m_sumP = VectorType::Zero();
}

template<class DataPoint, class _WFunctor, typename T>
bool
MeanPosition<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint &_nei) {
    VectorType q = _nei.pos() - Base::basisCenter();
    // compute weight
    Scalar w = m_w.w(q, _nei);

    if (w > Scalar(0.)) {
        m_sumP += w * q;
        m_sumW += w;

        ++(Base::m_nbNeighbors);
        return true;
    }
    return false;
}


template<class DataPoint, class _WFunctor, typename T>
FIT_RESULT
MeanPosition<DataPoint, _WFunctor, T>::finalize() {
    // handle specific configurations
    // We need to have at least one neighbor to compute the mean
    if (m_sumW == Scalar(0.) || Base::m_nbNeighbors < 1) {
        Base::resetPrimitive();
        return Base::m_eCurrentState = UNDEFINED;
    }

    return Base::m_eCurrentState = STABLE;
}

template < class DataPoint, class _WFunctor, typename T>
void
MeanNormal<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);

    // Setup fitting internal values
    m_sumN = VectorType::Zero();
}

template < class DataPoint, class _WFunctor, typename T>
bool
MeanNormal<DataPoint, _WFunctor, T>::addNeighbor(const DataPoint& _nei)
{
    if (Base::addNeighbor(_nei))
    {
        /// \todo Avoid to compute the weight multiple times
        // centered basis
        VectorType q = _nei.pos() - Base::basisCenter();

        // compute weight
        Scalar w = Base::m_w.w(q, _nei);
        m_sumN += _nei.normal() * w;

        return true;
    }
    return false;
}