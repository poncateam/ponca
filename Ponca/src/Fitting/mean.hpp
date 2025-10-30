/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
 Copyright (C) 2015 Gael Guennebaud <gael.guennebaud@inria.fr>

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

template<class DataPoint, class _NFilter, typename T>
void
MeanPosition<DataPoint, _NFilter, T>::init() {
    Base::init();
    m_sumP = VectorType::Zero();
}

template<class DataPoint, class _NFilter, typename T>
bool
MeanPosition<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                        const VectorType &localQ,
                                                        const DataPoint &attributes) {
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        m_sumP += w * localQ;
        return true;
    }
    return false;
}

template < class DataPoint, class _NFilter, typename T>
void
MeanNormal<DataPoint, _NFilter, T>::init()
{
    Base::init();
    m_sumN = VectorType::Zero();
}


template<class DataPoint, class _NFilter, typename T>
bool
MeanNormal<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                      const VectorType &localQ,
                                                      const DataPoint &attributes) {
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        m_sumN += w * attributes.normal();
        return true;
    }
    return false;
}

template<class DataPoint, class _NFilter, int DiffType, typename T>
void
MeanPositionDer<DataPoint, _NFilter, DiffType, T>::init() {
    Base::init();
    m_dSumP.setZero();
}


template<class DataPoint, class _NFilter, int DiffType, typename T>
bool
MeanPositionDer<DataPoint, _NFilter, DiffType, T>::addLocalNeighbor(Scalar w,
                                                                 const VectorType &localQ,
                                                                 const DataPoint &attributes,
                                                                 ScalarArray &dw) {
    if (Base::addLocalNeighbor(w, localQ, attributes, dw)) {
        m_dSumP += localQ * dw;
        return true;
    }

    return false;
}


template<class DataPoint, class _NFilter, int DiffType, typename T>
void
MeanNormalDer<DataPoint, _NFilter, DiffType, T>::init() {
    Base::init();
    m_dSumN.setZero();
}

template<class DataPoint, class _NFilter, int DiffType, typename T>
bool
MeanNormalDer<DataPoint, _NFilter, DiffType, T>::addLocalNeighbor(Scalar w,
                                                                 const VectorType &localQ,
                                                                 const DataPoint &attributes,
                                                                 ScalarArray &dw) {
    if (Base::addLocalNeighbor(w, localQ, attributes, dw)) {
        m_dSumN += attributes.normal() * dw;
        return true;
    }

    return false;
}
