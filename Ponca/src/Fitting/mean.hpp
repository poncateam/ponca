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
    Base::init(_evalPos);
    m_sumP = VectorType::Zero();
}

template<class DataPoint, class _WFunctor, typename T>
bool
MeanPosition<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
                                                        const VectorType &localQ,
                                                        const DataPoint &attributes) {
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        m_sumP += w * localQ;
        return true;
    }
    return false;
}

template < class DataPoint, class _WFunctor, typename T>
void
MeanNormal<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init(_evalPos);
    m_sumN = VectorType::Zero();
}


template<class DataPoint, class _WFunctor, typename T>
bool
MeanNormal<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
                                                      const VectorType &localQ,
                                                      const DataPoint &attributes) {
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        m_sumN += w * attributes.normal();
        return true;
    }
    return false;
}


namespace internal {

    template<class DataPoint, class _WFunctor, typename T, int Type>
    void
    MeanPositionDer<DataPoint, _WFunctor, T, Type>::init(const VectorType &_evalPos) {
        Base::init(_evalPos);
        m_dSumP.setZero();
    }


    template<class DataPoint, class _WFunctor, typename T, int Type>
    bool
    MeanPositionDer<DataPoint, _WFunctor, T, Type>::addLocalNeighbor(Scalar w,
                                                                     const VectorType &localQ,
                                                                     const DataPoint &attributes,
                                                                     ScalarArray &dw) {
        if (Base::addLocalNeighbor(w, localQ, attributes, dw)) {
            int spaceId = (Type & FitScaleDer) ? 1 : 0;

            m_dSumP += localQ * dw;
            return true;
        }

        return false;
    }
}
