/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


template<class DataPoint, class _WFunctor, typename T>
void RimlsPlaneFitImpl<DataPoint, _WFunctor, T>::init(const VectorType& _evalPos)
{
    Base::init( _evalPos );

    // Setup internal values
    m_minConvergence = Scalar(1e-3);
    m_maxIteration = 10;
    m_sigmaN = Scalar( 0.5 );
}

template<class DataPoint, class _WFunctor, typename T>
bool RimlsPlaneFitImpl<DataPoint, _WFunctor, T>::addLocalNeighbor(Scalar w,
                                                                    const VectorType &localQ,
                                                                    const DataPoint &attributes)
{

    //if valid, the neighbor is added to the neighbors vector
    if( Base::addLocalNeighbor( w, localQ, attributes ) ){
        m_neighbors.push_back( attributes ); 
        return true;
    }
    return false;
}

template<class DataPoint, class _WFunctor, typename T>
FIT_RESULT RimlsPlaneFitImpl<DataPoint, _WFunctor, T>::finalize()
{

    // handle UNDEFINED cases
    if( Base::finalize() != STABLE ) {
        return Base::m_eCurrentState = UNDEFINED;
    }

    // handle conflict error
    if( Base::plane().isValid() ){
        return Base::m_eCurrentState = CONFLICT_ERROR_FOUND;
    }

    unsigned int iteration{ 0 };
    Scalar convergence = Scalar(1);

    Scalar f = Scalar(0.);
    VectorType gradF = VectorType::Zero();
    VectorType prevGrad = VectorType::Zero();
    Scalar alpha = Scalar(1.);

    VectorType px;
    Scalar fx;
    Scalar w;
    VectorType gradW;

    while( iteration < m_maxIteration && convergence > m_minConvergence ) {

        // initialise values
        Scalar sumF = Scalar(0);
        Scalar sumW = Scalar(0);
        VectorType sumGF = VectorType::Zero();
        VectorType sumGW = VectorType::Zero();
        VectorType sumN = VectorType::Zero();

        for(DataPoint& nei : m_neighbors) {

            px = Base::m_w.evalPos() - nei.pos();
            fx = px.dot( nei.normal() );

            // use gaussian weights
            if( iteration > 0 ) {
                alpha = gaussian( ( nei.normal() - gradF ).squaredNorm(), m_sigmaN );
                alpha *= gaussian( (fx - f) * (fx - f), 0.5 * Base::m_w.evalScale() );
            }

            w = alpha * Base::m_w.w( nei.pos(), nei ).first;
            gradW = alpha * 2 * px * Base::m_w.scaledw(px, nei);

            sumW += w;
            sumGW += gradW;
            sumF += w * fx;
            sumGF += gradW * fx;
            sumN += w * nei.normal();
        }

        if( sumW == 0 ){
            return Base::m_eCurrentState = UNDEFINED;
        }

        prevGrad = gradF;

        // update the actual potential and gradient
        f = sumF / sumW;
        gradF = ( sumGF - f * sumGW + sumN ) / sumW;
        convergence = ( prevGrad - gradF ).squaredNorm();
        ++iteration;
		
    }

    Base::setPlane( -gradF, -f * gradF );
    return Base::m_eCurrentState = STABLE;
}
