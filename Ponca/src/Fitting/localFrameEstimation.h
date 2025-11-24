#pragma once

#include "./defines.h"

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include PONCA_MULTIARCH_INCLUDE_STD(limits)

namespace Ponca
{
/*!
    \brief Class provides a local Frame Estimation

    This class can be used to estimate a pair of orthogonal vectors defining a local frame, given a normal vector.

    \see LocalFrame
*/
template < class DataPoint, class _NFilter, typename T >
class LocalFrameEstimator : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        Check= Base::PROVIDES_LOCAL_FRAME /*!< \brief Check LocalFrame */
    };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(LocalFrameEstimator, localFrameEstimator)

    /*!
     * \brief Given a normal vector, this function computes a local coordinate frame (orthogonal basis).
     *
     * It starts by generating a non-collinear vector to the normal (to ensure that the their cross product is non-zero). 
     * The vector, obtained via the cross product, becomes the first axis (m_u) of the frame.
     * The second axis (m_v) is computed by taking the cross product between the normal and the first axis.
     * At the end, we obtain the basis B = [Base::primitiveGradient(), m_u, m_v].
     * \param _norm Normal vector for which the local frame is computed.
     */
    PONCA_MULTIARCH inline void computeFrameFromNormalVector(const VectorType& _norm){
        // Creation of the vector 'a' non-collinear to the normal vector
        VectorType a;
        if (std::abs(_norm.x()) > std::abs(_norm.z())) {
            a = VectorType(-_norm.y(), _norm.x(), 0);
        } else {
            a = VectorType(0, -_norm.z(), _norm.y());
        }
        a.normalize();
        // Creation of the two vectors of the local frame (m_u and m_v) thanks to the cross product
        VectorType m_u = _norm.cross(a);
        VectorType m_v = _norm.cross(m_u);
        m_u.normalize();
        m_v.normalize();
        Base::setFrameUV (m_u, m_v);
    }

}; //class LocalFrameEstimator

} // namespace Ponca
