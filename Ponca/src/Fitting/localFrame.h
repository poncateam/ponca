#pragma once

#include "./defines.h"

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include PONCA_MULTIARCH_INCLUDE_STD(limits)

namespace Ponca
{
/*!
    \brief Class provides a local Frame

    This class can be used to do base change from the global frame to a local frame described by its axis m_u and m_v.
    The local frame is defined by the three vectors m_u, m_v and Base::primitiveGradient().

    It could be used to do the translation/rotation of a point.

    \see worldToLocalFrame
    \see localFrameToWorld
*/
template < class DataPoint, class _WFunctor, typename T >
class LocalFrame : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

private:

    // local frame
    VectorType m_u; 
    VectorType m_v;

protected:
    enum
    {
        PROVIDES_LOCAL_FRAME           /*!< \brief Provides LocalFrame */
    };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(LocalFrame,localFrame)

    /*! \brief Set the vectors of the local frame to zero
     */
    PONCA_MULTIARCH inline void init()
    {
        Base::init();
        m_u = VectorType::Zero();
        m_v = VectorType::Zero();
    }

    /*!
     * \brief Return the first axis of the local frame
     */
    PONCA_MULTIARCH inline const VectorType& getFrameU() const { return m_u; }

    /*!
     * \brief Return the second axis of the local frame
     */
    PONCA_MULTIARCH inline const VectorType& getFrameV() const { return m_v; }

    /*!
     * \brief Set the axis of the local frame
     *
     * The frame will be defined by the two vectors m_u and m_v and the normal given by Base::primitiveGradient().
     * At the end, we have the basis B = [Base::primitiveGradient(), _u, _v].
     * \param _u First axis of the local frame
     * \param _v Second axis of the local frame
     */
    PONCA_MULTIARCH inline void setFrameUV(const VectorType& _u, const VectorType& _v) { 
        m_u = _u;
        m_v = _v;
    }

    /*!
     * \brief Express a point in ambient space relatively to the local frame.
     *
     * Output vector is: [h, u, v]^T, where u, v are 2d coordinates on the plane,
     * and h the height of the sample.
     * \tparam ignoreTranslation must be set to true when passing vectors instead of points
     * \param _q Point coordinates expressed in ambient space
     * \return Point coordinates expressed in local frame
     */
    template <bool ignoreTranslation = false>
    PONCA_MULTIARCH inline VectorType worldToLocalFrame(const VectorType &_q) const;

    /*!
     * \brief Transform a point from the local frame [h, u, v]^T to ambient space
     *
     * \tparam ignoreTranslation must be set to true when passing vectors instead of points
     * \param _q Point coordinates expressed in local frame
     * \return Point coordinates expressed in ambient space
     */
    template <bool ignoreTranslation = false>
    PONCA_MULTIARCH inline VectorType localFrameToWorld(const VectorType &_q) const;

}; //class LocalFrame

#include "localFrame.hpp"

} // namespace Ponca
