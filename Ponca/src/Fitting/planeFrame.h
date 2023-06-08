#pragma once

#include "./defines.h"

#include PONCA_MULTIARCH_INCLUDE_STD(cmath)
#include PONCA_MULTIARCH_INCLUDE_STD(limits)

namespace Ponca
{
/*!
    \brief Class extending Plane to provide a local frame

    This class can be used to do base change from the global frame to the local frame of the plane.

    \see Plane
    \see worldToLocalFrame
    \see localFrameToWorld
*/
template < class DataPoint, class _WFunctor, typename T >
class PlaneFrame : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        check = Base::PROVIDES_PLANE,  /*!< \brief Requires PrimitiveBase */
        PROVIDES_PLANE_FRAME           /*!< \brief Provides PlaneFrame */
    };

protected:
    // local frame
    VectorType m_u; 
    VectorType m_v;

public:
    PONCA_EXPLICIT_CAST_OPERATORS(PlaneFrame,planeFrame)

    /*! \brief Set the vectors of the local frame to zero 
     *  \param _basisCenter Center of the local frame
     */
    PONCA_MULTIARCH inline void init(const VectorType& _basisCenter = VectorType::Zero())
    {
        Base::init(_basisCenter);
        m_u = VectorType::Zero();
        m_v = VectorType::Zero();
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

}; //class PlaneFrame

#include "planeFrame.hpp"

} // namespace Ponca
