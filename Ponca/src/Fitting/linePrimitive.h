
/*
 Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>
 
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma onces

#include "./defines.h"
#include "./primitive.h" // PrimitiveBase
#include <Eigen/Geometry>
#include <Eigen/Core>

namespace Ponca
{

/*!
    \brief A parametrized line is defined by an origin point o and a unit direction vector d such that the line corresponds to the set l(t)=o+td, t∈R.

    In 3-dimensionnal space, the line has an orgin and a direction vector.

    This class inherits Eigen::ParametrizedLine.

    This primitive requires the definition of n-dimensionnal vectors
    (VectorType) in Concept::PointConcept.

    This primitive provides:
    \verbatim PROVIDES_LINE \endverbatim

    \ingroup fitting
*/

template < class DataPoint, class _WFunctor, typename T = void  >
class Line : public PrimitiveBase<DataPoint, _WFunctor>,
             public Eigen::ParametrizedLine<typename DataPoint::Scalar, DataPoint::Dim >
{
private:

   typedef PrimitiveBase<DataPoint, _WFunctor> Base;

public:
    /// \brief Specialization of Eigen::ParametrizedLine inherited by Ponca::Line
    using EigenBase = Eigen::ParametrizedLine<typename DataPoint::Scalar, DataPoint::Dim >;

protected:

    enum
    {
        PROVIDES_LINE /*!< \brief Provides  Line */
    };

public:

    /*! \brief Scalar type inherited from DataPoint */
    typedef typename DataPoint::Scalar      Scalar;
    /*! \brief Vector type inherited from DataPoint */
    typedef typename DataPoint::VectorType  VectorType;
    /*! \brief Matrix type inherited from DataPoint */
    typedef typename DataPoint::MatrixType  MatrixType;
    /*! \brief Weight Function */
    typedef _WFunctor                       WFunctor;

private:

    /*! \brief Evaluation position (needed for centered basis) */
    VectorType m_p;

public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline Line()
        : Base()
    {
        m_p = VectorType::Zero();
        resetPrimitive();
    }


    /*! \brief Explicit conversion to Line, to access methods potentially hidden by inheritage */
    PONCA_MULTIARCH inline
    Line<DataPoint, WFunctor, T>& line()
    { return * static_cast<Line<DataPoint, WFunctor, T>*>(this); }

    /*! \brief Set the scalar field values to 0 and reset the distance() and origin() status

        \FIXME Set and use Base::m_state to handle invalid configuration
    */
    PONCA_MULTIARCH inline void resetPrimitive()
    {
        Base::resetPrimitive();
        EigenBase* cc = static_cast<EigenBase*>(this);
        *cc = EigenBase();

    }
    /*! \brief Comparison operator */
    PONCA_MULTIARCH inline bool operator==(const Line<DataPoint, WFunctor, T>& other) const{
        return EigenBase::isApprox(other);
    }

    /*! \brief Comparison operator, convenience function */
    PONCA_MULTIARCH inline bool operator!=(const Line<DataPoint, WFunctor, T>& other) const{
        return ! ((*this) == other);
    }
      /*! \brief Reading access to the basis center (evaluation position) */
    PONCA_MULTIARCH inline const VectorType& basisCenter () const { return m_p; }
    /*! \brief Writing access to the (evaluation position) */
    PONCA_MULTIARCH inline       VectorType& basisCenter ()       { return m_p; }


    /* \brief Init the line from a direction and a position
       \param _dir Orientation of the line, does not need to be normalized
       \param _pos Position of the line
    */
    PONCA_MULTIARCH inline void setLine (const VectorType& _origin,
                                    const VectorType& _direction)
    {
        EigenBase* cc = static_cast<EigenBase*>(this);
        *cc = EigenBase(_origin, _direction);
    }

     //! \brief direction of the  lines
    PONCA_MULTIARCH inline VectorType direction () const
    {
        // Uniform direction defined only by the orientation of the line
        return EigenBase::direction();
    }
     //! \brief point on the fitting line
    PONCA_MULTIARCH inline VectorType point () const
    {
        // Uniform point defined only by the orientation of the line
        return EigenBase::origin();
    }

    //! \brief Value of the scalar field at the location \f$ \mathbf{q} \f$
    PONCA_MULTIARCH inline Scalar distance (const VectorType& _q) const
    {
        // The potential is the distance from a point to the line
        return EigenBase::distance(_q - m_p);
    }

    //! \brief Project a point on the line
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        // Project on the normal vector and add the offset value
        return EigenBase::projection(_q - m_p) + m_p;
    }


    /*!
        \brief Used to know if the fitting result to a line
        \return true if finalize() have been called and the fitting result to a line
    */
    PONCA_MULTIARCH inline bool isLine() const
    {
        
        bool bReady    = Base::isReady();
        return bReady;
    }

}; //class Line


}
