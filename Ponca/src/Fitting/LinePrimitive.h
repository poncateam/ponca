/*
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



template < class DataPoint, class _WFunctor, typename T = void  >
class Line : public PrimitiveBase<DataPoint, _WFunctor>
{
private:

   typedef PrimitiveBase<DataPoint, _WFunctor> Base;

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
    /*! \brief  The equation of a line with a point vectorType point 
    and direction vectorType direction */

    VectorType _direction;
    VectorType _point; 

private:

    //! \brief Evaluation position (needed for centered basis)
    VectorType m_p;

public:

    /*! \brief Default constructor */
    PONCA_MULTIARCH inline Line()
        : Base()
    {
        resetPrimitive();
    }


    /*! \brief Explicit conversion to Line, to access methods potentially hidden by inheritage */
    PONCA_MULTIARCH inline
    Line<DataPoint, WFunctor, T>& line()
    { return * static_cast<Line<DataPoint, WFunctor, T>*>(this); }

    /*! \brief Set the scalar field values to 0 and reset the isNormalized() status

        \warning Set a_y to Zero(), which leads to nans in OrientedLine::normal()
        \FIXME Set and use Base::m_state to handle invalid configuration
    */
    PONCA_MULTIARCH inline void resetPrimitive()
    {
        Base::resetPrimitive();
        _direction = VectorType::Zero();
        _point = VectorType::Zero();

    }

    PONCA_MULTIARCH inline bool operator==(const Line<DataPoint, WFunctor, T>& other) const{
        return direction == other.direction && point == other.point;
              
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
                                    const VectorType& _axis)
    {
        _point = _origin;
        _direction = _axis;
    }

     /*! \brief direction of the  lines */
    PONCA_MULTIARCH inline VectorType direction () const
    {
        /* Uniform direction defined only by the orientation of the line*/
        return _direction;
    }
     /*! \brief point on the fitting line */
    PONCA_MULTIARCH inline VectorType point () const
    {
        /* Uniform point defined only by the orientation of the line */
        return _point;
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
