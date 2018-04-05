



#ifndef _GRENAILLE_MLS_SPHERE_FIT_DER_
#define _GRENAILLE_MLS_SPHERE_FIT_DER_



namespace Grenaille
{

//TODO(thib) get Type from *SphereDer class
//TODO(thib) or should it be put it in internal:: namespace ?
//           and create other public class as OrientedSphereDer does

template < class DataPoint, class _WFunctor, typename T>
class MlsSphereFitDer : public T
{
private:
    typedef T Base;

protected:
    enum
    {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE_DERIVATIVE,
        PROVIDES_NORMAL_SPACE_DERIVATIVE
    };

public:
    typedef typename Base::Scalar     Scalar;     /*!< \brief Inherited scalar type*/
    typedef typename Base::VectorType VectorType; /*!< \brief Inherited vector type*/

    typedef typename Base::VectorArray VectorArray; /*!< \brief Inherited vector array type */
    typedef typename Base::ScalarArray ScalarArray; /*!< \brief Inherited scalar array type */

    //TODO(thib) redundant with OrientedSphereDer macro, which cannot be used here due to internal namespace...
#define DER_NB_DERIVATIVES(TYPE,DIM) ((TYPE & internal::FitScaleDer) ? 1 : 0 ) + ((TYPE & internal::FitSpaceDer) ? DIM : 0)

    //TODO(thib) use appropriate storage order
    //#define GLS_DER_STORAGE_ORDER(TYPE)      ((TYPE & FitSpaceDer) ? Eigen::RowMajor : Eigen::ColMajor )

    typedef Eigen::Matrix< Scalar,
        DER_NB_DERIVATIVES(Base::Type,DataPoint::Dim),
        DER_NB_DERIVATIVES(Base::Type,DataPoint::Dim)> Matrix;

    typedef Eigen::Matrix< Scalar,
        DER_NB_DERIVATIVES(Base::Type,DataPoint::Dim),
        DataPoint::Dim*DER_NB_DERIVATIVES(Base::Type,DataPoint::Dim)> MatrixArray;

protected:
    // computation data
    Matrix m_d2SumDotPN, /*!< \brief Sum of the dot product betwen relative positions and normals with twice differenciated weights */
           m_d2SumDotPP, /*!< \brief Sum of the squared relative positions with twice differenciated weights */
           m_d2SumW;     /*!< \brief Sum of queries weight with twice differenciated weights */

public:
    // results
    Matrix m_d2Uc,      /*!< \brief Second derivative of the hyper-sphere constant term  */
           m_d2Uq;      /*!< \brief Second derivative of the hyper-sphere quadratic term */
    MatrixArray m_d2Ul; /*!< \brief Second derivative of the hyper-sphere linear term    */

public:
    /************************************************************************/
    /* Initialization                                                       */
    /************************************************************************/
    /*! \see Concept::FittingProcedureConcept::init() */
    MULTIARCH void init(const VectorType &evalPos);

    /************************************************************************/
    /* Processing                                                           */
    /************************************************************************/
    /*! \see Concept::FittingProcedureConcept::addNeighbor() */
    MULTIARCH bool addNeighbor(const DataPoint &nei);

    /*! \see Concept::FittingProcedureConcept::finalize() */
    MULTIARCH FIT_RESULT finalize();

    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/

    /*! \brief Returns the derivatives of the scalar field at the evaluation point */
    MULTIARCH inline ScalarArray dPotential() const;

    /*! \brief Returns the second derivatives of the scalar field at the evaluation point */
    MULTIARCH inline VectorArray dNormal() const;



}; //class MlsSphereFitDer

#include "mlsSphereFitDer.hpp"

} //namespace Grenaille

#endif
