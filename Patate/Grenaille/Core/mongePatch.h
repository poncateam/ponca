/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef _GRENAILLE_MONGE_PATCH_
#define _GRENAILLE_MONGE_PATCH_

#include <iostream>


namespace Grenaille
{



/*!
 * \brief Extension to compute the best fit quadric on 3d points expressed as f(u,v)=h
 *
 * \note This procedure requires at least two passes, the first one for plane fitting,
 * the second one for quadric fitting.
 * \warning This class is valid only in 3D.
 */
template < class DataPoint, class _WFunctor, typename T>
class MongePatch : public T
{
private:
    using Base = T;

protected:
    enum
    {
        Check = Base::PROVIDES_PLANE && Base::PROVIDES_TANGENT_PLANE_BASIS
    };

public:
    typedef typename Base::Scalar          Scalar;      /*!< \brief Inherited scalar type*/
    typedef typename Base::VectorType      VectorType;  /*!< \brief Inherited vector type*/
    typedef typename DataPoint::MatrixType MatrixType;  /*!< \brief Matrix type inherited from DataPoint*/
    typedef _WFunctor                       WFunctor;   /*!< \brief Weight Function */

    typedef Eigen::Matrix<Scalar,2,1> Vector2;
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> VectorX;
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> MatrixX;

protected:
    MatrixX m_A; /*!< \brief Quadric input samples */
    MatrixX m_x; /*!< \brief Quadric parameters */
    VectorX m_b;         /*!< \brief Obervations */

    int  m_neiIdx;       /*!< \brief Counter of observations, used in addNeighhor() */
    bool m_planeIsReady;
public:

    /*! \brief Explicit conversion to MongePatch, to access methods potentially hidden by inheritage */
    MULTIARCH inline
    MongePatch<DataPoint, WFunctor, T>& mongePatch()
    { return * static_cast<MongePatch<DataPoint, WFunctor, T>*>(this); }
    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::init() */
    MULTIARCH inline void init (const VectorType& _evalPos);

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::addNeighbor() */
    MULTIARCH inline bool addNeighbor(const DataPoint &_nei);

    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    MULTIARCH inline FIT_RESULT finalize();


    //! \brief Returns an estimate of the mean curvature
    MULTIARCH inline Scalar kMean() const;

    //! \brief Returns an estimate of the Gaussian curvature
    MULTIARCH inline Scalar GaussianCurvature() const;

    MULTIARCH inline Scalar evalUV(Scalar u, Scalar v) const {
      return h_uu()*u*u + h_vv()*v*v + h_uv()*u*v + h_u()*u + h_v()*v + h_c();
    }

    /*! \brief Value of the scalar field at the evaluation point */
    MULTIARCH inline Scalar potential(const VectorType& _q) const {
      VectorType x = Base::worldToTangentPlane(_q);
      return evalUV(*(x.data()+1),*(x.data()+2)) - *(x.data());
    }

    //! \brief Orthogonal projecting on the patch, such that h = f(u,v)
    MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        VectorType x = Base::worldToTangentPlane(_q);
        *(x.data()) = evalUV(*(x.data()+1),*(x.data()+2));
        return Base::tangentPlaneToWorld(x);
    }

    MULTIARCH inline const Scalar & h_uu () const { return *(m_x.data()); }
    MULTIARCH inline const Scalar & h_vv () const { return *(m_x.data()+1); }
    MULTIARCH inline const Scalar & h_uv () const { return *(m_x.data()+2); }
    MULTIARCH inline const Scalar & h_u  () const { return *(m_x.data()+3); }
    MULTIARCH inline const Scalar & h_v  () const { return *(m_x.data()+4); }
    MULTIARCH inline const Scalar & h_c  () const { return *(m_x.data()+5); }

};

#include "mongePatch.hpp"

} //namespace Grenaille

#endif
