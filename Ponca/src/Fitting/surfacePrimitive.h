/*
 Copyright (C) 2014 Nicolas Mellado <nmellado0@gmail.com>
Copyright (C) 2021 aniket agarwalla <aniketagarwalla37@gmail.com>


 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./primitive.h"
#include <Eigen/Geometry>

namespace Ponca
{

    /*!
        \brief Implicit surface defined by an homogeneous vector \f$\mathbf{p}\f$.

        In n-dimensionnal space, the surface is defined as
        the \f$0\f$-isosurface of the scalar field

        \f$ s_\mathbf{u}(\mathbf{x}) =
        \left[ \mathbf{x}^T \; 1 \;\right]^T \cdot \mathbf{p} \f$.


        This primitive requires the definition of n-dimensionnal vectors
        (VectorType) in Concept::PointConcept.

        \ingroup fitting

    */
    template < class DataPoint, class _WFunctor, typename T = void  >
    class Surface : public PrimitiveBase<DataPoint, _WFunctor>
    {
    private:

        using Base      = PrimitiveBase<DataPoint, _WFunctor>;

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

        //! \brief Evaluation position (needed for centered basis)
        VectorType m_p;
    // results
    private:

        Eigen::Matrix<Scalar, 9, 1>  q_para;
        VectorType q_point;   


    public:

        /*! \brief Default constructor */
        PONCA_MULTIARCH inline Surface()
            : Base()
        {
            m_p = VectorType::Zero();
            resetPrimitive();
        }

        /*! \brief Explicit conversion to Surface, to access methods potentially hidden by inheritage */
        PONCA_MULTIARCH inline
        Surface<DataPoint, WFunctor, T>& surface(){ }

        /*! \brief Set the scalar field values to 0
            status */
        PONCA_MULTIARCH inline void resetPrimitive()
        {
            Base::resetPrimitive();
            q_para = Eigen::Matrix<Scalar, 9, 1> ::Zero();
            q_point= VectorType::Zero();
        }

        PONCA_MULTIARCH inline bool operator==(const Surface<DataPoint, WFunctor, T>& other) const{
            return q_para == other.getParamters();
        }

        /*! \brief Comparison operator, convenience function */
        PONCA_MULTIARCH inline bool operator!=(const Surface<DataPoint, WFunctor, T>& other) const{
            return ! ((*this) == other);
        }

        /*! \brief Reading access to the basis center (evaluation position) */
        PONCA_MULTIARCH inline const VectorType& basisCenter () const { return m_p; }
        /*! \brief Writing access to the (evaluation position) */
        PONCA_MULTIARCH inline       VectorType& basisCenter ()       { return m_p; }

        PONCA_MULTIARCH inline       Eigen::Matrix<Scalar, 9, 1>& getParameters ()  { return q_para; }
        PONCA_MULTIARCH inline       VectorType& getPoint ()       { return q_point; }

        /* \brief Init the surface from a direction and a position
        \param _dir Orientation of the surface, does not need to be normalized
        \param _pos Position of the surface
        */
        PONCA_MULTIARCH inline void setSurface (const Eigen::Matrix<Scalar, 9, 1> & _dir,
                                        const VectorType& _pos)
        {
        q_para = _dir;
        q_point = _pos;
        }
        
    }; //class Surface

}