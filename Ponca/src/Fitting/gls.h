/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"
#include "./curvature.h"

namespace Ponca
{

/*!
    \brief Growing Least Squares reparemetrization of the OrientedSphereFit
    \inherit Concept::FittingExtensionConcept

    Method published in \cite Mellado:2012:GLS

    This class assumes that the WeightFunc defines the accessor
    \code
        w.evalScale();
    \endcode
    in order to access to the evaluation scale, needed to compute the
    scale invariant GLS reparametrization (all *_normalized methods).

    Computed values:
    - tau(), eta() and kappa(): the GLS descriptor
    \f$ \left[ \tau \; \eta \; \kappa \right]\f$
    - tau_normalized(), eta_normalized() and kappa_normalized():
    the scale invariant GLS descriptor
    \f$ \left[ \frac{\tau}{t} \; \eta \; t\kappa \right]\f$

    Requirements:
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim
    Provides:
    \verbatim PROVIDES_GLS_PARAMETRIZATION \endverbatim
*/
template < class DataPoint, class _WFunctor, typename T>
class GLSParam : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES

//! [Requirements]
protected:
    enum
    {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE,
        PROVIDES_GLS_PARAMETRIZATION
    };
//! [Requirements]

protected:
    Scalar m_fitness {0};   /*!< \brief Save the fitness value to avoid side effect with Pratt normalization*/

public:
    PONCA_EXPLICIT_CAST_OPERATORS(GLSParam,glsParam)
    PONCA_FITTING_DECLARE_FINALIZE

    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
    /*! \brief Compute and return \f$ \tau \f$ */
    PONCA_MULTIARCH inline Scalar tau() const
    {
        return Base::isNormalized() ? Base::m_uc : Base::m_uc / Base::prattNorm();
    }

    /*! \brief Compute and return \f$ \eta \f$ */
    PONCA_MULTIARCH inline VectorType eta() const { return Base::primitiveGradient(); }

    /*! \brief Compute and return \f$ \kappa \f$ */
    PONCA_MULTIARCH inline Scalar kappa() const
    {
        return Scalar(2.) * (Base::isNormalized() ? Base::m_uq : Base::m_uq / Base::prattNorm());
    }

    /*! \brief Compute and return \f$ \frac{\tau}{t} \f$ */
    PONCA_MULTIARCH inline Scalar tau_normalized() const { return tau() / Base::m_w.evalScale(); }

    /*! \brief Compute and return \f$ \eta \f$ */
    PONCA_MULTIARCH inline VectorType eta_normalized() const { return eta(); }

    /*! \brief Compute and return \f$ t \kappa \f$ */
    PONCA_MULTIARCH inline Scalar kappa_normalized() const { return kappa() * Base::m_w.evalScale(); }

    /*! \brief Return the fitness, e.g. the pratt norm of the initial scalar field */
    PONCA_MULTIARCH inline Scalar fitness() const { return m_fitness; }

    /*!
    \brief Compare current instance with other.
    \return a distance between two fits (0 correspond to two similar fits)
    \warning Use the same scale to have a useful comparison (normalized value are used)
    */
    PONCA_MULTIARCH inline Scalar compareTo (const GLSParam<DataPoint, _WFunctor, T>& _other,
                                        bool _useFitness = true) const
    {
        Scalar nTau     = this->tau_normalized()   - _other.tau_normalized();
        Scalar nKappa   = this->kappa_normalized() - _other.kappa_normalized();
        Scalar nFitness = _useFitness ? this->fitness() - _other.fitness() : Scalar(0.);

        return nTau * nTau + nKappa * nKappa + nFitness * nFitness;
    }

}; //class GLSParam


/*!
    \brief Differentiation of GLSParam
    \inherit Concept::FittingExtensionConcept

    Method published in \cite Mellado:2012:GLS
*/
template < class DataPoint, class _WFunctor, int DiffType, typename T>
class GLSDer : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES

protected:
    enum
    {
        Check = Base::PROVIDES_GLS_PARAMETRIZATION &
                Base::PROVIDES_PRIMITIVE_DERIVATIVE &
                Base::PROVIDES_ALGEBRAIC_SPHERE_DERIVATIVE,
        PROVIDES_GLS_DERIVATIVE,
        PROVIDES_GLS_GEOM_VAR
    };

public:
    PONCA_EXPLICIT_CAST_OPERATORS_DER(GLSDer,glsDer)

    PONCA_MULTIARCH inline ScalarArray dtau()   const; /*!< \brief Compute and return \f$ \tau \f$ derivatives */
    PONCA_MULTIARCH inline VectorArray deta()   const; /*!< \brief Compute and return \f$ \eta \f$ derivatives */
    PONCA_MULTIARCH inline ScalarArray dkappa() const; /*!< \brief Compute and return \f$ \kappa \f$ derivatives */

    PONCA_MULTIARCH inline ScalarArray dtau_normalized()   const; /*!< \brief Compute and return \f$ \tau \f$ derivatives */
    PONCA_MULTIARCH inline VectorArray deta_normalized()   const; /*!< \brief Compute and return \f$ t * d\eta \f$ */
    PONCA_MULTIARCH inline ScalarArray dkappa_normalized() const; /*!< \brief Compute and return \f$ d\kappa * t^{2} \f$ */

    /*!
       \brief The Geometric Variation is computed as the weighted sum of the GLS scale-invariant partial derivatives
       \f[
        \nu(\mathbf{p},t) =
        w_\tau   \left(\frac{\delta\tau}{\delta t}\right)^2 +
        w_\eta   \left( t   \frac{\delta\eta}{\delta t}\right)^2 +
        w_\kappa \left( t^2 \frac{\delta\kappa}{\delta t}\right)^2
       \f]

       Method published in \cite Mellado:2012:GLS
       @param wtau Weight applied to \f$ \tau \f$
       @param weta Weight applied to \f$ \eta \f$
       @param wkappa Weight applied to \f$ \kappa \f$
       @return
     */
    PONCA_MULTIARCH inline Scalar geomVar(Scalar wtau   = Scalar(1),
                                          Scalar weta   = Scalar(1),
                                          Scalar wkappa = Scalar(1)) const;
}; //class GLSScaleDer


#include "gls.hpp"

} //namespace Ponca
