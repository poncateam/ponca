/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#pragma once

/**
  *
  * \defgroup fitting Fitting module
  * \brief This modules includes classes and methods for primitive fitting
  * @{
  * @}
  *
  */
#include "../Common/defines.h"

#define PONCA_EXPLICIT_CAST_OPERATORS(CLASSNAME,CONVERTER)                                                         \
/*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */                   \
PONCA_MULTIARCH inline                                                                                             \
CLASSNAME<DataPoint, _WFunctor, T>& CONVERTER()                                                                    \
{ return * static_cast<CLASSNAME<DataPoint, _WFunctor, T>*>(this); }                                               \
/*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */                   \
PONCA_MULTIARCH inline                                                                                             \
const CLASSNAME<DataPoint, _WFunctor, T>& CONVERTER() const                                                        \
{ return * static_cast<const CLASSNAME<DataPoint, _WFunctor, T>*>(this); }

// CAST OPERATORS

#define PONCA_EXPLICIT_CAST_OPERATORS_DER(CLASSNAME,CONVERTER)                                                     \
/*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */                   \
PONCA_MULTIARCH inline                                                                                             \
CLASSNAME<DataPoint, _WFunctor, DiffType, T>& CONVERTER()                                                          \
{ return * static_cast<CLASSNAME<DataPoint, _WFunctor, DiffType, T>*>(this); }                                     \
/*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */                   \
PONCA_MULTIARCH inline                                                                                             \
const CLASSNAME<DataPoint, _WFunctor, DiffType, T>& CONVERTER() const                                              \
{ return * static_cast<const CLASSNAME<DataPoint, _WFunctor, DiffType, T>*>(this); }


// FIT DEFAULT TYPES

/// Declare the following defaults types: Base, Scalar, VectorType, WFunctor
#define PONCA_FITTING_DECLARE_DEFAULT_TYPES                                                                         \
private:                                                                                                            \
using Base = T;  /*!< \brief Base class of the procedure*/                                                          \
public:                                                                                                             \
using Scalar     = typename DataPoint::Scalar; /*!< \brief Alias to scalar type*/                                   \
using VectorType = typename Base::VectorType;  /*!< \brief Alias to vector type*/                                   \
using WFunctor   = typename Base::WFunctor;    /*!< \brief Alias to weight function*/

/// Declare the following defaults types: Base, Scalar, VectorType, WFunctor
#define PONCA_FITTING_DECLARE_MATRIX_TYPE                                                                           \
public:                                                                                                             \
using MatrixType  = typename DataPoint::MatrixType; /*!< \brief Alias to matrix type*/                              \

/// Declare the following defaults types: Base, Scalar, VectorType, WFunctor
#define PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES                                                                     \
public:                                                                                                             \
using ScalarArray = typename Base::ScalarArray;     /*!< \brief Alias to scalar derivatives array */                \
using VectorArray = typename Base::VectorArray;     /*!< \brief Alias to vector derivatives array */


// FIT API DECLARATION

/// Declare Concept::FittingProcedureConcept::init()
#define PONCA_FITTING_DECLARE_INIT                                                                                 \
/*! \copydoc Concept::FittingProcedureConcept::init() */                                                           \
PONCA_MULTIARCH inline void init (const VectorType& _evalPos);

/// Declare Concept::FittingProcedureConcept::addLocalNeighbor
#define PONCA_FITTING_DECLARE_ADDNEIGHBOR                                                                          \
/*! \copydoc Concept::FittingProcedureConcept::addLocalNeighbor() */                                               \
PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes);

/// Declare Concept::FittingProcedureConcept::addLocalNeighbor
#define PONCA_FITTING_DECLARE_ADDNEIGHBOR_DER                                                                      \
/*! \see Concept::FittingProcedureConcept::addLocalNeighbor() */                                                   \
PONCA_MULTIARCH inline bool                                                                                        \
addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes, ScalarArray &dw);

/// Declare Concept::FittingProcedureConcept::finalize
#define PONCA_FITTING_DECLARE_FINALIZE                                                                             \
/*! \copydoc Concept::FittingProcedureConcept::finalize() */                                                       \
PONCA_MULTIARCH inline FIT_RESULT finalize();

#define PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE                                                                    \
PONCA_FITTING_DECLARE_INIT                                                                                         \
PONCA_FITTING_DECLARE_ADDNEIGHBOR                                                                                  \
PONCA_FITTING_DECLARE_FINALIZE

#define PONCA_FITTING_DECLARE_INIT_ADDDER_FINALIZE                                                                 \
PONCA_FITTING_DECLARE_INIT                                                                                         \
PONCA_FITTING_DECLARE_ADDNEIGHBOR_DER                                                                              \
PONCA_FITTING_DECLARE_FINALIZE
