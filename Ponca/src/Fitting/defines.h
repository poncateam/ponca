/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <type_traits> // std::true_type

/**
 *
 * \defgroup fitting Fitting module
 * \brief This modules includes classes and methods for primitive fitting
 * \{
 * \}
 *
 */
#include "../Common/defines.h"

#define PONCA_EXPLICIT_CAST_OPERATOR(CONVERTER, ...)                                           \
    /*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */ \
    PONCA_MULTIARCH inline __VA_ARGS__& CONVERTER()                                                    \
    {                                                                                                \
        return *static_cast<__VA_ARGS__*>(this);                                                       \
    }                                                                                                \
    /*! \brief Explicit conversion to CLASSNAME, to access methods potentially hidden by heritage */ \
    PONCA_MULTIARCH inline const __VA_ARGS__& CONVERTER() const                                        \
    {                                                                                                \
        return *static_cast<const __VA_ARGS__*>(this);                                                 \
    }


#define PONCA_EXPLICIT_BASKET_CAST_OPERATORS(CLASSNAME, CONVERTER) \
    PONCA_EXPLICIT_CAST_OPERATOR(CONVERTER, CLASSNAME<DataPoint, _NFilter, T>)

// CAST OPERATORS

#define PONCA_EXPLICIT_BASKET_CAST_OPERATORS_DER(CLASSNAME, CONVERTER) \
    PONCA_EXPLICIT_CAST_OPERATOR(CONVERTER, CLASSNAME<DataPoint, _NFilter, DiffType, T>)

// FIT DEFAULT TYPES

/// Declare the following defaults types: Base, Scalar, VectorType, NeighborFilter
#define PONCA_FITTING_DECLARE_DEFAULT_TYPES                                                  \
protected:                                                                                   \
    using Base = T; /*!< \brief Base class of the procedure*/                                \
public:                                                                                      \
    using Scalar         = typename DataPoint::Scalar;    /*!< \brief Alias to scalar type*/ \
    using VectorType     = typename Base::VectorType;     /*!< \brief Alias to vector type*/ \
    using NeighborFilter = typename Base::NeighborFilter; /*!< \brief Alias to the filter applied on the neighbors */

/// Declare the following defaults types: Base, Status, Scalar, VectorType, NeighborFilter
#define PONCA_FITTING_DECLARE_MULTIPASS_TYPES                                                         \
protected:                                                                                            \
    using Base   = MultipassStatus<DataPoint, _NFilter, T>; /*!< \brief Base class of the procedure*/ \
    using Status = Base;                                    /*!< \brief Status class for multipass */ \
public:                                                                                               \
    using Scalar         = typename DataPoint::Scalar;    /*!< \brief Alias to scalar type*/          \
    using VectorType     = typename Base::VectorType;     /*!< \brief Alias to vector type*/          \
    using NeighborFilter = typename Base::NeighborFilter; /*!< \brief Alias to the filter applied on the neighbors */

/// Declare the following defaults types: Base, Scalar, VectorType, NeighborFilter
#define PONCA_FITTING_DECLARE_MATRIX_TYPE \
public:                                   \
    using MatrixType = typename DataPoint::MatrixType; /*!< \brief Alias to matrix type*/

/// Declare the following defaults types: Base, Scalar, VectorType, NeighborFilter
#define PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES                                                      \
public:                                                                                              \
    using ScalarArray = typename Base::ScalarArray; /*!< \brief Alias to scalar derivatives array */ \
    using VectorArray = typename Base::VectorArray; /*!< \brief Alias to vector derivatives array */

// FIT API DOCUMENTATION
#define PONCA_FITTING_APIDOC_SETWFUNC
/*! Init the WeightFunc, without changing the other internal states. Calls #startNewPass internally. \warning Must be \
 * called be for any computation (and before #init). \see getWeightFunc */
#define PONCA_FITTING_APIDOC_INIT
/*! Set the evaluation position and reset the internal states. \warning Must be called be for any computation (but \
 * after #setNeighborFilter) */
#define PONCA_FITTING_APIDOC_ADDNEIGHBOR     /*! Add a neighbor to perform the fit */
#define PONCA_FITTING_APIDOC_ADDNEIGHBOR_DER /*! Add a neighbor to perform the fit */
#define PONCA_FITTING_APIDOC_FINALIZE
/*! Finalize the procedure \return Fitting Status \warning Must be called be for any use of the fitting output */

// FIT API DECLARATION

/// Declare Primitive::isiSigned()
#define PONCA_FITTING_IS_SIGNED(IS_SIGNED)                                                             \
    /*! \brief Is scalar field signed. If not, the method the sign of `potential()` must be ignored */ \
    PONCA_MULTIARCH inline constexpr bool isSigned() { return IS_SIGNED; }

/// Declare Concept::ComputationalObjectConcept::init()
#define PONCA_FITTING_DECLARE_INIT \
    PONCA_FITTING_APIDOC_INIT      \
    PONCA_MULTIARCH inline void init();

/// Declare Concept::ComputationalObjectConcept::addLocalNeighbor
#define PONCA_FITTING_DECLARE_ADDNEIGHBOR \
    PONCA_FITTING_APIDOC_ADDNEIGHBOR      \
    PONCA_MULTIARCH inline void addLocalNeighbor(Scalar w, const VectorType& localQ, const DataPoint& attributes);

/// Declare Concept::ComputationalDerivativesConcept::addLocalNeighbor
#define PONCA_FITTING_DECLARE_ADDNEIGHBOR_DER                                                                     \
    PONCA_FITTING_APIDOC_ADDNEIGHBOR_DER                                                                          \
    PONCA_MULTIARCH inline void addLocalNeighbor(Scalar w, const VectorType& localQ, const DataPoint& attributes, \
                                                 ScalarArray& dw);

/// Declare Concept::ComputationalObjectConcept::finalize
#define PONCA_FITTING_DECLARE_FINALIZE \
    PONCA_FITTING_APIDOC_FINALIZE      \
    PONCA_MULTIARCH inline FIT_RESULT finalize();

#define PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE \
    PONCA_FITTING_DECLARE_INIT                  \
    PONCA_FITTING_DECLARE_ADDNEIGHBOR           \
    PONCA_FITTING_DECLARE_FINALIZE

#define PONCA_FITTING_DECLARE_INIT_ADDDER_FINALIZE \
    PONCA_FITTING_DECLARE_INIT                     \
    PONCA_FITTING_DECLARE_ADDNEIGHBOR_DER          \
    PONCA_FITTING_DECLARE_FINALIZE

namespace Ponca
{
    /// \FIXME create a macro to automatically generate the testing functions
    template <typename T, typename = void>
    struct hasFirstFundamentalForm : std::false_type
    {
    };

    template <typename T>
    struct hasFirstFundamentalForm<T, std::void_t<decltype(std::declval<T>().firstFundamentalForm())>> : std::true_type
    {
    };

} // namespace Ponca

