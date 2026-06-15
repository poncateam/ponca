/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "../../../Fitting"
#include "typelist.h"

namespace Ponca
{

// List of concept that can be filtered by the factory
// We can not template on concept. Here, we trick this
// by creating a class that extract the boolean from the
// given concept
#define DECLARE_FACTORY_CONCEPT_FULL(_conceptname, _newname)             \
    template <typename T>                                                \
    struct _newname : std::bool_constant<_conceptname<typename T::type>> \
    {                                                                    \
    };
#define DECLARE_FACTORY_CONCEPT(_name) DECLARE_FACTORY_CONCEPT_FULL(Provides##_name, _name##Provider)

    /**
     * \brief Class to check for if an entry has a corresponding Id
     *
     * Instead of adding more code specialized for this, we create
     * a predicated that indicates if any factory entry matches the
     * given id
     *
     * \tparam Id The Id to check for
     */
    template <int Id>
    struct MethodProvider
    {
        template <typename T>
        static constexpr bool value = (Id == T::Id);

        /**
         * \brief A type that behaves as value
         *
         * This will be used in conjunction which requires
         * templated types on the type.
         */
        template <typename T>
        struct pred : std::bool_constant<(Id == T::Id)>
        {
        };
    };

    DECLARE_FACTORY_CONCEPT(ProjectionOperator)
    DECLARE_FACTORY_CONCEPT(ImplicitPrimitive)
    DECLARE_FACTORY_CONCEPT(AlgebraicSphere)
    DECLARE_FACTORY_CONCEPT(MeanPosition)
    DECLARE_FACTORY_CONCEPT(MeanNormal)
    DECLARE_FACTORY_CONCEPT(GLSParam)
    DECLARE_FACTORY_CONCEPT(GeomVar)
    DECLARE_FACTORY_CONCEPT(Plane)
    DECLARE_FACTORY_CONCEPT(Line)
    DECLARE_FACTORY_CONCEPT(HeightField)
    DECLARE_FACTORY_CONCEPT(MongePatch)
    DECLARE_FACTORY_CONCEPT(PositionCovariance)
    DECLARE_FACTORY_CONCEPT(PrincipalCurvatures)
    DECLARE_FACTORY_CONCEPT(TangentPlaneBasis)
    DECLARE_FACTORY_CONCEPT(MeanCurvature)

    /**
     * List of methods supported by the factory
     */
    enum class Method
    {
        UNNAMED = 0, // 0 is considered as an unnamed, and is the default
        APSS,
    };

    // We disable format here in order to better align lists items
    // clang-format off
    template<typename P, typename NF, int DerType>
    using FactoryGenericList = std::tuple<
        FactoryEntry<"CovariancePlaneFit", Basket<P, NF, CovariancePlaneFit>>, 
        FactoryEntry<"SphereFit"         , Basket<P, NF, SphereFit>>, 
        FactoryEntry<"MeanPlaneFit"      , Basket<P, NF, MeanPlaneFit>>, 
        FactoryEntry<"APSS"              , Basket<P, NF, OrientedSphereFit, GLSParam>, (unsigned int)Method::APSS>
    >;

    template<typename P, typename NF, int DerType>
    using Factory3DList = std::tuple<
        FactoryEntry<"MongePatchQuadratic"          , Basket<P, NF, MongePatchQuadraticFit>>, 
        FactoryEntry<"MongePatchRestrictedQuadratic", Basket<P, NF, MongePatchRestrictedQuadraticFit>>
    >;
    // clang-format on
} // namespace Ponca
