#pragma once

#include <memory>

#include "defineEnum.h"
#include "defineEstimators.h"
#include "../basket.h"


namespace Ponca::Estimators {

//     template <typename DataType, typename WeightFunc>
//     class EstimatorFactory {
//
//         std::shared_ptr< Ponca::BasketBase< DataType, WeightFunc > > estimatorsList[NUMBER_OF_FIT_TYPES];
//
//     public:
//
//         EstimatorFactory() {
// #define ENUM_FIT(name) \
//             estimatorsList[Estimators::FitType::name] = std::make_shared<Estimators::Fit ## _ ## name<WeightFunc>>(#name);
// ENUM_FITS
// #undef ENUM_FIT
//         }
//
//         std::shared_ptr< Ponca::BasketBase< DataType, WeightFunc > > getEstimator(FitType name) {
//             if (estimatorsList[name] == nullptr)
//                 throw std::runtime_error("Estimator type not found");
//             return estimatorsList[name];
//
//         }
//     };
//
//
//     template <typename DataType, typename WeightFunc>
//     std::shared_ptr<EstimatorFactory< DataType, WeightFunc > > getEstimatorFactory() {
//         static auto factory = std::make_shared<EstimatorFactory< DataType, WeightFunc >>();
//         return factory;
//     }

}

namespace Ponca
{

    template <typename DataType, typename WeightFunc>
    std::shared_ptr<BasketBase< DataType, WeightFunc > > getFit(const Estimators::FitType name) {
        switch (name) {
#define ENUM_FIT(name) \
case Estimators::FitType::name :   \
return std::make_shared<Estimators::Fit ## _ ## name<WeightFunc>>(#name); \
break;
            ENUM_FITS
            #undef ENUM_FIT
        }
        throw std::runtime_error("Unknown fit type");
    }
}