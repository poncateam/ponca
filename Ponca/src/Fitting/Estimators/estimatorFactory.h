#pragma once

#include <memory>

#include "defineEnum.h"
#include "defineEstimators.h"
#include "estimator.h"


namespace Ponca::Estimators {

    template <typename DataType, typename WeightFunc>
    class EstimatorFactory {

        std::shared_ptr< BaseEstimator< DataType > > estimatorsList[NUMBER_OF_FIT_TYPES];

    public:


        EstimatorFactory() {
#define ENUM_FIT(name) \
            estimatorsList[Estimators::FitType::name] = std::make_shared<Estimators::Estimator ## _ ## name<WeightFunc>>(#name);
ENUM_FITS
#undef ENUM_FIT
        }

        std::shared_ptr< BaseEstimator<DataType> > getEstimator(FitType name) {
            if (estimatorsList[name] == nullptr)
                throw std::runtime_error("Estimator type not found");
            return estimatorsList[name];

        }
    };


    template <typename DataType, typename WeightFunc>
    std::shared_ptr<EstimatorFactory< DataType, WeightFunc > > getEstimatorFactory() {
        static auto factory = std::make_shared<EstimatorFactory< DataType, WeightFunc >>();
        return factory;
    }
}