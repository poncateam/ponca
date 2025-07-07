#pragma once

#include "differentialQuantities.hpp"
#include "estimatorFactory.h"

namespace Ponca::Estimators {

    /// Generic processing function: traverse point cloud and compute mean, first and second curvatures + their direction
    /// \tparam FitT Defines the type of estimator used for computation
    template< typename DataType >
    Quantity<Scalar> estimateDifferentialQuantitiesImpl( BaseEstimator<DataType>& estimator, const DataType& query , const std::vector<DataType>& neighbors, const Scalar& radius, const int& nbMLSiter ) {
        if (!estimator)
            throw std::runtime_error("Invalid estimator");

        // Properties
        if (! estimator->isFixedMLS())
            estimator->setMLSMax(nbMLSiter);
        estimator->setRadius(radius);

        // Compute the fitting
        Quantity<Scalar> q;
        try {
            (*estimator)(query, neighbors, q);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error in estimator: " + std::string(e.what()));
        }

        return q;
    }

    template<typename DataType, typename WeightFunc>
    Quantity<Scalar> estimateDifferentialQuantities( FitType name, const DataType& query, const std::vector<DataType>& neighbors, const Scalar& radius, const int& nbMLSiter ) {

        EstimatorFactory< DataType, WeightFunc > factory = getEstimatorFactory<DataType, WeightFunc>();
        BaseEstimator<DataType> estimator = factory->getEstimator(name);

        Quantity<Scalar> quantity = estimateDifferentialQuantitiesImpl<DataType>( estimator, query, neighbors, radius, nbMLSiter );
        return quantity;
    }

    template <typename DataType, typename WeightFunc>
    bool isOriented ( FitType name ) {
        return getEstimatorFactory<DataType, WeightFunc>()->getEstimator(name)->isOriented();
    }
}
