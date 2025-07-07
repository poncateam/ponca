#pragma once

#include "baseType.h"
#include <utility>

#include "differentialQuantities.hpp"

namespace Ponca::Estimators {

    template < typename _Fit >
    class InterfaceProcessHandler {

        protected :
            using FitT = _Fit;
            using WeightFunc = typename FitT::WeightFunction;
            using DataType = typename FitT::DataPoint;

            VectorType               current_pos = VectorType::Zero();
            Ponca::FIT_RESULT        current_fitResult = Ponca::UNDEFINED;
            int                      current_NeighborsCount = 0;

        public:
            virtual ~InterfaceProcessHandler() = default;
            InterfaceProcessHandler() {
                current_pos = VectorType::Zero();
                current_fitResult = Ponca::UNDEFINED;
                current_NeighborsCount = 0;
            };

            // [TODO] Careful with kNN queries
            void initHandler( const DataType& init_data ) {
                current_pos = init_data.pos();
            }

            void initEstimator( FitT &fit, const Scalar radius ) {

                fit.setWeightFunc(WeightFunc(radius));
                fit.init(current_pos);

                current_NeighborsCount = 0;
            }

            virtual void applyNeighbors( FitT & fit, const std::vector<DataType>& neighbors ) = 0;

            Ponca::FIT_RESULT finalize ( FitT &fit ) {
                current_fitResult = fit.finalize();
                if ( current_fitResult == Ponca::STABLE )
                    current_pos = fit.project( current_pos );

                return current_fitResult;
            }

            void applyFunctor( FitT &fit, Quantity<Scalar> &diffQuantity ) {
                diffQuantity.neighbor_count = current_NeighborsCount;
                diffQuantity.non_stable = current_fitResult != Ponca::STABLE;

                if (diffQuantity.non_stable) return;

                diffQuantity.k1 = fit.kmin();
                diffQuantity.k2 = fit.kmax();
                diffQuantity.mean = fit.kMean();
                diffQuantity.gauss = fit.GaussianCurvature();
                diffQuantity.d1 = fit.kminDirection();
                diffQuantity.d2 = fit.kmaxDirection();
                diffQuantity.normal = fit.primitiveGradient();

                diffQuantity.projection = current_pos;
            }

            virtual Scalar potential ( FitT &fit, const VectorType& pos ) = 0;
    }; // InterfaceProcessHandler

    // Generic ProcessHandler
    template < typename _Fit >
    class ProcessHandler : public InterfaceProcessHandler< _Fit > {

        using Base = InterfaceProcessHandler< _Fit >;
        using FitT = typename Base::FitT;
        using DataType = typename FitT::DataPoint;

        public :

            ProcessHandler () : Base() {}
            ~ProcessHandler() = default;

            void applyNeighbors( FitT & fit, const std::vector<DataType>& neighbors ) override {
                Base::current_NeighborsCount = 0;
                fit.startNewPass();

                for (const auto& neighbor : neighbors) {
                    fit.addNeighbor( neighbor );
                    Base::current_NeighborsCount += 1;
                }
            }

            Scalar potential ( FitT &fit, const VectorType& pos ) override {
                return fit.potential( pos );
            }

    }; // ProcessHandler


    template <typename _DataType>
    class BaseEstimator {
        public :
            using DataType = _DataType;

            virtual ~BaseEstimator() = default;
            virtual void setName(std::string name) const = 0;
            [[nodiscard]] virtual std::string getName () const = 0;
            [[nodiscard]] virtual std::string toString() const = 0;
            [[nodiscard]] virtual bool isOriented() const = 0;
            [[nodiscard]] virtual bool isFixedMLS () const = 0;

            virtual void operator() ( const DataType &query,
                                        const std::vector<DataType>& neighbors,
                                        Quantity<Scalar> &quantity ) = 0;

            virtual Scalar potential ( const VectorType& pos ) = 0;

            virtual void setMLSMax( int mls_max ) = 0;

            virtual void setRadius( const Scalar& size ) = 0;

        private :

            Scalar radius = 0;

    }; // BaseEstimator


    // Basket fit wrapper. The Handler store the functors, called in the main loop
    template < typename _FitT, bool _isOriented = false, int _nFixedMLS = -1>
    class Estimator final : public BaseEstimator< typename _FitT::DataPoint > {
        using Self = Estimator<_FitT>;
        using DataType = typename _FitT::DataPoint;

        protected :
            int mls_max = 1;
            float radius = 0;

            bool fixedMLS = false;
            bool oriented = _isOriented;

            _FitT computedFit;
            bool isFinished = false;

            std::string name = "None";
        public :
            // Overload the FitT
            using FitT = _FitT;

            using WeightFunc = typename FitT::WeightFunction;

            Estimator() {
                if (_nFixedMLS != -1) {
                    fixedMLS = true;
                    mls_max = _nFixedMLS;
                }
            }
            Estimator(std::string name) : name(std::move(name)) {
                if (_nFixedMLS != -1) {
                    fixedMLS = true;
                    mls_max = _nFixedMLS;
                }
            }

            void setName (std::string name) const override { this.name = name; }

            [[nodiscard]] std::string getName () const override {  return name; }

            [[nodiscard]] bool isOriented() const override { return oriented; }

            [[nodiscard]] bool isFixedMLS () const override { return fixedMLS; }

            [[nodiscard]] std::string toString() const override { return name; }

            void setMLSMax( int mls_max ) override {
                if ( fixedMLS )
                    throw std::runtime_error("MLS was fixed and can't be changed for this estimator");
                this->mls_max = mls_max;
            }

            void setRadius( const Scalar& size ) override { radius = size; }

            void operator() ( const DataType &query, const std::vector<DataType>& neighbors, Quantity<Scalar> &quantity ) override {

                ProcessHandler<FitT> pHandler( );

                int mls_current = 0;

                pHandler.initHandler ( query );

                for ( mls_current = 0 ; mls_current < mls_max ; mls_current ++ ){
                    FitT fit;
                    pHandler.initEstimator( fit, radius );

                    Ponca::FIT_RESULT res;
                    do {
                        pHandler.applyNeighbors( fit, neighbors );
                        res = pHandler.finalize( fit );
                    } while ( res == Ponca::NEED_OTHER_PASS );

                    if ( mls_current == mls_max - 1 ){
                        pHandler.applyFunctor( fit, quantity );
                        isFinished = true;
                        computedFit = fit;
                    }
                }
            }

            Scalar potential ( const VectorType& pos ) override {
                if ( !isFinished )
                    return 0;
                ProcessHandler<FitT> pHandler( name );
                return pHandler.potential(computedFit, pos);
            }

    }; // Estimator

}