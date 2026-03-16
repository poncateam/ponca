// This file is expected to be compiled by doxygen only, for documentation purpose

namespace Ponca 
{
    namespace Concept 
    {
//! [PointConcept]
        class PointConcept 
        {
        public:
            /* \brief Defines the ambient space dimension, 3 in this example */
            enum { Dim = 3 };
            // \brief Defines the type used ton encode scalar values
            typedef float Scalar;
            // \brief Defines type used ton encode vector values
            // 
            // \note VectorType should have the same API than Eigen::Matrix API, 
            // and be implicitly convertible to Eigen::Matrix
            typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;

            // \brief Default constructor
            PONCA_MULTIARCH inline PointConcept(const VectorType& pos = VectorType::Zero());

            // \brief Read access to the position property
            // \note The return type should have the same API than VectorType (e.g. const Eigen::Map< const VectorType >&)
            // and be implicitly convertible to Eigen::Matrix
            PONCA_MULTIARCH inline const VectorType& pos() const;

            // \brief Write access to the position property
            // \note Same constraints on return type
            PONCA_MULTIARCH inline VectorType& pos();

        }; //class PointConcept
//! [PointConcept]

//! [ComputationalObjectConcept]
        // Naming of template parameters matters for macros 
        template < class DataPoint, class _NFilter, typename T>
        class ComputationalObjectConcept : public T
        {
            // Defines common types (Base, Scalar, VectorType, NeighborFilter)
            PONCA_FITTING_DECLARE_DEFAULT_TYPES
        protected:
            enum {
                Check = PROVIDES_CAPABILITY_1 && PROVIDES_CAPABILITY_2 //< List of required capabilities
                PROVIDES_CAPABILITY,                                   //< List of provided capabilities
            };
        public:
            /**************************************************************************/
            /* Initialization                                                         */
            /**************************************************************************/
            // Init the WeightFunc, without changing the other internal states.
            // \warning Must be called before any computation
            PONCA_MULTIARCH inline void setNeighborFilter (const NeighborFilter& _w);

            // Reset the internal states.
            // \warning Must be called before any computation
            PONCA_MULTIARCH inline void init();

            /**************************************************************************/
            /* Computations                                                           */
            /**************************************************************************/
            // Add a neighbor to perform the fit
            // \param w The weight of the point
            // \param lQ Local coordinates of the point (e.g.: the difference between the point and the eval position)
            // \param pt The point as defined by the user
            // \return false if param nei is not a valid neighbour (weight = 0)
            PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &lQ, const DataPoint & pt);
            // Finalize the fitting procedure. This function is called after all neighbors have been processed.
            // \return State of fitting
            // \warning Must be called be for any use of the fitting output
            PONCA_MULTIARCH inline FIT_RESULT finalize();
        }; //class ComputationalObjectConcept
//! [ComputationalObjectConcept]

//! [ComputationalDerivativesConcept]
        // Naming of template parameters matters for macros 
        template < class DataPoint, class _NFilter, int Type, typename T>
        class ComputationalDerivativesConcept
        {
            // Defines common types (Base, Scalar, VectorType, NeighborFilter)
            PONCA_FITTING_DECLARE_DEFAULT_TYPES
            PONCA_FITTING_DECLARE_DEFAULT_DER_TYPES
        protected:
            enum {
                Check = PROVIDES_CAPABILITY_1 && PROVIDES_CAPABILITY_2 //< List of required capabilities
                PROVIDES_CAPABILITY,                                   //< List of provided capabilities
            };
        public:
            /**************************************************************************/
            /* Initialization                                                         */
            /**************************************************************************/
            // Init the WeightFunc, without changing the other internal states.
            // \warning Must be called be for any computation
            PONCA_MULTIARCH inline void setNeighborFilter (const NeighborFilter& _w);

            // Set the evaluation position and reset the internal states.
            // \warning Must be called be for any computation
            PONCA_MULTIARCH inline void init();

            /**************************************************************************/
            /* Computations                                                           */
            /**************************************************************************/
            // Add a neighbor to perform the fit
            // \param w The weight of the point
            // \param localQ Local coordinates of the point (e.g.: the difference between the point and the eval position)
            // \param attributes The point as defined by the user
            // \param dw The derivatives of the weight function w.r.t prescirbed parameter
            // \return false if param nei is not a valid neighbour (weight = 0)
            PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w,
                                                         const VectorType &localQ,
                                                         const DataPoint &attributes,
                                                         ScalarArray &dw);
            // Finalize the fitting procedure. This function is called after all neighbors have been processed.
            // \return State of fitting
            // \warning Must be called be for any use of the fitting output
            PONCA_MULTIARCH inline FIT_RESULT finalize();
        }; //class ComputationalDerivativesConcept
//! [ComputationalDerivativesConcept]

//! [WeightKernelConcept]
        // \brief Concept of a 1D weighting function and its derivatives.
        template <typename _Scalar>
        class WeightKernelConcept{
        public:
            typedef _Scalar Scalar;

            // \brief Apply the weighting kernel to the scalar value f(x)
            PONCA_MULTIARCH inline Scalar f  (const Scalar& x) const {}
            // \brief Apply the first derivative of the weighting kernel to the scalar value f'(x)
            PONCA_MULTIARCH inline Scalar df (const Scalar& x) const {}
            // \brief Apply the second derivative of the weighting kernel to the scalar value f''(x)
            PONCA_MULTIARCH inline Scalar ddf(const Scalar& x) const {}
        };// class WeightKernelConcept
//! [WeightKernelConcept]

    } //namespace Concept
} //namespace Ponca
