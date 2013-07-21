/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/


#ifndef _GRENAILLE_GLS_
#define _GRENAILLE_GLS_


#include <Eigen/Eigenvalues>
#include <utility>


namespace Grenaille
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
          
    Requierement: 
    \verbatim PROVIDES_ALGEBRAIC_SPHERE \endverbatim
    Provide: 
    \verbatim PROVIDES_GLS_PARAMETRIZATION \endverbatim
   */
  template < class DataPoint, class _WFunctor, typename T>
  class GLSParam : public T{
  private:
    typedef T Base;

  protected:
    enum
      {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE,
        PROVIDES_GLS_PARAMETRIZATION
      };

  public:
    typedef typename Base::Scalar     Scalar;     /*!< \brief Inherited scalar type*/
    typedef typename Base::VectorType VectorType; /*!< \brief Inherited vector type*/
    typedef typename Base::WFunctor   WFunctor;   /*!< \brief Weight Function*/


  protected:
    Scalar _t; /*!< \brief Evaluation scale. Needed to computed the normalized values*/
    Scalar _fitness; /*!< \brief Save the fitness value to avoid side effect with Pratt normalization*/

  public:
    /*! \brief Default constructor */
    MULTIARCH inline GLSParam() : _t(0) {}

    
    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::setWeightFunc() */
    MULTIARCH inline void setWeightFunc (const WFunctor& w){
      Base::setWeightFunc(w);
      _t = w.evalScale();
    }
    
    
    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    MULTIARCH inline void finalize   (){
      Base::finalize();
      _fitness = Scalar(1.) - Base::prattNorm2();
    }


    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
    /*! \brief Compute and return \f$ \tau \f$ */
    MULTIARCH inline Scalar     tau()   const 
    {return Base::isNormalized() ? Base::_uc : Base::_uc / Base::prattNorm();}

    /*! \brief Compute and return \f$ \eta \f$ */
    MULTIARCH inline VectorType eta()   const 
    {return Base::_ul / Base::_ul.norm();}

    /*! \brief Compute and return \f$ \kappa \f$ */
    MULTIARCH inline Scalar     kappa() const 
    {return Scalar(2.) * (Base::isNormalized() ? Base::_uq : Base::_uq / Base::prattNorm());}
    
    /*! \brief Compute and return \f$ \frac{\tau}{t} \f$ */
    MULTIARCH inline Scalar     tau_normalized()   const {return tau()/_t;}
    
    /*! \brief Compute and return \f$ \eta \f$ */    
    MULTIARCH inline VectorType eta_normalized()   const {return eta();}
    
    /*! \brief Compute and return \f$ t \kappa \f$ */
    MULTIARCH inline Scalar     kappa_normalized() const {return kappa()*_t;}    
    
    /*! \brief Return the fitness, e.g. the pratt norm of the initial scalar field */
    MULTIARCH inline Scalar     fitness()          const {return _fitness;}

  }; //class GLSParam

  
  /*!
    \brief Differentiation of GLSParam
    \inherit Concept::FittingExtensionConcept

    Method published in \cite Mellado:2012:GLS
   */
  template < class DataPoint, class _WFunctor, typename T>
  class GLSDer : public T{
  private:
    typedef T Base;

  protected:
    enum
      {
        Check = Base::PROVIDES_GLS_PARAMETRIZATION,
        PROVIDES_GLS_DERIVATIVE
      };

  public:
    typedef typename Base::Scalar     Scalar;
    typedef typename Base::VectorType VectorType;
    typedef typename Base::WFunctor   WFunctor;

    typedef typename Base::VectorArray VectorArray;
    typedef typename Base::ScalarArray ScalarArray;   
    
    MULTIARCH inline ScalarArray dtau()   const; /*!< \brief Access to \f$ \tau \f$ derivatives */
    MULTIARCH inline VectorArray deta()   const; /*!< \brief Access to \f$ \eta \f$ derivatives */
    MULTIARCH inline ScalarArray dkappa() const; /*!< \brief Access to \f$ \kappa \f$ derivatives */
    
    MULTIARCH inline ScalarArray dtau_normalized()   const;
    MULTIARCH inline VectorArray deta_normalized()   const;
    MULTIARCH inline ScalarArray dkappa_normalized() const;
  }; //class GLSScaleDer


  /*!
    \brief Extension to compute the Geometric Variation of GLSParam
    \inherit Concept::FittingExtensionConcept
    
    Method published in \cite Mellado:2012:GLS
    \todo Add more details
   */
  template < class DataPoint, class _WFunctor, typename T>
  class GLSGeomVar : public T{
  private:
    typedef T Base;

  protected:
    enum
      {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE_SCALE_DERIVATIVE & Base::PROVIDES_GLS_DERIVATIVE,
	      PROVIDES_GLS_GEOM_VAR
      };

  public:
    typedef typename Base::Scalar Scalar;

    /*!
      \brief Compute and return the Geometric Variation
     */
    MULTIARCH inline Scalar geomVar(Scalar wtau   = Scalar(1), 
  			                            Scalar weta   = Scalar(1),
	                        			    Scalar wkappa = Scalar(1)) const;
  };
  
  
  /*!
    \brief Extension to compute the variation tensor of GLSParam
    \inherit Concept::FittingExtensionConcept
    
    Method published in \cite Mellado:2013:Thesis
    
    \warning This class is valid only in dimension 3.
    
    \todo Add more details
   */
  template < class DataPoint, class _WFunctor, typename T>
  class GLSSpatialVariation : public T{
  private:
    typedef T Base;

  protected:
    enum
      {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE_SPACE_DERIVATIVE & Base::PROVIDES_GLS_DERIVATIVE,
	      PROVIDES_GLS_SPACE_VARIATION
      };

  public:
    typedef typename Base::Scalar Scalar;
    typedef typename Base::VectorType VectorType;
    typedef typename Base::MatrixType MatrixType;
    // \todo Change this by using cuda compatible structure (or check compatibility)
    typedef typename std::pair< Eigen::Matrix<Scalar, 1, DataPoint::Dim-1>,
                                Eigen::Matrix<Scalar, DataPoint::Dim, DataPoint::Dim-1> > GLSSpatialEigenDecomposition;

    /*!
      \brief Compute and return the eigendecomposition of the spatial variation tensor projected on the tangent plane
      
      The tangent plane is defined by \f$\eta\f$.
      
      \warning This method is valid only in dimension 3.
      \todo Add a way to generate an error at compilation time for dimension check
      
      \param wtau Influence of \f$ \frac{\delta\tau}{\delta\mathbf{x}} \f$ in variation analysis. Set to 0 to ignore.
      \param weta Influence of \f$ \frac{\delta\eta}{\delta\mathbf{x}} \f$ in variation analysis. Set to 0 to ignore.
      \param wkappa Influence of \f$ \frac{\delta\kappa}{\delta\mathbf{x}} \f$ in variation analysis. Set to 0 to ignore.
      
      Use weighting parameters to measure different properties:
      \code
      // return the variation of the complete sphere
      fit.projectedVariationDecomposition();
      
      // return the variation of tau
      fit.projectedVariationDecomposition(1,0,0);
      
      // return the variation of eta, that is an approximation of the curvature tensor
      fit.projectedVariationDecomposition(0,1,0);
      
      // return the variation of kappa and get the mean curvature tensor
      fit.projectedVariationDecomposition(0,0,1);      
      \endcode
      
      \return Eigen values and vectors describing the spatial variations, stored
      in a pair<values,vectors>. Since the decomposition is computed in the 
      tangent plane, this function returns Dim-1 values and vectors:
      \code
      MyFit fit;
      
      // ... Compute fit
      
      // compute the decomposition
      MyFit::GLSSpatialEigenDecomposition decomp = fit.projectedVariationDecomposition();
      
      // Print eigen values and associated eigen vectors
      cout << "Eigen values:\n\t"  << decomp.first << endl;
      cout << "Eigen vectors:\n\t" << decomp.second << endl;
      \endcode
     */
    MULTIARCH inline GLSSpatialEigenDecomposition
    projectedVariationDecomposition( Scalar wtau   = Scalar(1), 
  			                             Scalar weta   = Scalar(1),
	                        			     Scalar wkappa = Scalar(1)) const;
  };

  #include "gls.hpp"

} //namespace Grenaille

#endif
