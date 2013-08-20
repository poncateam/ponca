/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/

#include <Eigen/Eigenvalues>

#ifndef _GRENAILLE_GLS_
#define _GRENAILLE_GLS_

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
    
    /*! 
      \brief Compare current instance with other 
    */
    MULTIARCH inline Scalar compareTo (const GLSParam<DataPoint, _WFunctor, T>& other,
                                       bool useFitness = true) const
    { 
      Scalar nTau     = this->tau_normalized()   - other.tau_normalized();
      Scalar nKappa   = this->kappa_normalized() - other.kappa_normalized();
      Scalar nFitness = useFitness ? this->fitness() - other.fitness() : 0.;
      
      return nTau * nTau + nKappa * nKappa + nFitness * nFitness;
    }

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
    \brief Extension to compute curvature values from \f$ \frac{\delta\tau}{\tau\mathbf{x}} \f$
    \inherit Concept::FittingExtensionConcept
    
    This class requires an Eigendecomposition of the Jacobian matrix of 
    \f$ \tau \f$. We use an Eigen::SelfAdjointEigenSolver, compatible with 
    Eigen-nvcc for 2D and 3D spaces on cuda, and in arbitrary dimension for CPU
    compilation (cuda version also use a fast but less accurate decomposition
    provided by computeDirect. See <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html#a85cda7e77edf4923f3fc0512c83f6323" target="_blank">Eigen documentation</a> for mor details).
   */
  template < class DataPoint, class _WFunctor, typename T>
  class GLSCurvatureHelper : public T{
  private:
    typedef T Base;

  protected:
    enum
      {
        Check = Base::PROVIDES_ALGEBRAIC_SPHERE_SPACE_DERIVATIVE & Base::PROVIDES_GLS_DERIVATIVE,
	      PROVIDES_PRINCIPALE_CURVATURES
      };

  public:
    typedef typename Base::Scalar     Scalar;
    typedef typename Base::VectorType VectorType;
    typedef typename Base::MatrixType MatrixType;
      
  private:
    Scalar _k1, _k2;
    VectorType _v1, _v2;

    
  public:
    /*! \brief Default constructor */
    MULTIARCH inline GLSCurvatureHelper() : _k1(0), _k2(0) {}
    
    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    MULTIARCH inline void finalize   (){
	    MULTIARCH_STD_MATH(sqrt);
	      
      Base::finalize();
      
      // Extract the spatial variations of eta
      MatrixType jacobian = Base::deta().template middleCols<DataPoint::Dim>(Base::isScaleDer() ? 1: 0);
      
      // Use a simple solver with 2x2 and 3x3 closed forms compatible with eigen-nvcc
      // This solver requires a triangular matrix
      Eigen::SelfAdjointEigenSolver<MatrixType> eig;
#ifdef __CUDACC__
      eig.computeDirect(jacobian.transpose()*jacobian);
#else
      eig.compute(jacobian.transpose()*jacobian);
#endif
      
      // Need sqrt because we solve eigendecomposition of JT * J.
      _k1 = sqrt(eig.eigenvalues()(2)); 
      _k2 = sqrt(eig.eigenvalues()(1)); 
      
      _v1 = eig.eigenvectors().col(2);
      _v2 = eig.eigenvectors().col(1);
      
      // Now check the sign of the mean curvature to detect if we need to change
      // the sign of the principal curvature values:
      // the eigen decomposition return k1 and k2 wrt k1*k1 > k2*k2
      
      // Compare with the values of the mean curvature computed without k1 and k2
      Scalar H2 = Scalar(2)*Base::kappa(); // we could also use the trace of the
                                           // jacobian matrix to get mean curvature
      
      // Change sign and flip values if needed
      if( H2 == Scalar(0)){
        _k2 = -_k2;        
        
      } else if ( H2 > Scalar(0) ) {
        if( H2 < _k1 )
          _k2 = -_k2;
      }else { // 2H < 0. In this case, we have k1<0, and k1 < k2
        if( H2 > _k1 )
          _k2 = -_k2;
        _k1 = -_k1;
        
        // need to invert k1 and k2, and get the corresponding vectors
        Scalar tmp = _k1; _k1 = _k2; _k2 = tmp;        
        _v1 = eig.eigenvectors().col(1);
        _v2 = eig.eigenvectors().col(2);        
      }
    }
    
    /**************************************************************************/
    /* Use results                                                            */
    /**************************************************************************/
    //! \brief Return the first principal curvature value
    MULTIARCH inline Scalar GLSk1() const { return _k1; }
    
    //! \brief Return the second principal curvature value
    MULTIARCH inline Scalar GLSk2() const { return _k2; }
    
    //! \brief Return the first principal curvature direction
    MULTIARCH inline VectorType GLSk1Direction() const { return _v1; }
    
    //! \brief Return the second principal curvature direction
    MULTIARCH inline VectorType GLSk2Direction() const { return _v2; }
    
    //! \brief Return the Gaussian Curvature
    MULTIARCH inline Scalar GLSGaussianCurvature() const { return _k1*_k2;}    
  };
  
  #include "gls.hpp"

} //namespace Grenaille

#endif
