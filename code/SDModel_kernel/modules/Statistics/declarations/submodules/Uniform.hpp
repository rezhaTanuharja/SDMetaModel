/**
  * @file Uniform.hpp 
  *
  * @brief
  * declares template functions to generate RVs from uniform distributions 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef UNIFORM_DECLARATIONS 
#define UNIFORM_DECLARATIONS 

#ifndef STATISTICS_DECLARATIONS 
  #include "Statistics.hpp" 
#endif

namespace Statistics {

/**
  * @namespace Uniform 
  * Contain functions to generate RVs from uniform distributions 
  *
  * @anchor _Uniform_ 
  */
namespace Uniform {


template < typename Z, typename R, typename C >
/**
  * @brief 
  * Perform random sampling from a @ref _Uniform_ distribution. @n 
  * Implemented in @ref _Uniform_imp_hpp_ 
  *
  * @tparam Z non-negative integer e.g. size_t 
  * @tparam R float number e.g. double 
  * @tparam C complex float number e.g. std::complex<double>  
  *
  * @param nPoints    number of sampled points 
  * @param Dim        dimension of each point
  * @param LowerBound lowest possible random variable value 
  * @param UpperBound higher possible random variable value 
  */
Vector<C> RandomSampling (
  const Z nPoints, const Z Dim, 
  const R LowerBound, const R UpperBound
);


} // namespace Uniform 

} // namespace Statistics 

#ifndef UNIFORM_IMPLEMENTATIONS 
  #include "Uniform_imp.hpp" 
#endif

#endif // UNIFORM_DECLARATIONS 

