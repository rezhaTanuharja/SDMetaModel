/**
  * @file Normal.hpp 
  *
  * @brief
  * declares template functions to generate RVs from normal distributions 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef NORMAL_DECLARATIONS 
#define NORMAL_DECLARATIONS 

#ifndef STATISTICS_DECLARATIONS 
  #include "Statistics.hpp" 
#endif

namespace Statistics {

/**
  * @namespace Normal 
  * Contain functions to generate RVs from normal distributions 
  *
  * @anchor _Normal_ 
  */
namespace Normal {


template < typename Z, typename R, typename C >
/**
  * @brief 
  * Perform random sampling from a @ref _Normal_ distribution. @n 
  * Implemented in @ref _Normal_imp_hpp_ 
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
  const R Mean, const R StdDev
);


} // namespace Normal 

} // namespace Statistics 

#ifndef NORMAL_IMPLEMENTATIONS 
  #include "Normal_imp.hpp" 
#endif 

#endif // NORMAL_DECLARATIONS 

