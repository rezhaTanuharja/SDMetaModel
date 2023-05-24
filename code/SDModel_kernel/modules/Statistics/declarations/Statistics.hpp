/**
  * @file Statistics.hpp 
  *
  * @brief
  * declares template functions for statistics and samplings 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef STATISTICS_DECLARATIONS 
#define STATISTICS_DECLARATIONS 

#ifndef LIBRARIES_LOADER_ST 
  #include "LibrariesLoader_ST.hpp" 
#endif 

/**
  * @namespace Statistics 
  * Contain functions for statistics and samplings 
  *
  * @anchor _Statistics_ 
  */
namespace Statistics {

template < typename T > 
using Vector = std::vector<T>; 

template < typename T > 
using MatrixXT = Eigen::Matrix < T, Eigen::Dynamic, Eigen::Dynamic >;

template < typename T > 
using VectorXT = Eigen::Vector < T, Eigen::Dynamic >;

} // Statistics 

#ifndef UNIFORM_DECLARATIONS 
  #include "Uniform.hpp" 
#endif 

#ifndef NORMAL_DECLARATIONS 
  #include "Normal.hpp" 
#endif 

#endif // STATISTICS_DECLARATIONS 

