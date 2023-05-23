/**
  * @file BasisFunctions.hpp 
  *
  * @brief
  * declares template functions to generate Hermite basis functions 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
#define BASIS_FUNCTIONS_DECLARATIONS 

#ifndef LIBRARIES_LOADER_BF 
  #include "LibrariesLoader_BF.hpp" 
#endif 

/** 
  * @namespace BasisFunctions 
  * 
  * @brief 
  * Contains functions to generate Hermite basis functions 
  * 
  * @anchor _BasisFunctions_ 
  */
namespace BasisFunctions {

template < typename T > 
using Vector = std::vector<T>; 


template < typename Z > 
/**
  * @brief 
  * Generate sets of indices of Hermite polynomials. 
  * Sorted with increasing sum of indices in each set. @n 
  * Implemented in @ref _GenerateIndices_imp_hpp_ 
  *
  * @tparam Z a type of non-negative integer e.g. size_t 
  *
  * @param SetSize number of indices in each set 
  * @param MaxSum  largest allowable sum of indices in each set 
  *
  * @return vector of indices, size of which is a multiple of SetSize 
  */
Vector<Z> GenerateIndices ( const Z SetSize, const Z MaxSum ); 


template < typename Z > 
/**
  * @brief 
  * Remove sets of indices that contain large index. 
  * Operations are performed in-place, no new vector is generated. @n
  * Implemented in @ref _GenerateIndices_imp_hpp_ 
  *
  * @tparam Z a type of non-negative integer e.g. size_t 
  *
  * @param Indices vector of indices e.g. from function GenerateIndices 
  * @param SetSize number of indices in each set 
  * @param iMax    largest allowable index 
  */
void RemoveLargeIndices ( Vector<Z>& Indices, const Z SetSize, const Z iMax );


} // namespace BasisFunctions 

#ifndef GENERATE_INDICES_IMPLEMENTATIONS 
  #include "GenerateIndices_imp.hpp"
#endif 

#endif // BASIS_FUNCTIONS_DECLARATIONS 

