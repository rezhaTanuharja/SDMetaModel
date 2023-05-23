/**
  * @file GenerateIndices_imp.hpp 
  *
  * @brief
  * implements functions to generate indices of Hermite polynomials 
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

#ifndef GENERATE_INDICES_IMPLEMENTATIONS 
#define GENERATE_INDICES_IMPLEMENTATIONS 

#ifndef BASIS_FUNCTIONS_DECLARATIONS 
  #include "BasisFunctions.hpp" 
#endif 

namespace BasisFunctions {


template < typename Z > 
/**
  * @private 
  *
  * @brief  
  * Computes binomial coefficient nCk recursively. @n 
  * Used instead of boost implementation which return floating number. 
  *
  * @tparam Z a type of non-negative integer e.g. size_t 
  *
  * @param n size of set to sample from 
  * @param k size of subset 
  *
  * @return number of distinct k-subset of a set with size n 
  */
Z BinomialCoefficient ( const Z n, const Z k );


template < typename Z >
/**
  * @private 
  *
  * @brief 
  * Generate next set of indices based on the last generated set. 
  *
  * @tparam Z a type of non-negative integer e.g. size_t 
  *
  * @param Indices vector to store all sets of indices 
  * @param LastSet boolean vector with 'true' to mark indices in last set 
  * @param n       number of positive integers used to generate indices 
  *
  */
void GenerateNextIndices ( 
  Vector<Z>& Indices, Vector<bool>& LastSet, const Z n 
);


} // namespace BasisFunctions 


namespace BasisFunctions {


template < typename Z > 
Vector<Z> GenerateIndices ( const Z SetSize, const Z MaxSum ) {

  Vector<Z> Indices; 

  Indices.reserve (
    SetSize * BinomialCoefficient<Z> ( MaxSum + SetSize, SetSize )
  );

  for ( auto i = 0; i < SetSize; i++ ) {
    Indices.push_back ( 0 );
  }

  for ( auto i = 0; i < MaxSum; i++ ) {

    auto jMax = BinomialCoefficient<Z> ( SetSize + i, SetSize - 1 );

    Vector<bool> LastSet ( SetSize + i ); 

    // Initiate marker s.t. it's only one permuation away from the first set 
    
    std::transform (
      LastSet.begin(), LastSet.begin() + SetSize - 1, 
      LastSet.begin(), 

      [=]( auto m ) { return true; }
    );

    // Generate indices and store in Indices vector 

    for ( auto j = 0; j < jMax; j++ ) {
      GenerateNextIndices<Z> ( Indices, LastSet, SetSize + i );
    }

  }

  return Indices; 

} // GenerateIndices 


template < typename Z > 
void RemoveLargeIndices ( Vector<Z>& Indices, const Z SetSize, const Z iMax ) {

  auto nSet = Indices.size(); 

  auto i = 0; 

  while ( true ) {

    if ( i >= nSet ) { break; }

    auto MaxIt = std::max_element (
      Indices.begin() + i, 
      Indices.begin() + i + SetSize 
    );

    if ( *MaxIt <= iMax ) { 
      i += SetSize; 
      continue; 
    }

    Indices.erase (
      Indices.begin() + i, 
      Indices.begin() + i + SetSize
    );

    nSet -= SetSize; 

  }

} // RemoveLargeIndices 


} // namespace BasisFunctions 


namespace BasisFunctions {


template < typename Z > 
Z BinomialCoefficient ( const Z n, const Z k ) {

  if ( n < k ) { return 0; } 

  if ( k == 0 || k == n ) { return 1; } 

  return BinomialCoefficient<Z> ( n - 1, k - 1 ) +
         BinomialCoefficient<Z> ( n - 1, k     );

} // BinomialCoefficient 


template < typename Z >
void GenerateNextIndices ( 
  Vector<Z>& Indices, Vector<bool>& LastSet, const Z n 
) {

  Vector<Z> NextIndices; 

  Vector<Z> Integers ( n ); 
  std::iota ( Integers.begin(), Integers.end(), 1 ); 

  std::next_permutation ( LastSet.begin() , LastSet.end() ); 

  Z m = 0; 
  for ( const auto& Active : LastSet ) {

    if ( Active ) {
      NextIndices.push_back ( Integers[m] );
    }

    m++; 

  }

  NextIndices.push_back ( n + 1 ); 

  std::adjacent_difference (
    NextIndices.begin(), NextIndices.end(), 
    NextIndices.begin()
  );

  std::transform (
    NextIndices.begin(), NextIndices.end(), 
    NextIndices.begin(), 

    [=]( auto m ) { return m - 1; }
  );

  Indices.insert( Indices.end(), NextIndices.begin(), NextIndices.end() );

} // GenerateNextIndices 


} // namespace BasisFunctions 


#endif // GENERATE_INDICES_IMPLEMENTATIONS 

