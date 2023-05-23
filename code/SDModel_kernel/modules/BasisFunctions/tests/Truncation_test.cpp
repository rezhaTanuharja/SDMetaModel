#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 


TEST ( RemoveLargeIndices, NoTruncation ) {

  std::vector<size_t> Indices {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  };

  std::vector<size_t> expected {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  };

  size_t SetSize = 3; 
  size_t iMax    = 1; 

  BasisFunctions::RemoveLargeIndices<size_t> ( Indices, SetSize, iMax );

  ASSERT_EQ ( Indices.size(), expected.size() );

  for ( auto i = 0; i < Indices.size(); i++ ) {

    EXPECT_EQ ( Indices[i], expected[i] );

  }

} // RemoveLargeIndices, NoTruncation 


TEST ( RemoveLargeIndices, MaxConstant ) {

  std::vector<size_t> Indices {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  };

  std::vector<size_t> expected {
    0, 0, 0
  };

  size_t SetSize = 3; 
  size_t iMax    = 0; 

  BasisFunctions::RemoveLargeIndices<size_t> ( Indices, SetSize, iMax );

  ASSERT_EQ ( Indices.size(), expected.size() );

  for ( auto i = 0; i < Indices.size(); i++ ) {

    EXPECT_EQ ( Indices[i], expected[i] );

  }

} // RemoveLargeIndices, MaxConstant  


TEST ( RemoveLargeIndices, MaxLinear ) {

  std::vector<size_t> Indices {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    2, 0, 0,
    1, 1, 0,
    1, 0, 1,
    0, 2, 0, 
    0, 1, 1,
    0, 0, 2
  };

  std::vector<size_t> expected {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    1, 1, 0,
    1, 0, 1,
    0, 1, 1
  };

  size_t SetSize = 3; 
  size_t iMax    = 1; 

  BasisFunctions::RemoveLargeIndices<size_t> ( Indices, SetSize, iMax );

  ASSERT_EQ ( Indices.size(), expected.size() );

  for ( auto i = 0; i < Indices.size(); i++ ) {

    EXPECT_EQ ( Indices[i], expected[i] );

  }

} // RemoveLargeIndices, MaxLinear  

