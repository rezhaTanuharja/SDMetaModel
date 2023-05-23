#include "BasisFunctions.hpp" 
#include <gtest/gtest.h> 


TEST ( GenerateIndices, Linear1D ) {

  std::vector<size_t> expected {
    0, 1,
  };

  auto result = BasisFunctions::GenerateIndices<size_t> ( 1, 1 );

  ASSERT_EQ ( result.size(), expected.size() );

  for ( auto i = 0; i < result.size(); i++ ) {

    EXPECT_EQ ( result[i], expected[i] );

  }

} // GenerateIndices, Linear1D


TEST ( GenerateIndices, Linear2D ) {

  std::vector<size_t> expected {
    0, 0,
    1, 0,
    0, 1,
  };

  auto result = BasisFunctions::GenerateIndices<size_t> ( 2, 1 );

  ASSERT_EQ ( result.size(), expected.size() );

  for ( auto i = 0; i < result.size(); i++ ) {

    EXPECT_EQ ( result[i], expected[i] );

  }

} // GenerateIndices, Linear2D


TEST ( GenerateIndices, Linear3D ) {

  std::vector<size_t> expected {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  };

  auto result = BasisFunctions::GenerateIndices<size_t> ( 3, 1 );

  ASSERT_EQ ( result.size(), expected.size() );

  for ( auto i = 0; i < result.size(); i++ ) {

    EXPECT_EQ ( result[i], expected[i] );

  }

} // GenerateIndices, Linear3D


TEST ( GenerateIndices, Quadratic1D ) {

  std::vector<size_t> expected {
    0, 1, 2
  };

  auto result = BasisFunctions::GenerateIndices<size_t> ( 1, 2 );

  ASSERT_EQ ( result.size(), expected.size() );

  for ( auto i = 0; i < result.size(); i++ ) {

    EXPECT_EQ ( result[i], expected[i] );

  }

} // GenerateIndices, Quadratic1D


TEST ( GenerateIndices, Quadratic2D ) {

  std::vector<size_t> expected {
    0, 0,
    1, 0,
    0, 1,
    2, 0,
    1, 1,
    0, 2
  };

  auto result = BasisFunctions::GenerateIndices<size_t> ( 2, 2 );

  ASSERT_EQ ( result.size(), expected.size() );

  for ( auto i = 0; i < result.size(); i++ ) {

    EXPECT_EQ ( result[i], expected[i] );

  }

} // GenerateIndices, Quadratic2D


TEST ( GenerateIndices, Quadratic3D ) {

  std::vector<size_t> expected {
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

  auto result = BasisFunctions::GenerateIndices<size_t> ( 3, 2 );

  ASSERT_EQ ( result.size(), expected.size() );

  for ( auto i = 0; i < result.size(); i++ ) {

    EXPECT_EQ ( result[i], expected[i] );

  }

} // GenerateIndices, Quadratic3D

