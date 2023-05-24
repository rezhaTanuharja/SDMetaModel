/**
  * @file GenerateIndices_test.cpp 
  *
  * @brief
  * test functions to generate polynomial indices, 
  * including errors thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

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


#ifdef DEBUG_BASIS_FUNCTIONS 


TEST ( GenerateIndices, NegativeSetSize ) {

  EXPECT_THROW ({

    try {

      auto result = BasisFunctions::GenerateIndices<size_t> ( -1, 1 );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "GenerateIndices: SetSize not in [1,100] range", e.what()
        
      );
      
      throw; 
      
    }

  }, std::runtime_error );

} // GenerateIndices, NegativeSetSize 


TEST ( GenerateIndices, ZeroSetSize ) {

  EXPECT_THROW ({

    try {

      auto result = BasisFunctions::GenerateIndices<size_t> ( 0, 1 );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "GenerateIndices: SetSize not in [1,100] range", e.what()
        
      );
      
      throw; 
      
    }

  }, std::runtime_error );

} // GenerateIndices, ZeroSetSize  


TEST ( GenerateIndices, HugeSetSize ) {

  EXPECT_THROW ({

    try {

      auto result = BasisFunctions::GenerateIndices<size_t> ( 1000, 1 );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "GenerateIndices: SetSize not in [1,100] range", e.what()
        
      );
      
      throw; 
      
    }

  }, std::runtime_error );

} // GenerateIndices, HugeSetSize   


TEST ( GenerateIndices, NegativeMaxSum ) {

  EXPECT_THROW ({

    try {

      auto result = BasisFunctions::GenerateIndices<size_t> ( 5, -1 );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "GenerateIndices: MaxSum not in [0,100] range", e.what() 
        
      );
      
      throw; 
      
    }

  }, std::runtime_error );

} // GenerateIndices, NegativeMaxSum 


TEST ( GenerateIndices, HugeMaxSum ) {

  EXPECT_THROW ({

    try {

      auto result = BasisFunctions::GenerateIndices<size_t> ( 5, 1000 );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "GenerateIndices: MaxSum not in [0,100] range", e.what() 
        
      );
      
      throw; 
      
    }


  }, std::runtime_error );

} // GenerateIndices, HugeMaxSum 


#endif // DEBUG_BASIS_FUNCTIONS 

