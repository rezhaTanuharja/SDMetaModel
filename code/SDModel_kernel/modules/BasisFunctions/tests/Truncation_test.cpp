/**
  * @file Truncation_test.cpp 
  *
  * @brief
  * test functions to truncate polynomial indices, 
  * including erros thrown in DEBUG mode
  *
  * @author
  * Tanuharja R.A. @n 
  * contact: rezha.tanuharja@tum.de 
  */

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


#ifdef DEBUG_BASIS_FUNCTIONS 


TEST ( RemoveLargeIndices, NegativeSetSize ) {

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

  EXPECT_THROW ({

    try {

      BasisFunctions::RemoveLargeIndices<size_t> ( 
        Indices, -1, 2
      );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "RemoveLargeIndices: SetSize not in [1,100] range", e.what() 
        
      );
      
      throw; 
      
    }


  }, std::runtime_error );

} // RemoveLargeIndices, NegativeSetSize 


TEST ( RemoveLargeIndices, ZeroSetSize ) {

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

  EXPECT_THROW ({

    try {

      BasisFunctions::RemoveLargeIndices<size_t> ( 
        Indices, 0, 2
      );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "RemoveLargeIndices: SetSize not in [1,100] range", e.what() 
        
      );
      
      throw; 
      
    }


  }, std::runtime_error );

} // RemoveLargeIndices, ZeroSetSize  


TEST ( RemoveLargeIndices, HugeSetSize ) {

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

  EXPECT_THROW ({

    try {

      BasisFunctions::RemoveLargeIndices<size_t> ( 
        Indices, 1000, 2
      );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "RemoveLargeIndices: SetSize not in [1,100] range", e.what() 
        
      );
      
      throw; 
      
    }


  }, std::runtime_error );

} // RemoveLargeIndices, HugeSetSize  


TEST ( RemoveLargeIndices, NegativeIMax ) {

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

  EXPECT_THROW ({

    try {

      BasisFunctions::RemoveLargeIndices<size_t> ( 
        Indices, 3, -1
      );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "RemoveLargeIndices: iMax not in [0,100] range", e.what() 
        
      );
      
      throw; 
      
    }


  }, std::runtime_error );

} // RemoveLargeIndices, NegativeIMax 


TEST ( RemoveLargeIndices, HugeIMax ) {

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

  EXPECT_THROW ({

    try {

      BasisFunctions::RemoveLargeIndices<size_t> ( 
        Indices, 3, 1000
      );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "RemoveLargeIndices: iMax not in [0,100] range", e.what() 
        
      );
      
      throw; 
      
    }


  }, std::runtime_error );

} // RemoveLargeIndices, HugeIMax 


TEST ( RemoveLargeIndices, MismatchSetSize ) {

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

  EXPECT_THROW ({

    try {

      BasisFunctions::RemoveLargeIndices<size_t> ( 
        Indices, 4, 1
      );

    } catch ( const std::exception& e ) {

      EXPECT_STREQ (
        
        "RemoveLargeIndices: Indices vector size not multiple of SetSize",
        e.what()
        
      );
      
      throw; 
      
    }


  }, std::runtime_error );

} // RemoveLargeIndices, MismatchSetSize  


#endif 

