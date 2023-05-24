### ===================================================================
### @file CMakeLists.txt 
### @brief Config and Build script for BasisFunctions module 
### @authors    Tanuharja, R.A. 
### ===================================================================

cmake_minimum_required ( VERSION 3.25 ) 
project ( BasisFunctions CXX ) 

set ( CMAKE_CXX_STANDARD 17 ) 


# ----- Options -------------------------------------------------------

option ( 
  TEST_BASIS_FUNCTIONS "Test module BasisFunctions" ON 
)


# ----- Add templated library -----------------------------------------

add_library ( BasisFunctions INTERFACE 

  declarations/BasisFunctions.hpp 

  implementations/GenerateIndices_imp.hpp 
  implementations/ComputeBasisFunctions_imp.hpp 

  utilities/LibrariesLoader_BF.hpp 
  
)

target_include_directories ( BasisFunctions INTERFACE 
  
  ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
  ${CMAKE_CURRENT_SOURCE_DIR}/implementations  
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities
  
)

target_link_libraries ( BasisFunctions INTERFACE 

)


# ----- Enable additional checks when in DEBUG build mode ------------- 

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )

  add_compile_definitions ( DEBUG_BASIS_FUNCTIONS )

endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )


# ----- Configure Google Test for BasisFunctions ----------------------

if ( TEST_BASIS_FUNCTIONS ) 

  find_package ( GTest REQUIRED ) 

  enable_testing () 

  add_executable ( BasisFunctions_testrunner

    tests/GenerateIndices_test.cpp 
    tests/Truncation_test.cpp 
    tests/HermitePolynomial_test.cpp 
    tests/HermiteBasis_test.cpp 
    
  )

  target_include_directories ( BasisFunctions_testrunner PUBLIC

    ${GTEST_INCLUDE_DIRS} 

  )

  target_link_libraries ( BasisFunctions_testrunner PUBLIC 

    BasisFunctions 
    GTest::gtest_main 

  )

  gtest_discover_tests ( BasisFunctions_testrunner )
  
endif ( TEST_BASIS_FUNCTIONS )
