### ===================================================================
### @file     CMakeLists.txt 
### @brief    Config and Build script for SimpleModels module 
### @author   Tanuharja, R.A. 
### ===================================================================

cmake_minimum_required ( VERSION 3.25 ) 
project ( SimpleModels CXX ) 

set ( CMAKE_CXX_STANDARD 17 ) 


# ----- Options -------------------------------------------------------

option ( 
  TEST_SIMPLE_MODELS "Test module SimpleModels" ON 
)


# ----- Add package "Eigen3" and "Boost" ------------------------------

find_package ( Eigen3 3.4 REQUIRED NO_MODULE ) 


# ----- Add templated library -----------------------------------------

add_library ( SimpleModels INTERFACE 

  declarations/SimpleModels.hpp 

  declarations/submodules/MassSpringDamper.hpp 

  implementations/MassSpringDamper_imp.hpp 

  utilities/LibrariesLoader_SM.hpp 

)

target_include_directories ( SimpleModels INTERFACE 

  ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
  ${CMAKE_CURRENT_SOURCE_DIR}/declarations/submodules 
  ${CMAKE_CURRENT_SOURCE_DIR}/implementations  
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities
  
)

target_link_libraries ( SimpleModels INTERFACE 

  Eigen3::Eigen 
  
)


# ----- Enable additional checks when in DEBUG build mode ------------- 

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )

  add_compile_definitions ( DEBUG_SIMPLE_MODELS )
  
endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )


# ----- Configure Google Test for Statistics --------------------------

if ( TEST_SIMPLE_MODELS ) 

  find_package ( GTest REQUIRED ) 

  enable_testing () 

  add_executable ( SimpleModels_testrunner 

    tests/massSpringDamper/MassMatrix_test.cpp 
    tests/massSpringDamper/StiffnessMatrix_test.cpp 
    tests/massSpringDamper/DynamicStiffnessMatrix_test.cpp 
    
  )

  target_include_directories ( SimpleModels_testrunner PUBLIC

    ${GTEST_INCLUDE_DIRS} 
    
  )

  target_link_libraries ( SimpleModels_testrunner PUBLIC 

    SimpleModels 
    GTest::gtest_main 
    
  )

  gtest_discover_tests ( SimpleModels_testrunner )
  
endif ( TEST_SIMPLE_MODELS )

