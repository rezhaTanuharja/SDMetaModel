### ===================================================================
### @file     CMakeLists.txt 
### @brief    Config and Build script for Statistics module 
### @author   Tanuharja, R.A. 
### ===================================================================

cmake_minimum_required ( VERSION 3.25 ) 
project ( Statistics CXX ) 

set ( CMAKE_CXX_STANDARD 17 ) 


# ----- Options -------------------------------------------------------

option ( 
  TEST_STATISTICS "Test module Statistics" ON 
)


# ----- Add package "Eigen3" and "Boost" ------------------------------

find_package ( Eigen3 3.4 REQUIRED NO_MODULE ) 
find_package ( Boost 1.81 REQUIRED ) 


# ----- Add templated library -----------------------------------------

add_library ( Statistics INTERFACE 

  declarations/Statistics.hpp 

  declarations/submodules/Uniform.hpp 
  declarations/submodules/Normal.hpp 

  implementations/Uniform_imp.hpp 
  implementations/Normal_imp.hpp 

  utilities/LibrariesLoader_ST.hpp 
  
)

target_include_directories ( Statistics INTERFACE 

  ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
  ${CMAKE_CURRENT_SOURCE_DIR}/declarations/submodules 
  ${CMAKE_CURRENT_SOURCE_DIR}/implementations  
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities

  ${Boost_INCLUDE_DIR} 
  
)

target_link_libraries ( Statistics INTERFACE 

  Eigen3::Eigen 
  ${Boost_LIBRARIES} 

)


# ----- Enable additional checks when in DEBUG build mode ------------- 

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )

  add_compile_definitions ( DEBUG_STATISTICS )
  
endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )


# ----- Configure Google Test for Statistics --------------------------

if ( TEST_STATISTICS ) 

  find_package ( GTest REQUIRED ) 

  enable_testing () 

  add_executable ( Statistics_testrunner

    tests/uniform/RandomSampling_test.cpp 
    
    tests/normal/RandomSampling_test.cpp 
    
  )

  target_include_directories ( Statistics_testrunner PUBLIC

    ${GTEST_INCLUDE_DIRS} 
    
  )

  target_link_libraries ( Statistics_testrunner PUBLIC 

    Statistics 
    GTest::gtest_main 
    
  )

  gtest_discover_tests ( Statistics_testrunner )
  
endif ( TEST_STATISTICS )

