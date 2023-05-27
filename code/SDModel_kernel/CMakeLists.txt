### ===================================================================
### @file     CMakeLists.txt 
### @brief    Config and Build script for RSMD project
###
### >> ------------------------------------------------------------- >>
###
###             ######   #####  #     #  #####  ######  
###             #     # #     # ##   ## #     # #     # 
###             #     # #       # # # # #       #     # 
###             ######   #####  #  #  #  #####  #     # 
###             #   #         # #     #       # #     # 
###             #    #  #     # #     # #     # #     # 
###             #     #  #####  #     #  #####  ###### 
###
### >> ------------------------------------------------------------- >>
### 
### @copyright  TUM Copyright (C)
### @version    1.0
### @date       2023-05-02
### @author     Tanuharja, R.A. 
###
### Some description here
###
### ===================================================================

cmake_minimum_required ( VERSION 3.25 ) 
project ( ReducedSurrogateModel CXX ) 

set ( CMAKE_CXX_STANDARD 17 ) 


# ----- Options -------------------------------------------------------

option ( 
  TEST_RSMSD "Test RSMSD library" ON 
)


# ----- Add package "Eigen3" and "Boost" ------------------------------

find_package ( Eigen3 3.4 REQUIRED NO_MODULE ) 
find_package ( Boost 1.81 REQUIRED ) 


# ----- Add compile options based on compiler ID ----------------------

if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" ) 

  add_compile_options ( 
    -O2 -pedantic -Wall -Werror -fPIC -Wno-deprecated-declarations 
  )

elseif ( CMAKE_CXX_COMPILER_ID MATCHES "GNU" ) 

  add_compile_options ( 
    -O2 -pedantic -Wall -Werror -fPIC -Wno-deprecated-declarations 
  )

endif ()


# ----- Set up install path -------------------------------------------

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )

  set ( 
    CMAKE_INSTALL_PREFIX 
    "${CMAKE_BINARY_DIR}/install" 
    CACHE 
    PATH 
    "default install path" 
    FORCE 
  )

endif ()

set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}" )


# ----- Enable additional checks when in DEBUG build mode ------------- 

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )

  add_compile_definitions ( DEBUG_RSMSD )
  
endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )


# ----- Add modules ---------------------------------------------------

add_subdirectory ( modules/BasisFunctions )
add_subdirectory ( modules/Statistics )
add_subdirectory ( modules/SimpleModels )


# ----- Add shared library "ReducedSurrogateModel" --------------------

add_library ( ReducedSurrogateModel SHARED 

  declarations/Stochastic.hpp 
  declarations/submodules/DirectMC/MassSpringDamper_DMC.hpp 

  declarations/submodules/Surrogate/NIRPCE.hpp 

  implementations/DirectMC/MassSpringDamper_DMC.cpp 

  implementations/Surrogate/NIRPCE.cpp 

  utilities/LibrariesLoader_RS.hpp 
  
)

target_include_directories ( ReducedSurrogateModel PUBLIC 

  ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
  ${CMAKE_CURRENT_SOURCE_DIR}/declarations/submodules/DirectMC
  ${CMAKE_CURRENT_SOURCE_DIR}/declarations/submodules/Surrogate 
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities
  ${Boost_INCLUDE_DIR} 
  
)

target_link_libraries ( ReducedSurrogateModel PUBLIC 

  BasisFunctions 
  SimpleModels 
  Eigen3::Eigen 
  ${Boost_LIBRARIES} 

)


# ----- Configure Google Test -----------------------------------------

if ( TEST_RSMSD ) 

  find_package ( GTest REQUIRED ) 

  enable_testing () 

  add_executable ( RSMSD_testrunner 

    tests/DirectMC/MassSpringDamper_DMC_test.cpp 
    tests/Surrogate/NIRPCE_test.cpp 
    
  )

  target_include_directories ( RSMSD_testrunner PUBLIC

    ${GTEST_INCLUDE_DIRS} 
    
  )

  target_link_libraries ( RSMSD_testrunner PUBLIC 

    ReducedSurrogateModel 
    GTest::gtest_main 
    
  )

  gtest_discover_tests ( RSMSD_testrunner )
  
endif ( TEST_RSMSD )

