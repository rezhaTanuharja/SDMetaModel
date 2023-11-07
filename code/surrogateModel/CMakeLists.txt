### ---------------------------------------------------------------------------
### @file     CMakeLists.txt 
### @brief    Configuration script for RPCE model project
### @author   Rezha Adrian Tanuharja
### @date     2023-11-02
### ---------------------------------------------------------------------------

cmake_minimum_required ( VERSION 3.27 )
project ( surrogateModel CXX )

set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_EXPORT_COMPILE_COMMANDS ON )

option ( HIGH_PRECISION "Use higher float and integer precision" OFF )
option ( ADDITIONAL_DEBUG "Additional check and debug messages" OFF )


# ----- Add packages ----------------------------------------------------------

# Eigen for linear algebra
# OpenMP for shared memory parallelization

find_package ( Eigen3 REQUIRED NO_MODULE )
find_package ( OpenMP )

if ( OPENMP_FOUND ) 
  SET ( CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} )
	SET ( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} )
endif ()


# ----- Add suitable compile options for common compilers ---------------------

# Maximum optimization for speed
# Ignore warning about deprecated declarations from Eigen & Spectra
# Enable OpenMP

if ( CMAKE_CXX_COMPILER_ID MATCHES "GNU" ) 

  # GNU Compiler (g++)
  add_compile_options ( -O3 -Wno-deprecated-declarations -fopenmp ) 

elseif ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" ) 

  # Clang Compiler
  add_compile_options ( -O3 -Wno-deprecated -fopenmp ) 

elseif ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" ) 

  # Microsoft Visual C++ (MSVC) Compiler
  add_compile_options ( /O2 /EHsc /openmp )

elseif ( CMAKE_CXX_COMPILER_ID MATCHES "INTEL" ) 

  # Intel C++ Compiler
  add_compile_options ( -O2 -no-ftz -openmp )

endif ()


# ----- Precision settings ----------------------------------------------------

if ( HIGH_PRECISION ) 

  add_compile_definitions ( HIGH_PRECISION )

endif ()


# ----- Additional debug settings ---------------------------------------------

if ( ADDITIONAL_DEBUG ) 

  add_compile_definitions ( ADDITIONAL_DEBUG )

endif ()


# ----- Create shared library minFEM ------------------------------------------

# A minimalist FE library

add_library ( surrogateModel SHARED

  inc/basisFunctions.hpp 
  inc/models.hpp 
  inc/regression.hpp 
  inc/utilities.hpp 

  src/basisFunctions.cpp 
  src/models.cpp 
  src/regression.cpp 
  src/utilities.cpp 

)

target_include_directories ( surrogateModel PUBLIC 

  ${CMAKE_CURRENT_SOURCE_DIR}/inc 
  ${EIGEN3_INCLUDE_DIRS}

)

target_link_libraries ( surrogateModel

  Eigen3::Eigen

)


# ----- Add drivers -----------------------------------------------------------

add_executable ( regularRPCE driver/regularRPCE.cpp ) 
add_executable ( sparseRPCE driver/sparseRPCE.cpp ) 

target_include_directories ( regularRPCE PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/inc 
)
target_link_libraries ( regularRPCE surrogateModel )

target_include_directories ( sparseRPCE PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/inc 
)
target_link_libraries ( sparseRPCE surrogateModel )
