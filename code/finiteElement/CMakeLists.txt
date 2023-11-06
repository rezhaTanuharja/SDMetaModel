### ---------------------------------------------------------------------------
### @file		  CMakeLists.txt 
### @brief		Project configuration script for simple FE simulations
### @author		Rezha Adrian Tanuharja
### @date     2023-11-02
### ---------------------------------------------------------------------------

cmake_minimum_required ( VERSION 3.27 )
project ( finiteElement CXX )

set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_EXPORT_COMPILE_COMMANDS ON )

option ( HIGH_PRECISION "Use higher float and integer precision" OFF )
option ( ADDITIONAL_DEBUG "Additional check and debug messages" OFF )


# ----- Add packages ----------------------------------------------------------

# Eigen for linear algebra
# Spectra to compute eigenvectors and eigenvalues
# MPI for distributed memory parallelization
# OpenMP for shared memory parallelization

find_package ( Eigen3 REQUIRED NO_MODULE )
find_package ( Spectra REQUIRED NO_MODULE )
find_package ( MPI REQUIRED ) 
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

add_library ( minFEM SHARED

  inc/alias.hpp 
  inc/basisFunctions.hpp 
  inc/materials.hpp 
  inc/plates.hpp 
  inc/solvers.hpp 
  inc/substructuring.hpp
  inc/utilities.hpp 

  src/basisFunctions.cpp 
  src/materials.cpp 
  src/plates.cpp 
  src/solvers.cpp 
  src/substructuring.cpp 
  src/utilities.cpp 

)

target_include_directories ( minFEM PUBLIC 

  ${CMAKE_CURRENT_SOURCE_DIR}/inc 
  ${EIGEN3_INCLUDE_DIRS}
  ${SPECTRA_INCLUDE_DIRS}
  ${MPI_INCLUDE_PATH}

)

target_link_libraries ( minFEM

  Eigen3::Eigen
  Spectra::Spectra
  ${MPI_LIBRARIES}

)


# ----- Add drivers -----------------------------------------------------------

add_executable ( completeModel driver/completeModel.cpp ) 
add_executable ( substructuredModel driver/substructuredModel.cpp ) 

target_include_directories ( completeModel PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/inc 
)
target_link_libraries ( completeModel minFEM )

target_include_directories ( substructuredModel PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/inc 
)
target_link_libraries ( substructuredModel minFEM ${MPI_LIBRARIES} ) 
