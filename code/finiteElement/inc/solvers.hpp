/**
 * @file solvers.hpp
 *
 * @brief
 * Contains declarations of functions to solve system of linear eqs.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef SOLVERS_DECLARATIONS
#define SOLVERS_DECLARATIONS

#ifndef ALIAS_DECLARATIONS
#include "alias.hpp"
#endif

#include <Eigen/IterativeLinearSolvers>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymGEigsShiftSolver.h>


namespace mFEM {


/**
 * @brief
 * Compute undamped dynamic stiffness matrix.
 *
 * @param massMatrix        sparse mass matrix
 * @param stiffnessMatrix   sparse stiffness matrix
 * @param omega             angular speed
 *
 * @return sparse undamped dynamic stiffness matrix
*/
GlobalMatrix computeDynamicStiffness (
  const GlobalMatrix& massMatrix,
  const GlobalMatrix& stiffnessMatrix, 
  const FLOAT omega
);


/**
 * @brief
 * Apply homogenous Dirichlet boundary condition to system of linear eqs.
 *
 * @param dynamicStiffness  sparse dynamic stiffness matrix
 * @param loadVector        dense load vector
 * @param boundaryDofIds    ids of boundary DOFs
*/
void applyHomogenousDirichletBC (
  GlobalMatrix& dynamicStiffness,
  GlobalVector& loadVector,
  const std::vector<INT>& boundaryDofIds
);


/**
 * @brief
 * Solve sparse dynamic system of eqs.
 *
 * @details
 * Solve a system of linear equations in which the matrix is sparse
 * and the vector is dense. Uses preconditioned conjugate gradient
 * solver with incomplete Cholesky preconditioner.
 *
 * @param dynamicStiffness  sparse dynamic stiffness matrix
 * @param loadVector        dense load vector
 * @param maxIteration      maximum number of iterations
 * @param tolerance         relative tolerance
 *
 * @return pair: ( solution, convergence status )
*/
std::pair<GlobalVector,bool> solveDynamicSystem (
  const GlobalMatrix dynamicStiffness,
  const GlobalVector loadVector,
  INT maxIteration,
  FLOAT tolerance
);


/**
 * @brief
 * Extract subsets of rows and columns of a sparse matrix.
 *
 * @param matrix  sparse matrix to extract from
 * @param rows    vector of row indices to extract
 * @param cols    vector of col indices to extract
 *
 * @return a smaller sparse matrix
*/
GlobalMatrix extract (
  const GlobalMatrix& matrix,
  const std::vector<INT>& rows,
  const std::vector<INT>& cols
);


/**
 * @brief
 * Extract subsets of rows of a vector.
 *
 * @param vector  a vector to extract from
 * @param rows    vector of row indices to extract
 *
 * @return a smaller vector
*/
GlobalVector extract (
  const GlobalVector& vector,
  const std::vector<INT>& rows
);


/**
 * @brief
 * Compute eigenvectors from a generalized eigenproblem K - lambda x M.
 *
 * @details
 * Compute eigenvectors from a generalized eigenproblem where the matrices
 * are symmetric. Uses shift and invert technique, with shift equals to
 * squared omega. The function returns n eigenvectors that has eigenvalues
 * closest to squared omega.
 *
 * @param M       sparse mass matrix
 * @param K       sparse stiffness matrix
 * @param omega   angular frequency
 * @param n       number of eigenvectors
 *
 * @return eigenvectors stored in matrix's columns
*/
MatrixXF computeEigenvectors (
  const GlobalMatrix& M,
  const GlobalMatrix& K, 
  const FLOAT omega,
  const INT n
);


} // namespace mFEM


#endif // SOLVERS_DECLARATIONS

