/**
 * @file solvers.cpp
 *
 * @brief
 * Contains implementations of functions to solve system of linear eqs.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef SOLVERS_DECLARATIONS
#include "solvers.hpp"
#endif


namespace mFEM {


GlobalMatrix computeDynamicStiffness (
  const GlobalMatrix& massMatrix,
  const GlobalMatrix& stiffnessMatrix, 
  const FLOAT omega
) {
  return stiffnessMatrix - omega * omega * massMatrix;
} 


void applyHomogenousDirichletBC (
  GlobalMatrix& dynamicStiffness,
  GlobalVector& loadVector,
  const std::vector<INT>& boundaryDofIds
) {

  // create mask matrix: identity matrix with diagonals associated with 
  // boundary DOFs set to zero

  Eigen::Vector<INT, Eigen::Dynamic> mask = 
    Eigen::Vector<INT, Eigen::Dynamic>::Ones ( 
      dynamicStiffness.cols() 
    );

  for ( INT i = 0; i < boundaryDofIds.size(); i++ ) {
    mask ( boundaryDofIds[i] ) = 0;
  }

  GlobalMatrix maskMatrix (
    dynamicStiffness.rows(), 
    dynamicStiffness.cols()
  );

  for ( INT i = 0; i < mask.size(); i++ ) {
    if ( mask(i) == 1 ) {
      maskMatrix.insert ( i, i ) = 1.0;
    }
  }


  // set rows and columns of boundary DOFs to zero

  dynamicStiffness = maskMatrix * dynamicStiffness * maskMatrix;
  loadVector       = maskMatrix * loadVector;


  // set diagonals associated with boundary DOFs in the matrix to ones

  for ( INT i = 0; i < boundaryDofIds.size(); i++ ) {
    INT dof = boundaryDofIds[i];
    dynamicStiffness.coeffRef ( dof, dof ) = 1.0;
  }

} // applyHomogenousDirichletBC


std::pair<GlobalVector,bool> solveDynamicSystem (
  const GlobalMatrix dynamicStiffness,
  const GlobalVector loadVector,
  INT maxIteration,
  FLOAT tolerance
) {

  Eigen::ConjugateGradient <
    GlobalMatrix,
    Eigen::Lower | Eigen::Upper,
    Eigen::IncompleteCholesky<FLOAT>
  > solver;

  solver.setMaxIterations ( maxIteration );
  solver.setTolerance ( tolerance );

  solver.preconditioner().compute ( dynamicStiffness );
  solver.compute ( dynamicStiffness );

  GlobalVector U = solver.solve ( loadVector );

  return { U, solver.info() == Eigen::Success };

} // solveDynamicSystem


GlobalMatrix extract (
  const GlobalMatrix& matrix,
  const std::vector<INT>& rows,
  const std::vector<INT>& cols
) {

  GlobalMatrix extractedMatrix ( rows.size(), cols.size() );
  extractedMatrix.reserve ( rows.size() * cols.size() );

  for ( INT j = 0; j < cols.size(); j++ ) {

    INT col = cols[j];

    extractedMatrix.startVec ( j );

    for ( GlobalMatrix::InnerIterator it ( matrix, col ); it; ++it ) {

      for ( INT i = 0; i < rows.size(); i++ ) {

        INT row = rows[i];

        if ( it.row() == row ) {
          extractedMatrix.insertBack ( i, j ) = it.value();
          break;
        }

      }
    }
  }

  extractedMatrix.finalize();

  return extractedMatrix;

} // extract


GlobalVector extract (
  const GlobalVector& vector,
  const std::vector<INT>& rows
) {

  GlobalVector extractedVector = GlobalVector::Zero ( rows.size() );

  for ( INT i = 0; i < rows.size(); i++ ) {
    extractedVector(i) = vector(rows[i]);
  }

  return extractedVector;

} // extract


MatrixXF computeEigenvectors (
  const GlobalMatrix& M,
  const GlobalMatrix& K, 
  const FLOAT omega,
  const INT n
) {

  using OpType  = Spectra::SymShiftInvert<FLOAT,Eigen::Sparse,Eigen::Sparse>;
  using BopType = Spectra::SparseSymMatProd<FLOAT>;

  OpType   op ( K, M );
  BopType bop (    M );

  INT acc = ( n + M.rows() ) / 2;

  Spectra::SymGEigsShiftSolver<
    OpType, BopType, Spectra::GEigsMode::ShiftInvert
  > solver ( op, bop, n, acc, omega * omega );

  solver.init ();
  solver.compute ();

  return solver.eigenvectors();

} // computeEigenvectors


} // namespace mFEM

