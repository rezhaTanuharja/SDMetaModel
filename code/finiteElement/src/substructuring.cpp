/**
 * @file substructuring.cpp
 *
 * @brief
 * Contains implementations of functions for CB method.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef SUBSTRUCTURING_DECLARATIONS
#include "substructuring.hpp"
#endif


namespace mFEM {


std::tuple<
  MatrixXF, VectorXF, MatrixXF, MatrixXF, VectorXF
> reduceComponent (
  const GlobalMatrix dynamicStiffness,
  const GlobalVector loadVector,
  const std::vector<INT>& internalDofIds,
  const std::vector<INT>& boundaryDofIds,
  const MatrixXF& Psin,
  const MatrixXF& Psic
) {


  // Partition dynamic stiffness matrix and load vector

  GlobalMatrix Dii = extract ( 
    dynamicStiffness, internalDofIds, internalDofIds 
  ); 

  GlobalMatrix Dic = extract ( 
    dynamicStiffness, internalDofIds, boundaryDofIds 
  ); 

  GlobalMatrix Dcc = extract ( 
    dynamicStiffness, boundaryDofIds, boundaryDofIds 
  ); 

  GlobalVector Fi = extract ( loadVector, internalDofIds );
  GlobalVector Fc = extract ( loadVector, boundaryDofIds );


  // Rayleigh-Ritz transformation into a modal space

  MatrixXF Dnn = Psin.transpose() * Dii * Psin;

  MatrixXF Dnb = Psin.transpose() * Dii * Psic + Psin.transpose() * Dic;
  MatrixXF Dbb = Psic.transpose() * Dii * Psic + Psic.transpose() * Dic +
    Dic.transpose() * Psic + Dcc;

  VectorXF Fn = Psin.transpose() * Fi;
  VectorXF Fb = Psic.transpose() * Fi + Fc;


  // Condense dynamic stiffness matrix and load vector to boundary DOFs

  MatrixXF S = Dbb - Dnb.transpose() * Dnn.ldlt().solve(Dnb);
  VectorXF F = Fb  - Dnb.transpose() * Dnn.ldlt().solve(Fn);

  return { S, F, Dnn, Dnb, Fn };

} // reduceComponent


} // namespace mFEM

