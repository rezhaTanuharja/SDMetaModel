/**
 * @file substructuring.hpp
 *
 * @brief
 * Contains declarations of functions for CB method.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef SUBSTRUCTURING_DECLARATIONS
#define SUBSTRUCTURING_DECLARATIONS

#ifndef ALIAS_DECLARATIONS
#include "alias.hpp"
#endif

#ifndef SOLVERS_DECLARATIONS
#include "solvers.hpp"
#endif

#include <Eigen/Dense>


namespace mFEM {


/**
 * @brief
 * Reduce component matrix and vector using the CB transformation.
 *
 * @param dynamicStiffness  the component's sparse dynamic stiffness matrix
 * @param loadVector        the component's dense load vector
 * @param internalDofIds    the component's internal DOFs' IDs
 * @param boundaryDofIds    the component's boundary DOFs' IDs
 * @param Psin              the component's dominant internal modes matrix
 * @param Psic              the component's constraint modes matrix
 *
 * @return { S, F, Dii, Dic, Fi }
*/
std::tuple<
  MatrixXF, VectorXF, MatrixXF, MatrixXF, VectorXF
> reduceComponent (
  const GlobalMatrix dynamicStiffness,
  const GlobalVector loadVector,
  const std::vector<INT>& internalDofIds,
  const std::vector<INT>& boundaryDofIds,
  const MatrixXF& Psin,
  const MatrixXF& Psic
);


} // namespace mFEM


#endif // SUBSTRUCTURING_DECLARATIONS

