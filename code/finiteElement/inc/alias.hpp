/**
 * @file alias.hpp
 *
 * @brief
 * Contains declarations convenient aliases.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef ALIAS_DECLARATIONS
#define ALIAS_DECLARATIONS

#include <array>
#include <functional>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace mFEM {


#ifdef HIGH_PRECISION
#define FLOAT long double
#define INT   long int
#else
#define FLOAT double
#define INT   int
#endif // HIGH_PRECISION


using ElementMatrix = Eigen::Matrix<FLOAT, 27, 27>;
using ElementVector = Eigen::Vector<FLOAT, 27>;

using ElementSystem = std::tuple<
  const std::vector<FLOAT>, const std::vector<FLOAT>, const std::vector<FLOAT>
>;

using GlobalMatrix = Eigen::SparseMatrix<FLOAT>;
using GlobalVector = Eigen::VectorXd;

using GlobalSystem = std::tuple<
  const GlobalMatrix, const GlobalMatrix, const GlobalVector
>;

using IntegrationPoints         = std::array<std::vector<FLOAT>, 2>;
using IntegrationPointsProvider = std::function<IntegrationPoints(INT)>;

using LocationMap  = std::vector<INT>;
using LocationMaps = std::vector<LocationMap>;

using MatrixXF = Eigen::Matrix<FLOAT, Eigen::Dynamic, Eigen::Dynamic>;

using SpatialFunction = std::function<FLOAT(FLOAT, FLOAT)>;

using VectorXF = Eigen::Vector<FLOAT, Eigen::Dynamic>;


} // namespace mFEM


#endif // ALIAS_DECLARATIONS

