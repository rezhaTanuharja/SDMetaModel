/**
 * @file models.hpp
 *
 * @brief
 * Contains declarations of surrogate models.
 *
 * @author Rezha Adrian Tanuharja <rezhadr@outlook.com>
 * @date 2023-11-02
*/

#ifndef MODELS_DECLARATIONS
#define MODELS_DECLARATIONS

#ifndef UTILITIES_DECLARATIONS
#include "utilities.hpp"
#endif

#ifndef BASIS_FUNCTIONS_DECLARATIONS
#include "basisFunctions.hpp"
#endif

#include <omp.h>


namespace surrogate {


/**
 * @class RPCE
 *
 * @brief
 * Surrogate model in the form of a ratio between two PCEs.
*/
class RPCE {

  std::vector<INT> numIndices_;
  std::vector<INT> denIndices_;

  std::vector<FLOAT> numCoeffs_;
  std::vector<FLOAT> denCoeffs_;

  INT nDOFs_, dimension_;

public: // class RPCE


  /**
   * @brief
   * Constructor for RPCE surrogate model.
   *
   * @param nDOFs     number of DOFs in the model to approximate
   * @param dimension number of input parameters
  */
  RPCE ( const INT nDOFs, const INT dimension );


  /**
   * @brief
   * Set the numerator indices for RPCE model.
   *
   * @param indices the new indices of basis functions
  */
  void setNumIndices ( const std::vector<INT>& indices );


  /**
   * @brief
   * Set the denominator indices for RPCE model.
   *
   * @param indices the new indices of basis functions
  */
  void setDenIndices ( const std::vector<INT>& indices );


  /**
   * @brief
   * Compute the coefficients in the RPCE model.
   *
   * @param args      the training input parameters
   * @param responses the training data set
  */
  void train ( 
    const std::vector<FLOAT>& args, 
    const std::vector<FLOAT>& responses 
  );


  /**
   * @brief
   * Compute the coefficients in the RPCE model.
   *
   * @param args      the training input parameters
   * @param responses the training data set
   * @param epsilon   the tuning parameter
  */
  void trainSparse (
    const std::vector<FLOAT>& args,
    const std::vector<FLOAT>& responses, 
    const FLOAT epsilon
  );


  /**
   * @brief
   * Approximate responses using RPCE model.
   *
   * @param args  the input parameters
   * @return vector of approximate responses
  */
  std::vector<FLOAT> computeResponse ( const std::vector<FLOAT>& args ) const;


}; // class RPCE


} // namespace surrogate


#endif // MODELS_DECLARATIONS

