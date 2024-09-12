#include "Basics.h"

#include <aby3/Circuit/CircuitLibrary.h>
#include "../aby3-RTR/debug.h"

/**
 * Shuffles the elements of the input vector T using the specified parameters.
 * 
 * @param T The input vector to be shuffled.
 * @param pIdx The index of the permutation to be used for shuffling.
 * @param Tres The output vector to store the shuffled elements.
 * @param enc The Sh3Encryptor object for encryption operations.
 * @param eval The Sh3Evaluator object for evaluation operations.
 * @param runtime The Sh3Runtime object for runtime operations.
 * @return The number of elements successfully shuffled.
 */
int efficient_shuffle(std::vector<aby3::sbMatrix> &T, int pIdx, std::vector<aby3::sbMatrix> &Tres, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int efficient_shuffle(aby3::sbMatrix &T, int pIdx, aby3::sbMatrix &Tres, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int efficient_shuffle_with_random_permutation(std::vector<aby3::sbMatrix> &T, int pIdx, std::vector<aby3::sbMatrix> &Tres, std::vector<aby3::si64> &Pi, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int switch_given_bits(aby3::sbMatrix &left, aby3::sbMatrix &right, int pIdx, aby3::sbMatrix &switch_bit, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int random_switch(aby3::sbMatrix &left, aby3::sbMatrix &right, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int random_switch(std::vector<aby3::sbMatrix>& left, std::vector<aby3::sbMatrix>& right, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int permutation_network(std::vector<aby3::sbMatrix> &T, int pIdx, std::vector<aby3::sbMatrix> &Tres, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);