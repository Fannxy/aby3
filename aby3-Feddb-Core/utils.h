#ifndef _ABY3_FEDDB_UTILS_H_
#define _ABY3_FEDDB_UTILS_H_

#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-GORAM-Core/Basics.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3-GORAM-Core/Sort.h"

void prefixsum(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &result);

void prefixsum_with_initial_elements(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &result, aby3::si64Matrix &initial_elements);

void prefixsum_inv(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &result);

void permutation_inverse(aby3::i64Matrix& rsigma_plain, aby3::i64Matrix& rsigma_inv_plain);

// argsort not used
void fed_argsort(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &result, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void plain_argsort(aby3::i64Matrix& v, aby3::i64Matrix& result);

void permutate(int pIdx, aby3::si64Matrix &data, aby3::si64Matrix &res, std::vector<size_t> &permutation);

void permutate(int pIdx, std::vector<aby3::si64Matrix> &data, std::vector<aby3::si64Matrix> &res, std::vector<size_t> &permutation);

void permutate(int pIdx, aby3::si64Matrix &data, aby3::si64Matrix &res, aby3::i64Matrix  &permutation);

//void set_null_share(int pIdx, aby3::si64Matrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void set_const_share(int pIdx, aby3::i64 const_value, aby3::si64Matrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

#endif