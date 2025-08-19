#ifndef _ABY3_FEDDB_UTILS_H_
#define _ABY3_FEDDB_UTILS_H_

#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-GORAM-Core/Basics.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3-GORAM-Core/Sort.h"

void prefixsum(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &result);

void prefixsum_inv(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &result);

void permutation_inverse(aby3::i64Matrix& rsigma_plain, aby3::i64Matrix& rsigma_inv_plain);

void fed_argsort(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &result, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

#endif