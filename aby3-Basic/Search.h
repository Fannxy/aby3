#pragma once
#include "Basics.h"
#include "../aby3-RTR/BuildingBlocks.h"

#include <aby3/Circuit/CircuitLibrary.h>
#include "../aby3-RTR/debug.h"


// the search functions correspond to the Blanton NDSS22 designs.
int mcompBS(std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int compBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int mtagBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int tagBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int subHBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);    