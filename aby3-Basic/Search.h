#pragma once
#include <aby3/Circuit/CircuitLibrary.h>
#include "Basics.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/PtATasks.h"


// the search functions correspond to the Blanton NDSS22 designs.
int mcompBS(std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int compBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int mtagBS(std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int tagBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int subHBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, int alpha, int threshold = 1<<8);

int ptaBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, ABY3MPITask<aby3::sb64, aby3::sb64, aby3::si64, aby3::si64, PtABS>* ptaTask);

int ptaMBS(std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, ABY3MPIPairOnlyTask<aby3::sb64, aby3::sb64, aby3::sb64, aby3::sb64, PtAMBS>* ptaTask);

int subHBS_with_PtA(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, int alpha, int threshold, ABY3MPITask<aby3::sb64, aby3::sb64, aby3::si64, aby3::si64, PtABS>* ptaBSTask, ABY3MPIPairOnlyTask<aby3::sb64, aby3::sb64, aby3::sb64, aby3::sb64, PtAMBS>* ptaMBSTask);