#ifndef _ABY3_FEDDB_GENPERM_H_
#define _ABY3_FEDDB_GENPERM_H_

#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-GORAM-Core/Basics.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3-GORAM-Core/Sort.h"


// Function declarations
//void reshare(int pIdx, aby3::si64& x, int targetPartyIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime);
aby3::si64Matrix reshare_matrix(int pIdx, aby3::si64Matrix& x, int targetPartyIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime);
void getBitKey(int pIdx, aby3::sbMatrix &k_bool, aby3::si64Matrix &k_j_arith, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void genBitPerm(int pIdx, aby3::si64Matrix &k_j, aby3::si64Matrix &perm, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void applyPerm(int pIdx, aby3::si64Matrix &perm, aby3::si64Matrix &k_j, aby3::si64Matrix &k_j_prime, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void composePerm(int pIdx, aby3::si64Matrix &perm_a, aby3::si64Matrix &perm_b, aby3::si64Matrix &perm, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void genPerm(int pIdx, aby3::si64Matrix &k, aby3::si64Matrix &perm, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void genPerm_bool(int pIdx, aby3::sbMatrix &k, aby3::si64Matrix &perm, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void concat_k_v(int pIdx, aby3::sbMatrix& k, aby3::sbMatrix& v, aby3::sbMatrix& res);
//void genPerm_kv(int pIdx, aby3::si64Matrix &k,aby3::si64Matrix &v, aby3::si64Matrix &perm, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);


#endif