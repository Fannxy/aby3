#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/RepeatThenReduce.h"
#include "../aby3-RTR/CipherIndex.h"
#include "../aby3-GORAM-Core/Basics.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3/sh3/Sh3Converter.h"


#ifndef _ABY3_FEDDB_DAGOPNODE_H_
#define _ABY3_FEDDB_DAGOPNODE_H_


// 函数声明
// void concatRows_sbMatrix(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
//     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
//     aby3::Sh3Runtime &runtime);

// void concatRows_i64Matrix(int pIdx, aby3::i64Matrix &sharedA, aby3::i64Matrix &sharedB,
//     aby3::i64Matrix &res, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
//     aby3::Sh3Runtime &runtime);

// // void shuffle_sbMatrix(int pIdx, aby3::sbMatrix &T, aby3::sbMatrix &Tres, 
// //     aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

// void project_sbMatrix(int pIdx, aby3::sbMatrix &T, std::vector<int> &cols, aby3::sbMatrix &Tres, 
//     aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

// void project_i64Matrix(int pIdx, aby3::i64Matrix &T, std::vector<int> &cols, aby3::i64Matrix &Tres, 
//     aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

// template<typename MatrixType>
// void persist_cipher(int pIdx, std::string &table_name, MatrixType &T);

// void persist_plain(int pIdx, std::string &table_name, aby3::i64Matrix &T);

void oblivious_idx_select(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &idx, aby3::si64Matrix &result,              
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

#endif