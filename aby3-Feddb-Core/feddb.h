#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/RepeatThenReduce.h"
#include "../aby3-RTR/CipherIndex.h"
#include "../aby3-GORAM-Core/Basics.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3/sh3/Sh3Converter.h"
#include "../aby3/aby3-Feddb-Core/genperm.h"



#ifndef _ABY3_FEDDB_DAGOPNODE_H_
#define _ABY3_FEDDB_DAGOPNODE_H_


// 函数声明
// void concatRows_sbMatrix(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
//     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
//     aby3::Sh3Runtime &runtime);

// void concatRows_i64Matrix(int pIdx, aby3::i64Matrix &sharedA, aby3::i64Matrix &sharedB,
//     aby3::i64Matrix &res, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
//     aby3::Sh3Runtime &runtime);

void shuffle(int pIdx, aby3::si64Matrix& T, aby3::si64Matrix &Tres, 
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void shuffle(int pIdx, aby3::sbMatrix& T, aby3::sbMatrix &Tres, 
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void shuffle(int pIdx, std::vector<aby3::si64Matrix>& T, std::vector<aby3::si64Matrix>& Tres,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

template<typename MatrixType>
void persist_cipher(int pIdx, const std::string &table_name, MatrixType &T){
    std::string filename = "/root/GORAM-ABY3/aby3/aby3-Feddb-tmpfile/cipher/" + table_name + "_" + std::to_string(pIdx) + ".txt";
    
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing" << std::endl;
        return;
    }
    
    int rows = T.rows();
    //默认：n行1列
    
    outFile << "Matrix rows: " << rows  << std::endl;

    outFile << "Matrix shares[0]:" << std::endl;
    //T.mShares[0]和T.mShares[1]分别写入
    for (int i = 0; i < rows; i++) {
        outFile << T.mShares[0](i, 0) << " "<< std::endl;
    }
    outFile << "Matrix shares[1]:" << std::endl;
    for (int i = 0; i < rows; i++) {
        outFile << T.mShares[1](i, 0) << " "<< std::endl;
    }

    outFile.close();

    return;
}

template<typename MatrixType>
void read_cipher(int pIdx, const std::string &table_name, MatrixType &T){
    std::string filename = "/root/GORAM-ABY3/aby3/aby3-Feddb-tmpfile/cipher/" + table_name + "_" + std::to_string(pIdx) + ".txt";
    
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for reading" << std::endl;
        return;
    }

    std::string dummy;
    int rows;
    inFile >> dummy >> dummy >> rows;  
    T.resize(rows, 1);
    
    inFile >> dummy >> dummy; 
    for (int i = 0; i < rows; i++) {
        inFile >> T.mShares[0](i, 0);
    }
    
    inFile >> dummy >> dummy;  
    for (int i = 0; i < rows; i++) {
        inFile >> T.mShares[1](i, 0);
    }

    inFile.close();
    return;
}

void persist_plain(int pIdx, const std::string &table_name, aby3::i64Matrix &T);

void read_plain(int pIdx, const std::string &table_name, aby3::i64Matrix &T);

void oblivious_idx_select(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &idx, aby3::si64Matrix &result,              
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void index_agg(int pIdx, aby3::si64Matrix &equalFlag, aby3::i64Matrix &idx,std::vector<aby3::si64Matrix> &data,std::vector<aby3::si64Matrix> &finalRes,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

#endif