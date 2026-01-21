#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/RepeatThenReduce.h"
#include "../aby3-RTR/CipherIndex.h"
#include "../aby3-GORAM-Core/Basics.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3/sh3/Sh3Converter.h"
#include "../aby3/aby3-Feddb-Core/genperm.h"
#include <cmath>



#ifndef _ABY3_FEDDB_DAGOPNODE_H_
#define _ABY3_FEDDB_DAGOPNODE_H_


void both2cipher(int pIdx, std::vector<int> &plain_cols_idx, std::vector<int> &cipher_cols_idx, 
    std::vector<aby3::i64Matrix> &input_plain_cols, std::vector<aby3::si64Matrix> &input_cipher_cols,
    std::vector<aby3::si64Matrix> &input_cols,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void shuffle(int pIdx, aby3::si64Matrix& T, aby3::si64Matrix &Tres, 
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void shuffle(int pIdx, aby3::sbMatrix& T, aby3::sbMatrix &Tres, 
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void shuffle(int pIdx, std::vector<aby3::si64Matrix>& T, std::vector<aby3::si64Matrix>& Tres,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

//si64 & sb
template<typename MatrixType>
void persist_cipher(int pIdx, const std::string &table_name, std::vector<MatrixType> &T){
    std::string filename = "./aby3-Feddb-tmpfile/cipher/" + table_name + "_" + std::to_string(pIdx) + ".txt";
    
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing" << std::endl;
        return;
    }
    
    
    int cols = T.size();
    outFile << "Matrix cols: " << cols << std::endl;
    int rows = T[0].rows();
    outFile << "Matrix rows: " << rows << std::endl;

    for(size_t i = 0; i < T.size(); i++){
        outFile << "Matrix_" << i << " shares[0]:" << std::endl;
    //T.mShares[0]和T.mShares[1]分别写入
        for (int j = 0; j < rows; j++) {
            outFile << T[i].mShares[0](j, 0) << std::endl;
        }
        outFile << "Matrix_" << i << " shares[1]:" << std::endl;
        for (int j = 0; j < rows; j++) {
            outFile << T[i].mShares[1](j, 0) << std::endl;
        }
    }

    outFile.close();

    return;
}

template<typename MatrixType>
void read_cipher(int pIdx, const std::string &table_name, std::vector<MatrixType> &T){
    std::string filename = "./aby3-Feddb-tmpfile/cipher/" + table_name + "_" + std::to_string(pIdx) + ".txt";
    
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for reading" << std::endl;
        return;
    }

    std::string dummy;
    int cols;
    inFile >> dummy >> dummy >> cols;
    T.resize(cols);
    
    int rows;
    inFile >> dummy >> dummy >> rows;  
    for(size_t i = 0; i < T.size(); i++){
        T[i].resize(rows, 1);
        inFile >> dummy >> dummy;
        for (int j = 0; j < rows; j++) {
            inFile >> T[i].mShares[0](j, 0);
        }

        inFile >> dummy >> dummy;
        for (int j = 0; j < rows; j++) {
            inFile >> T[i].mShares[1](j, 0);
        }
    }

    inFile.close();
    return;
}

void persist_plain(int pIdx, const std::string &table_name, std::vector<aby3::i64Matrix> &T);

void read_plain(int pIdx, const std::string &table_name, std::vector<aby3::i64Matrix> &T);

void oblivious_idx_select(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &idx, aby3::si64Matrix &result,              
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void index_agg(int pIdx, aby3::si64Matrix &equalFlag, aby3::i64Matrix &idx,std::vector<aby3::si64Matrix> &data_key,aby3::si64Matrix &data_val,std::vector<aby3::si64Matrix> &finalRes,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void group_by_common(int pIdx, aby3::si64Matrix &key, aby3::si64Matrix &val, 
    aby3::si64Matrix &key_g, aby3::si64Matrix &val_g, aby3::si64Matrix &e, aby3::si64Matrix &perm_GN, aby3::si64Matrix &key_GN, aby3::si64Matrix &key_out,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void group_by_common_without_val(int pIdx, aby3::si64Matrix &key, 
    aby3::si64Matrix &key_g, aby3::si64Matrix &e, aby3::si64Matrix &perm_GN, aby3::si64Matrix &key_GN, aby3::si64Matrix &key_out,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void group_count(int pIdx, std::vector<aby3::si64Matrix> &key, 
    std::vector<aby3::si64Matrix> &key_out, aby3::si64Matrix &c, 
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void group_sum(int pIdx, std::vector<aby3::si64Matrix> &key, aby3::si64Matrix &val, 
    std::vector<aby3::si64Matrix> &key_out, aby3::si64Matrix &sum,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void group_max(int pIdx, std::vector<aby3::si64Matrix> &key, aby3::si64Matrix &val, 
    std::vector<aby3::si64Matrix> &key_out, aby3::si64Matrix &max,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void group_min(int pIdx, std::vector<aby3::si64Matrix> &key, aby3::si64Matrix &val, 
    std::vector<aby3::si64Matrix> &key_out, aby3::si64Matrix &min,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void augment_table(int pIdx, std::vector<aby3::si64Matrix> &T_1_key, std::vector<aby3::si64Matrix> &T_1_other, 
    std::vector<aby3::si64Matrix> &T_2_key, std::vector<aby3::si64Matrix> &T_2_other,
    std::vector<aby3::sbMatrix> &T_1_auged, std::vector<aby3::sbMatrix> &T_2_auged,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void oblivious_expand(int pIdx, std::vector<aby3::sbMatrix> &T, std::vector<aby3::sbMatrix> &A, aby3::i64 tid,
    aby3::i64Matrix &s_plain,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void oblivious_distribute(int pIdx, std::vector<aby3::sbMatrix> &T_prime, aby3::sbMatrix &flag, aby3::sbMatrix &fx, aby3::i64Matrix &s_plain,
    aby3::sbMatrix &A_vector, aby3::sbMatrix &flag_sorted_auged,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void align_table(int pIdx, std::vector<aby3::sbMatrix> &T, std::vector<aby3::sbMatrix> &T_aligned,
        aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void join(int pIdx, std::vector<aby3::si64Matrix> &T_1_key, std::vector<aby3::si64Matrix> &T_1_other,
    std::vector<aby3::si64Matrix> &T_2_key, std::vector<aby3::si64Matrix> &T_2_other,
    std::vector<aby3::si64Matrix> &T_joined,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);


void filter(int pIdx, std::vector<aby3::si64Matrix> &T, int filColIdx, int value, bool is_scalar, std::string op_str,
    std::vector<aby3::si64Matrix> &T_filtered,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
 
#endif