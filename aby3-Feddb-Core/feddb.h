#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/RepeatThenReduce.h"
#include "../aby3-RTR/CipherIndex.h"
#include "../aby3-GORAM-Core/Basics.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3/sh3/Sh3Converter.h"
#include "../aby3/aby3-Feddb-Core/genperm.h"
#include <cmath>
#include <cstdio>



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
    std::string filename = "./aby3-Feddb-tmpfile/cipher/" + table_name + "_" + std::to_string(pIdx) + ".bin";

    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error: Unable to open file " << filename << " for writing" << std::endl;
        return;
    }

    int cols = T.size();
    int rows = (cols > 0) ? static_cast<int>(T[0].rows()) : 0;
    fwrite(&cols, sizeof(int), 1, fp);
    fwrite(&rows, sizeof(int), 1, fp);

    for(size_t i = 0; i < T.size(); i++){
        fwrite(T[i].mShares[0].data(), sizeof(T[i].mShares[0](0, 0)), rows, fp);
        fwrite(T[i].mShares[1].data(), sizeof(T[i].mShares[1](0, 0)), rows, fp);
    }

    fclose(fp);
    return;
}

template<typename MatrixType>
void read_cipher(int pIdx, const std::string &table_name, std::vector<MatrixType> &T){
    std::string filename = "./aby3-Feddb-tmpfile/cipher/" + table_name + "_" + std::to_string(pIdx) + ".bin";

    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        std::cerr << "Error: Unable to open file " << filename << " for reading" << std::endl;
        return;
    }

    int cols, rows;
    fread(&cols, sizeof(int), 1, fp);
    fread(&rows, sizeof(int), 1, fp);
    if (cols <= 0 || rows <= 0) {
        T.clear();
        fclose(fp);
        return;
    }
    T.resize(cols);

    for(size_t i = 0; i < T.size(); i++){
        T[i].resize(rows, 1);
        fread(T[i].mShares[0].data(), sizeof(T[i].mShares[0](0, 0)), rows, fp);
        fread(T[i].mShares[1].data(), sizeof(T[i].mShares[1](0, 0)), rows, fp);
    }

    fclose(fp);
    return;
}

void persist_plain(int pIdx, const std::string &table_name, std::vector<aby3::i64Matrix> &T);

void read_plain(int pIdx, const std::string &table_name, std::vector<aby3::i64Matrix> &T);

void oblivious_idx_select(int pIdx, aby3::si64Matrix &v, aby3::si64Matrix &idx, aby3::si64Matrix &result,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void oblivious_idx_select(int pIdx, std::vector<aby3::si64Matrix> &v_vec, aby3::si64Matrix &idx, std::vector<aby3::si64Matrix> &results,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void index_agg(int pIdx, aby3::si64Matrix &equalFlag, aby3::i64Matrix &idx,std::vector<aby3::si64Matrix> &data_key,aby3::si64Matrix &data_val,std::vector<aby3::si64Matrix> &finalRes,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

// Segmented MAX/MIN aggregation over rows sorted by idx.
// equalFlag[i] = 1 means row i belongs to the same group as row i-1.
// is_max=true computes MAX, false computes MIN.
void index_agg_maxmin(int pIdx, aby3::si64Matrix &equalFlag, aby3::i64Matrix &idx,
    std::vector<aby3::si64Matrix> &data_key, aby3::si64Matrix &data_val, bool is_max,
    std::vector<aby3::si64Matrix> &finalRes,
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
    std::vector<aby3::si64Matrix> &T_1_auged, std::vector<aby3::si64Matrix> &T_2_auged,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void oblivious_expand(int pIdx, std::vector<aby3::si64Matrix> &T, std::vector<aby3::sbMatrix> &A, aby3::i64 tid,
    aby3::i64Matrix &s_plain,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void oblivious_distribute(int pIdx, std::vector<aby3::si64Matrix> &T_prime, aby3::sbMatrix &flag, aby3::si64Matrix &fx, aby3::i64Matrix &s_plain,
    aby3::sbMatrix &A_vector, aby3::sbMatrix &flag_sorted_auged,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void align_table(int pIdx, std::vector<aby3::sbMatrix> &T, std::vector<aby3::si64Matrix> &T_aligned,
        aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void join(int pIdx, std::vector<aby3::si64Matrix> &T_1_key, std::vector<aby3::si64Matrix> &T_1_other,
    std::vector<aby3::si64Matrix> &T_2_key, std::vector<aby3::si64Matrix> &T_2_other,
    std::vector<aby3::si64Matrix> &T_joined,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);


void filter(int pIdx, std::vector<aby3::si64Matrix> &T, int filColIdx, int value, bool is_scalar, std::string op_str,
    std::vector<aby3::si64Matrix> &T_filtered,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void semi_join(int pIdx, std::vector<aby3::si64Matrix> &T_1_key, std::vector<aby3::si64Matrix> &T_1_other,
    std::vector<aby3::si64Matrix> &T_2_key, std::vector<aby3::si64Matrix> &T_2_other,
    std::vector<aby3::si64Matrix> &T_semi_joined,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

// 内积：sum_res = SUM_i (mul_0[i] * mul_1[i])。
// 与 cipher_mul + 累加相比，本函数把元素级乘法份额的 reshare 阶段聚合成一次
// 标量通信，因而通信量从 O(n) 个 i64 降到 O(1) 个 i64。输出 sum_res 为 1×1。
void mul_and_sum(int pIdx, aby3::si64Matrix &mul_0, aby3::si64Matrix &mul_1, aby3::si64Matrix &sum_res,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

#endif