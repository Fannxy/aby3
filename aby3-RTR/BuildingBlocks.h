#include <cryptoTools/Network/IOService.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include "debug.h"

#ifndef _A_H_
#define _A_H_

/// @brief generate the synchronized time using the one in role 0.
/// @param pIdx 
/// @param time_slot 
/// @param runtime 
/// @return 
double synchronized_time(int pIdx, double& time_slot, aby3::Sh3Runtime &runtime);

/// setup functions.
void distribute_setup(aby3::u64 partyIdx, oc::IOService &ios, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
           aby3::Sh3Runtime &runtime);

// local setup funtion, each party is assigned to a thread.
void basic_setup(aby3::u64 partyIdx, oc::IOService &ios, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
           aby3::Sh3Runtime &runtime);

void multi_processor_setup(aby3::u64 partyIdx, int rank, oc::IOService &ios, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void splitted_setup(aby3::u64 partyIdx, int rank, oc::IOService &ios, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime, std::string p0_ip, std::string p1_ip);


/// basic building block functions.
/// all functions for aby3 data structures are : cipher_func
/// all functions for vectors are : vector_cipher_func

/// 1. multiplitcations, totally, we support: 1) plain_a and shared_b;) 2_shared_a and shared_b; 3) shared_a and shared_a; 4) shared_f and shared_b; 5) shared_f and shared_f.
int pi_cb_mul(int pIdx, const aby3::i64Matrix &plainA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime &runtime);


int cipher_mul_seq(int pIdx, const aby3::si64Matrix &sharedA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);


int cipher_mul_seq(int pIdx, const aby3::si64Matrix &sharedA, const aby3::si64Matrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);


template <aby3::Decimal D>
int cipher_mul_seq(int pIdx, const aby3::sf64Matrix<D> &sharedA, const aby3::sbMatrix &sharedB, aby3::sf64Matrix<D> &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime){
    return cipher_mul_seq(pIdx, sharedA.i64Cast(), sharedB, res.i64Cast(), eval, enc, runtime);
}


template <aby3::Decimal D>
int cipher_mul_seq(int pIdx, aby3::sf64Matrix<D> &sharedA, aby3::sf64Matrix<D> &sharedB, aby3::sf64Matrix<D> &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime){
    eval.asyncMul(runtime, sharedA, sharedB, res).get();
    return 0;
}

inline int cipher_mul_seq(int pIdx, const aby3::i64Matrix &plainA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime &runtime){
    return pi_cb_mul(pIdx, plainA, sharedB, res, eval, enc, runtime);
}

// standard interfaces of ciphertext multiplications.
int cipher_mul(int pIdx, const aby3::i64Matrix &plainA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime &runtime);


int cipher_mul(int pIdx, const aby3::si64Matrix &sharedA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);


int cipher_mul(int pIdx, const aby3::si64Matrix &sharedA, const aby3::si64Matrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);


template <aby3::Decimal D>
int cipher_mul(int pIdx, const aby3::sf64Matrix<D> &sharedA, const aby3::sbMatrix &sharedB, aby3::sf64Matrix<D> &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime){
    return cipher_mul(pIdx, sharedA.i64Cast(), sharedB, res.i64Cast(), eval, enc, runtime);
}


template <aby3::Decimal D>
int cipher_mul(int pIdx, aby3::sf64Matrix<D> &sharedA, aby3::sf64Matrix<D> &sharedB, aby3::sf64Matrix<D> &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime){
    eval.asyncMul(runtime, sharedA, sharedB, res).get();
    return 0;
}

int vec2mat(std::vector<aby3::si64>& vec, aby3::si64Matrix& mat);

int mat2vec(aby3::si64Matrix& mat, std::vector<aby3::si64>& vec);

// synchronized version of fetch_msb.
int fetch_msb(int pIdx, aby3::si64Matrix &diffAB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime, aby3::Sh3Task &task);

// asynchronized version of fetch_msb.
int fetch_msb(int pIdx, aby3::si64Matrix &diffAB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// sint greater-than.
int cipher_gt(int pIdx, aby3::si64Matrix &sharedA, aby3::si64Matrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// sint and plaintext greater-than.
int cipher_gt(int pIdx, aby3::si64Matrix &sharedA, std::vector<int> &plainB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// sfixed greater-than.
template <aby3::Decimal D>
int cipher_gt(int pIdx, aby3::sf64Matrix<D> &sharedA, aby3::sf64Matrix<D> &sharedB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
    return cipher_gt(pIdx, sharedA.i64Cast(), sharedB.i64Cast(), res, eval, runtime);
}

int vector_cipher_gt(int pIdx, std::vector<aby3::si64>& sintA, std::vector<aby3::si64>& sintB, aby3::sbMatrix& res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);

// vector greater-than.
int vector_cipher_gt(int pIdx, std::vector<aby3::si64>& sintA, std::vector<aby3::si64>& sintB, std::vector<aby3::si64>& res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);

int vector_cipher_ge(int pIdx, std::vector<aby3::si64>& sintA, std::vector<aby3::si64>& sintB, aby3::sbMatrix& res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);

int cipher_ge(int pIdx, aby3::si64Matrix& sintA, aby3::si64Matrix& sintB, aby3::sbMatrix& res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);

// eq implemeneted through two fetch_msb and an Nor gate.
int cipher_eq(int pIdx, aby3::si64Matrix &intA, aby3::si64Matrix &intB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

template <aby3::Decimal D>
int cipher_eq(int pIdx, aby3::sf64Matrix<D> &sharedA, aby3::sf64Matrix<D> &sharedB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
    return cipher_eq(pIdx, sharedA.i64Cast(), sharedB.i64Cast(), res, eval, runtime);
}

// eq implemented through int_eq circuit.
int circuit_cipher_eq(int pIdx, aby3::si64Matrix &intA, aby3::si64Matrix &intB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

int vector_cipher_eq(int pIdx, std::vector<aby3::si64>& intA, std::vector<int>& intB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

int vector_cipher_eq(int pIdx, std::vector<int>& intA, std::vector<aby3::si64>& intB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// fetch_eq_res.
int fetch_eq_res(int pIdx, aby3::sbMatrix& circuitA, aby3::sbMatrix& circuitB, aby3::sbMatrix& res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

int vector_mean_square(int pIdx, const std::vector<aby3::si64>&sharedA, const std::vector<aby3::si64>&sharedB, std::vector<aby3::si64>&res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime &runtime);

// initiate functions.
int init_ones(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, aby3::si64Matrix &res, int n);

int init_zeros(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, aby3::si64Matrix &res, int n);

int vector_generation(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, std::vector<aby3::si64>& vecRes);

int vector_generation(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, std::vector<int>& vecRes);

int vector_generation(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, std::vector<std::vector<aby3::si64>>& vecRes);

int plain_mat_to_cipher_mat(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, aby3::i64Matrix& plain_mat, aby3::si64Matrix& cipher_mat);

// basic encryption methods.
int vector_to_matrix(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, std::vector<int>& plain_vec, aby3::si64Matrix& enc_mat);

template <aby3::Decimal D>
int vector_to_matrix(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, std::vector<double>& plain_vec, aby3::sf64Matrix<D>& enc_mat){
    aby3::u64 vec_len = plain_vec.size();
    aby3::f64Matrix<D> plain_mat; plain_mat.resize(vec_len, 1);
    for(int i=0; i<vec_len; i++) plain_mat(i, 0) = plain_vec[i];
    plain_mat_to_cipher_mat(pIdx, enc, runtime, plain_mat.i64Cast(), enc_mat.i64Cast());
    return 0;
}

int cipher_matrix_to_vector(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, aby3::si64Matrix& cipher_mat, std::vector<aby3::si64>& cipher_vec);

template <aby3::Decimal D>
int cipher_matrix_to_vector(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, aby3::sf64Matrix<D>& cipher_mat, std::vector<aby3::sf64<D>>& cipher_vec){
    size_t rows = cipher_mat.rows();
    size_t vec_len = cipher_vec.size();
    if(vec_len == rows){
        for(int i=0; i<vec_len; i++) cipher_vec[i] = cipher_mat(i, 0);
    }
    else if (vec_len == 0){
        for(int i=0; i<rows; i++) cipher_vec.push_back(cipher_mat(i, 0));
    }
    else{
        std::cerr << "not null vector" << std::endl;
        return 1;
    }
    return 0;
}

int vector_to_cipher_vector(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, std::vector<int>& plain_vec, std::vector<aby3::si64>& enc_vec);

template <aby3::Decimal D>
int vector_to_cipher_vector(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, std::vector<double>& plain_vec, std::vector<aby3::sf64<D>>& enc_vec){
    aby3::sf64Matrix<D> enc_mat;
    vector_to_matrix<D>(pIdx, enc, runtime, plain_vec, enc_mat);
    return cipher_matrix_to_vector<D>(pIdx, enc, runtime, enc_mat, enc_vec);
}

template <aby3::Decimal D>
int init_zeros(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, aby3::sf64Matrix<D> &res, int n){
    aby3::f64Matrix<D> plain_mat; plain_mat.resize(n, 1);
    for(int i=0; i<n; i++) plain_mat(i, 0) = 0;
    aby3::sf64Matrix<D> enc_mat;
    plain_mat_to_cipher_mat(pIdx, enc, runtime, plain_mat.i64Cast(), enc_mat.i64Cast());
}

template <aby3::Decimal D>
int init_zeros(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, std::vector<aby3::sf64<D>> &res, int n){
    aby3::f64Matrix<D> plain_mat; plain_mat.resize(n, 1);
    for(int i=0; i<n; i++) plain_mat(i, 0) = 0;
    aby3::sf64Matrix<D> enc_mat;
    plain_mat_to_cipher_mat(pIdx, enc, runtime, plain_mat.i64Cast(), enc_mat.i64Cast());
    if(res.size() == n){
        for(int i=0; i<n; i++) res[i] = enc_mat(i, 0);
    }
    else if(res.size() == 0){
        for(int i=0; i<n; i++) res.push_back(enc_mat(i, 0));
    }
    return 0;
}


#endif