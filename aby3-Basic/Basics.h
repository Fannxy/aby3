#include <cmath>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Network/IOService.h>

#include "../aby3-RTR/debug.h"
#include "assert.h"

#ifndef _ABY3_BASICS_H_
#define _ABY3_BASICS_H_

static int BITSIZE = 64;

inline size_t roundUpToPowerOfTwo(size_t num) {
    if (num == 0) {
        return 1;
    }
    if (checkPowerOfTwo(num)) {
        return num;
    }
    return (size_t) pow(2, ceil(log2(num)));
}

struct boolShare {
    std::array<bool, 2> bshares;

    boolShare(bool share0, bool share1) {
        bshares[0] = share0;
        bshares[1] = share1;
    }

    boolShare() {
        bshares[0] = false;
        bshares[1] = false;
    }

    boolShare(bool plain_val, int pIdx) {
        switch (pIdx) {
            case 0:
                bshares[0] = false;
                bshares[1] = false;
                break;
            case 1:
                bshares[0] = plain_val;
                bshares[1] = false;
                break;
            case 2:
                bshares[0] = false;
                bshares[1] = plain_val;
                break;
            default:
                THROW_RUNTIME_ERROR("boolShare: invalid pIdx");
        }
    }

    aby3::sbMatrix to_matrix() {
        aby3::sbMatrix res(1, 1);
        res.mShares[0](0, 0) = bshares[0];
        res.mShares[1](0, 0) = bshares[1];
        return res;
    }

    void from_matrix(aby3::i64 s1, aby3::i64 s2) {
        bshares[0] = s1 & 1;
        bshares[1] = s2 & 1;
    }
};

struct vecBoolShares {
    std::vector<boolShare> bshares;

    vecBoolShares(std::vector<boolShare> &share_data) { bshares = share_data; }

    vecBoolShares(size_t len) { bshares.resize(len); }

    aby3::sbMatrix to_matrix() {
        size_t len = bshares.size();
        aby3::sbMatrix res(len, 1);
        for (size_t i = 0; i < len; i++) {
            res.mShares[0](i, 0) = bshares[i].bshares[0];
            res.mShares[1](i, 0) = bshares[i].bshares[1];
        }
        return res;
    }

    void from_matrix(aby3::sbMatrix &mat) {
        size_t len = mat.rows();
        bshares.resize(len);
        for (size_t i = 0; i < len; i++) {
            bshares[i].bshares[0] = mat.mShares[0](i, 0) & 1;
            bshares[i].bshares[1] = mat.mShares[1](i, 0) & 1;
        }
    }
};

struct boolIndex {
    std::array<aby3::i64, 2> indexShares;

    boolIndex(aby3::i64 share0, aby3::i64 share1) {
        indexShares[0] = share0;
        indexShares[1] = share1;
    }

    boolIndex() {
        indexShares[0] = 0;
        indexShares[1] = 0;
    }

    boolIndex(aby3::sbMatrix &mat) { from_matrix(mat); }

    boolIndex(aby3::i64 plain_index, int pIdx) {
        switch (pIdx) {
            case 0:
                indexShares[0] = 0;
                indexShares[1] = 0;
                break;
            case 1:
                indexShares[0] = plain_index;
                indexShares[1] = 0;
                break;
            case 2:
                indexShares[0] = 0;
                indexShares[1] = plain_index;
                break;
            default:
                throw std::runtime_error("boolIndex: invalid pIdx");
        }
    }

    aby3::sbMatrix to_matrix() {
        aby3::sbMatrix res(1, BITSIZE);
        res.mShares[0](0, 0) = indexShares[0];
        res.mShares[1](0, 0) = indexShares[1];
        return res;
    }

    void from_matrix(aby3::sbMatrix &m) {
        assert(m.mShares[0].rows() == 1 && m.mShares[0].cols() == 1);
        indexShares[0] = m.mShares[0](0, 0);
        indexShares[1] = m.mShares[1](0, 0);
    }

    operator aby3::sbMatrix() const {
        aby3::sbMatrix res(1, BITSIZE);
        res.mShares[0](0, 0) = indexShares[0];
        res.mShares[1](0, 0) = indexShares[1];
        return res;
    }
};

struct vecBoolIndices {
    std::vector<boolIndex> indexShares;

    vecBoolIndices(std::vector<boolIndex> &share_data) {
        indexShares = share_data;
    }

    vecBoolIndices(aby3::sbMatrix &mat) { from_matrix(mat); }

    vecBoolIndices(size_t len) { indexShares.resize(len); }

    vecBoolIndices() {}

    aby3::sbMatrix to_matrix() {
        size_t len = indexShares.size();
        aby3::sbMatrix res(len, BITSIZE);
        for (size_t i = 0; i < len; i++) {
            res.mShares[0](i, 0) = indexShares[i].indexShares[0];
            res.mShares[1](i, 0) = indexShares[i].indexShares[1];
        }
        return res;
    }

    void from_matrix(aby3::sbMatrix &mat) {
        size_t len = mat.rows();
        indexShares.resize(len);
        for (size_t i = 0; i < len; i++) {
            indexShares[i].indexShares[0] = mat.mShares[0](i, 0);
            indexShares[i].indexShares[1] = mat.mShares[1](i, 0);
        }
    }
};

// if toNext is true, send to the next party, otherwise send to the previous party.
int large_data_sending(int pIdx, aby3::i64Matrix &sharedA, aby3::Sh3Runtime &runtime, bool toNext);

// if fromPrev is true, receive from the previous party, otherwise receive from the next party.
int large_data_receiving(int pIdx, aby3::i64Matrix &res, aby3::Sh3Runtime &runtime, bool fromPrev);

int large_data_encryption(int pIdx, aby3::i64Matrix &plainA, aby3::sbMatrix &sharedA,
                          aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);

int large_data_decryption(int pIdx, aby3::sbMatrix &sharedA, aby3::i64Matrix &plainA,
                          aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);

void bool_cipher_lt(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                    aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                    aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_eq(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                    aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                    aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_eq(int pIdx, aby3::sbMatrix &sharedA, aby3::i64Matrix &plainB,
                    aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                    aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_or(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                    aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                    aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_or(int pIdx, boolShare &sharedA, boolShare &sharedB,
                    boolShare &res, aby3::Sh3Encryptor &enc,
                    aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_add(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_sub(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_and(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_and(int pIdx, boolShare &sharedA, boolShare &sharedB,
                     boolShare &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_max(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_min(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_max_min_split(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res_max, aby3::sbMatrix &res_min, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_not(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &res);

void bool_cipher_not(int pIdx, std::vector<boolShare> &sharedA,
                     std::vector<boolShare> &res);

void bool_cipher_not(int pIdx, boolShare &sharedA, boolShare &res);

void bool_cipher_dot(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

void bool_cipher_dot(int pIdx, std::vector<aby3::sbMatrix> &sharedA,
                     aby3::sbMatrix &sharedB, aby3::sbMatrix &res,
                     aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                     aby3::Sh3Runtime &runtime);

void bool_cipher_selector(int pIdx, boolShare &flag, aby3::sbMatrix &trueVal,
                          aby3::sbMatrix &falseVal, aby3::sbMatrix &res,
                          aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                          aby3::Sh3Runtime &runtime);

void bool2arith(int pIdx, aby3::sbMatrix &boolInput, aby3::si64Matrix &res,
                aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                aby3::Sh3Runtime &runtime);

void bool_get_first_zero_mask(int pIdx, std::vector<boolShare> &inputA,
                              aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                              aby3::Sh3Evaluator &eval,
                              aby3::Sh3Runtime &runtime);

void bool_init_false(int pIdx, aby3::sbMatrix &res);

void bool_init_true(int pIdx, aby3::sbMatrix &res);

void bool_init_false(int pIdx, boolShare &res);

void bool_init_true(int pIdx, boolShare &res);

void bool_shift(int pIdx, boolIndex &sharedA, size_t shift_len,
                boolIndex &res_shift, bool right_flag = true);

void bool_shift_and_left(int pIdx, aby3::sbMatrix &sharedA, size_t shift_len,
                         aby3::sbMatrix &res_shift, aby3::sbMatrix &res_left);

void bool_shift_and_left(int pIdx, boolIndex &sharedA, size_t shift_len,
                         boolIndex &res_shift, boolIndex &res_left);


// TODO: support various functions.
void bool_aggregation(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &res,
                         aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                         aby3::Sh3Runtime &runtime, const std::string& func);

aby3::i64Matrix back2plain(int pIdx, aby3::sbMatrix &cipher_val,
                           aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                           aby3::Sh3Runtime &runtime);

bool back2plain(int pIdx, boolShare &cipher_val, aby3::Sh3Encryptor &enc,
                aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

aby3::i64 back2plain(int pIdx, boolIndex &cipher_val, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

std::vector<bool> back2plain(int pIdx, std::vector<boolShare> &cipher_val,
                             aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                             aby3::Sh3Runtime &runtime);

aby3::i64Matrix back2plain(int pIdx, std::vector<aby3::si64>& cipher_val, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                           aby3::Sh3Runtime &runtime);

void get_permutation(size_t len, std::vector<size_t> &permutation,
                     oc::block &seed);

void get_inverse_permutation(std::vector<size_t> &permutation,
                             std::vector<size_t> &inverse_permutation);

void combine_permutation(std::vector<std::vector<size_t>> &permutation_list,
                         std::vector<size_t> &final_permutation);

template <typename T>
void plain_permutate(std::vector<size_t> &permutation, std::vector<T> &data) {
    size_t len = data.size();
    std::vector<T> tmp(len);
    for (size_t i = 0; i < len; i++) {
        tmp[permutation[i]] = data[i];
    }
    data = tmp;
}

void plain_permutate(std::vector<size_t> &permutation, aby3::sbMatrix &data);

void plain_permutate(std::vector<size_t> &permutation, aby3::i64Matrix &data);

void get_random_mask(int pIdx, aby3::i64Matrix &res, oc::block &seed);

void arith_aggregation(int pIdx, aby3::si64Matrix &sharedA, aby3::si64Matrix &res,
                         aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                         aby3::Sh3Runtime &runtime, const std::string& func);

template <typename T>
void vector_tile(const std::vector<T>& v, int n, std::vector<T>& target, int start){
    for(int i=0; i<n; i++){
        std::copy(v.begin(), v.end(), target.begin() + start + (i * v.size()));
    }
    return;
}

template <typename T>
void vector_repeat(const std::vector<T>& v, int n, std::vector<T> &target, int start){
    for(int i=0; i<v.size(); i++){
        std::fill(target.begin() + start + i*n, target.begin() + start + (i + 1) * n, v[i]);
    }
    return;
}

std::vector<size_t> argwhere(aby3::i64Matrix& input, int target);

std::vector<size_t> argwhere(std::vector<std::vector<int>>& inputs, int target);

std::vector<size_t> argwhere(std::vector<size_t>& inputs, int target);

void left_shift_and_fill(int pIdx, aby3::sbMatrix &input, int tag_size, int tag_value);

void tag_append(int pIdx, std::vector<aby3::sbMatrix>& inputs);

void tag_remove(int pIdx, size_t tag_len, std::vector<aby3::sbMatrix>& inputs);

#endif