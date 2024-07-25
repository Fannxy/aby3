#include <aby3/Circuit/CircuitLibrary.h>
#include <aby3/sh3/Sh3BinaryEvaluator.h>

#include <bitset>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <cmath>

#include "debug.h"
#include "BuildingBlocks.h"
#include "Basics.h"

#define DEBUG_BASIC

using namespace oc;
using namespace aby3;

void bool_cipher_lt(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB,
                    sbMatrix &res, Sh3Encryptor &enc, Sh3Evaluator &eval,
                    Sh3Runtime &runtime) {
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_lt(bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedB);
    binEng.setInput(1, sharedA);

    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, 1);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_eq(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB,
                    sbMatrix &res, Sh3Encryptor &enc, Sh3Evaluator &eval,
                    Sh3Runtime &runtime) {
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_eq(bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);

    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, 1);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_eq(int pIdx, aby3::sbMatrix &sharedA, aby3::i64Matrix &plainB,
                    aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                    aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime) {
    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    aby3::sbMatrix sharedB(i64Size, bitSize);

    switch (pIdx) {
        case 0:
            for (int i = 0; i < i64Size; i++) {
                sharedB.mShares[0](i, 0) = 0;
                sharedB.mShares[1](i, 0) = 0;
            }
            break;
        case 1:
            for (int i = 0; i < i64Size; i++) {
                sharedB.mShares[0](i, 0) = plainB(i, 0);
                sharedB.mShares[1](i, 0) = 0;
            }
            break;
        case 2:
            for (int i = 0; i < i64Size; i++) {
                sharedB.mShares[0](i, 0) = 0;
                sharedB.mShares[1](i, 0) = plainB(i, 0);
            }
            break;
        default:
            // throw std::runtime_error("Error in file " + std::string(__FILE__)
            // + " at line " + std::to_string(__LINE__) + " bool_cipher_eq: pIdx
            // out of range.");
            THROW_RUNTIME_ERROR("bool_cipher_eq: pIdx out of range.");
    }

    bool_cipher_eq(pIdx, sharedA, sharedB, res, enc, eval, runtime);
    return;
}

void bool_cipher_or(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB,
                    sbMatrix &res, Sh3Encryptor &enc, Sh3Evaluator &eval,
                    Sh3Runtime &runtime) {
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_bitwiseOr(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);

    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_or(int pIdx, boolShare &sharedA, boolShare &sharedB,
                    boolShare &res, aby3::Sh3Encryptor &enc,
                    aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime) {
    bool cross1 = sharedA.bshares[0] && sharedB.bshares[0];
    bool cross2 = sharedA.bshares[0] && sharedB.bshares[1];
    bool cross3 = sharedA.bshares[1] && sharedB.bshares[0];
    bool share = (cross1 ^ cross2 ^ cross3);

    share = share ^ sharedA.bshares[0] ^ sharedB.bshares[0];

    runtime.mComm.mNext.asyncSendCopy(share);
    bool other_share;
    runtime.mComm.mPrev.recv(other_share);

    res.bshares = {share, other_share};

    return;
}

void bool_cipher_add(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB,
                     sbMatrix &res, Sh3Encryptor &enc, Sh3Evaluator &eval,
                     Sh3Runtime &runtime) {
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_add(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);

    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_sub(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB,
                     sbMatrix &res, Sh3Encryptor &enc, Sh3Evaluator &eval,
                     Sh3Runtime &runtime) {
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_subtract(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);

    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_and(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime) {
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_bitwiseAnd(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);

    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();

    return;
}

void bool_cipher_and(int pIdx, boolShare &sharedA, boolShare &sharedB,
                     boolShare &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime) {
    bool cross1 = sharedA.bshares[0] && sharedB.bshares[0];
    bool cross2 = sharedA.bshares[0] && sharedB.bshares[1];
    bool cross3 = sharedA.bshares[1] && sharedB.bshares[0];
    bool share = (cross1 ^ cross2 ^ cross3);

    runtime.mComm.mNext.asyncSendCopy(share);
    bool other_share;
    runtime.mComm.mPrev.recv(other_share);

    res.bshares = {share, other_share};

    return;
}

void bool_cipher_max(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
    
    aby3::sbMatrix comp;
    bool_cipher_lt(pIdx, sharedA, sharedB, comp, enc, eval, runtime);
    aby3::sbMatrix tmp1, tmp2;

    comp.resize(comp.rows(), BITSIZE);
    for(size_t i=0; i<comp.rows(); i++){
        comp.mShares[0](i, 0) = (comp.mShares[0](i, 0) == 1) ? -1 : -0;
        comp.mShares[1](i, 0) = (comp.mShares[1](i, 0) == 1) ? -1 : -0;
    }
    bool_cipher_and(pIdx, comp, sharedB, tmp1, enc, eval, runtime);
    bool_cipher_not(pIdx, comp, comp);
    bool_cipher_and(pIdx, comp, sharedA, tmp2, enc, eval, runtime);
    res.resize(sharedA.rows(), BITSIZE);
    for(size_t i=0; i<tmp1.rows(); i++){
        res.mShares[0](i, 0) = tmp1.mShares[0](i, 0) ^ tmp2.mShares[0](i, 0);
        res.mShares[1](i, 0) = tmp1.mShares[1](i, 0) ^ tmp2.mShares[1](i, 0);
    }
    return;
}

void bool_cipher_min(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
    aby3::sbMatrix comp;
    bool_cipher_lt(pIdx, sharedA, sharedB, comp, enc, eval, runtime);
    aby3::sbMatrix tmp1, tmp2;

    comp.resize(comp.rows(), BITSIZE);
    for(size_t i=0; i<comp.rows(); i++){
        comp.mShares[0](i, 0) = (comp.mShares[0](i, 0) == 1) ? -1 : -0;
        comp.mShares[1](i, 0) = (comp.mShares[1](i, 0) == 1) ? -1 : -0;
    }
    bool_cipher_and(pIdx, comp, sharedA, tmp1, enc, eval, runtime);
    bool_cipher_not(pIdx, comp, comp);
    bool_cipher_and(pIdx, comp, sharedB, tmp2, enc, eval, runtime);
    res.resize(sharedA.rows(), BITSIZE);
    for(size_t i=0; i<tmp1.rows(); i++){
        res.mShares[0](i, 0) = tmp1.mShares[0](i, 0) ^ tmp2.mShares[0](i, 0);
        res.mShares[1](i, 0) = tmp1.mShares[1](i, 0) ^ tmp2.mShares[1](i, 0);
    }
    return;
}

void bool_cipher_max_min_split(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res_max, aby3::sbMatrix &res_min, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
    aby3::sbMatrix comp;
    bool_cipher_lt(pIdx, sharedA, sharedB, comp, enc, eval, runtime);

    aby3::sbMatrix extend_comp(sharedA.rows()*2, BITSIZE);
    for(size_t i=0; i<sharedA.rows(); i++){
        extend_comp.mShares[0](i, 0) = (comp.mShares[0](i, 0) == 1) ? -1 : 0;
        extend_comp.mShares[1](i, 0) = (comp.mShares[1](i, 0) == 1) ? -1 : 0;
        extend_comp.mShares[0](i+sharedA.rows(), 0) = (comp.mShares[0](i, 0) == 1) ? -1 : 0;
        extend_comp.mShares[1](i+sharedA.rows(), 0) = (comp.mShares[1](i, 0) == 1) ? -1 : 0;
    }

    aby3::sbMatrix extend_AB(sharedA.rows()*2, BITSIZE);
    for(size_t i=0; i<sharedA.rows(); i++){
        extend_AB.mShares[0](i, 0) = sharedA.mShares[0](i, 0);
        extend_AB.mShares[1](i, 0) = sharedA.mShares[1](i, 0);
        extend_AB.mShares[0](i+sharedA.rows(), 0) = sharedB.mShares[0](i, 0);
        extend_AB.mShares[1](i+sharedA.rows(), 0) = sharedB.mShares[1](i, 0);
    }

    aby3::sbMatrix tmp1, tmp2;
    bool_cipher_and(pIdx, extend_comp, extend_AB, tmp1, enc, eval, runtime);
    bool_cipher_not(pIdx, extend_comp, extend_comp);
    bool_cipher_and(pIdx, extend_comp, extend_AB, tmp2, enc, eval, runtime);

    res_max.resize(sharedA.rows(), BITSIZE);
    res_min.resize(sharedA.rows(), BITSIZE);
    for(size_t i=0; i<sharedA.rows(); i++){
        res_min.mShares[0](i, 0) = tmp1.mShares[0](i, 0) ^ tmp2.mShares[0](i + sharedA.rows(), 0);
        res_min.mShares[1](i, 0) = tmp1.mShares[1](i, 0) ^ tmp2.mShares[1](i + sharedA.rows(), 0);
        res_max.mShares[0](i, 0) = tmp1.mShares[0](i + sharedA.rows(), 0) ^ tmp2.mShares[0](i, 0);
        res_max.mShares[1](i, 0) = tmp1.mShares[1](i + sharedA.rows(), 0) ^ tmp2.mShares[1](i, 0);
    }

    return;
}


void bool_cipher_not(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &res) {
    int i64Size = sharedA.i64Size();

    switch (pIdx) {
        case 0:
            for (size_t i = 0; i < i64Size; i++) {
                res.mShares[0](i, 0) = sharedA.mShares[0](i, 0);
                res.mShares[1](i, 0) = sharedA.mShares[1](i, 0);
            }
            break;
        case 1:
            for (size_t i = 0; i < i64Size; i++) {
                res.mShares[0](i, 0) = ~sharedA.mShares[0](i, 0);
                res.mShares[1](i, 0) = sharedA.mShares[1](i, 0);
            }
            break;
        case 2:
            for (size_t i = 0; i < i64Size; i++) {
                res.mShares[0](i, 0) = sharedA.mShares[0](i, 0);
                res.mShares[1](i, 0) = ~sharedA.mShares[1](i, 0);
            }
            break;
        default:
            throw std::runtime_error("bool_cipher_not: pIdx out of range.");
    }

    return;
}

void bool_cipher_not(int pIdx, std::vector<boolShare> &sharedA,
                     std::vector<boolShare> &res) {
    if (res.size() != sharedA.size()) res.resize(sharedA.size());
    size_t len = sharedA.size();
    switch (pIdx) {
        case 0:
            for (size_t i = 0; i < len; i++) {
                res[i].bshares[0] = sharedA[i].bshares[0];
                res[i].bshares[1] = sharedA[i].bshares[1];
            }
            break;
        case 1:
            for (size_t i = 0; i < len; i++) {
                res[i].bshares[0] = !sharedA[i].bshares[0];
                res[i].bshares[1] = sharedA[i].bshares[1];
            }
            break;
        case 2:
            for (size_t i = 0; i < len; i++) {
                res[i].bshares[0] = sharedA[i].bshares[0];
                res[i].bshares[1] = !sharedA[i].bshares[1];
            }
            break;
        default:
            throw std::runtime_error("bool_cipher_not: pIdx out of range.");
    }
    return;
}

void bool_cipher_not(int pIdx, boolShare &sharedA, boolShare &res) {
    switch (pIdx) {
        case 0:
            res.bshares[0] = sharedA.bshares[0];
            res.bshares[1] = sharedA.bshares[1];
            break;
        case 1:
            res.bshares[0] = !sharedA.bshares[0];
            res.bshares[1] = sharedA.bshares[1];
            break;
        case 2:
            res.bshares[0] = sharedA.bshares[0];
            res.bshares[1] = !sharedA.bshares[1];
            break;
        default:
            throw std::runtime_error("bool_cipher_not: pIdx out of range.");
    }
    return;
}

void bool_cipher_dot(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB,
                     aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime) {
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    aby3::sbMatrix mul_res(i64Size, bitSize);

    auto cir = lib.int_int_bitwiseAnd(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);

    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        mul_res.resize(i64Size, bitSize);
        binEng.getOutput(0, mul_res);
    });
    dep.get();

    res.resize(1, bitSize);
    res.mShares[0](0, 0) = mul_res.mShares[0](0, 0);
    res.mShares[1](0, 0) = mul_res.mShares[1](0, 0);
    for (size_t i = 1; i < i64Size; i++) {
        res.mShares[0](0, 0) ^= mul_res.mShares[0](i, 0);
        res.mShares[1](0, 0) ^= mul_res.mShares[1](i, 0);
    }
}

void bool_cipher_selector(int pIdx, boolShare &flag, aby3::sbMatrix &trueVal,
                          aby3::sbMatrix &falseVal, aby3::sbMatrix &res,
                          aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                          aby3::Sh3Runtime &runtime) {
    int bitSize = trueVal.bitCount();
    int i64Size = trueVal.i64Size();

    if (bitSize != 64) {
        THROW_RUNTIME_ERROR("The bitsize must be 64!");
    }

    boolShare nFlag;
    bool_cipher_not(pIdx, flag, nFlag);

    aby3::sbMatrix _flagMat(i64Size, 1);
    aby3::sbMatrix _nFlagMat(i64Size, 1);
    _flagMat.resize(i64Size, bitSize);
    _nFlagMat.resize(i64Size, bitSize);

    for (int i = 0; i < i64Size; i++) {
        _flagMat.mShares[0](i, 0) = flag.bshares[0] ? -1 : -0;
        _flagMat.mShares[1](i, 0) = flag.bshares[1] ? -1 : -0;
        _nFlagMat.mShares[0](i, 0) = nFlag.bshares[0] ? -1 : -0;
        _nFlagMat.mShares[1](i, 0) = nFlag.bshares[1] ? -1 : -0;
    }
    bool_cipher_and(pIdx, _flagMat, trueVal, _flagMat, enc, eval, runtime);
    bool_cipher_and(pIdx, _nFlagMat, falseVal, _nFlagMat, enc, eval, runtime);

    for (int i = 0; i < i64Size; i++) {
        res.mShares[0](i, 0) =
            _flagMat.mShares[0](i, 0) ^ _nFlagMat.mShares[0](i, 0);
        res.mShares[1](i, 0) =
            _flagMat.mShares[1](i, 0) ^ _nFlagMat.mShares[1](i, 0);
    }

    return;
}

void bool_cipher_dot(int pIdx, std::vector<aby3::sbMatrix> &sharedA,
                     aby3::sbMatrix &sharedB, aby3::sbMatrix &res,
                     aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                     aby3::Sh3Runtime &runtime) {
    size_t n = sharedA.size();
    if (n != sharedB.rows()) {
        THROW_RUNTIME_ERROR("The size of sharedA and sharedB does not match!");
    }

    size_t block_len = sharedA[0].i64Size();
    size_t bitSize = sharedA[0].bitCount();

    // expand the bitsize of sharedB to match with the sharedA.
    if (sharedB.bitCount() == 1) {
        for (size_t i = 0; i < n; i++) {
            sharedB.mShares[0](i, 0) =
                (sharedB.mShares[0](i, 0) == 1) ? -1 : -0;
            sharedB.mShares[1](i, 0) =
                (sharedB.mShares[1](i, 0) == 1) ? -1 : -0;
        }
    }

    // expand sharedB and sharedA to match the size.
    aby3::sbMatrix _expandA(n * block_len, bitSize);
    aby3::sbMatrix _expandB(n * block_len, bitSize);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < block_len; j++) {
            _expandA.mShares[0](i * block_len + j, 0) =
                sharedA[i].mShares[0](j, 0);
            _expandA.mShares[1](i * block_len + j, 0) =
                sharedA[i].mShares[1](j, 0);
            _expandB.mShares[0](i * block_len + j, 0) =
                sharedB.mShares[0](i, 0);
            _expandB.mShares[1](i * block_len + j, 0) =
                sharedB.mShares[1](i, 0);
        }
    }

    // compute the element-wise AND product.
    bool_cipher_and(pIdx, _expandA, _expandB, _expandA, enc, eval, runtime);

    // add the elements for the final result.
    res.resize(block_len, bitSize);
    for (size_t i = 0; i < block_len; i++) {
        res.mShares[0](i, 0) = _expandA.mShares[0](i, 0);
        res.mShares[1](i, 0) = _expandA.mShares[1](i, 0);
        for (size_t j = 1; j < n; j++) {
            res.mShares[0](i, 0) ^= _expandA.mShares[0](j * block_len + i, 0);
            res.mShares[1](i, 0) ^= _expandA.mShares[1](j * block_len + i, 0);
        }
    }
    return;
}

void bool2arith(int pIdx, aby3::sbMatrix &boolInput, aby3::si64Matrix &res,
                aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                aby3::Sh3Runtime &runtime){
    size_t len = boolInput.rows();
    size_t bitSize = boolInput.bitCount();

    if(bitSize == 1){
        aby3::i64Matrix plainInput(len, 1);
        for(size_t i=0; i<len; i++) plainInput(i, 0) = 1;
        pi_cb_mul(pIdx, plainInput, boolInput, res, eval, enc, runtime);
    }
    else{
        // 1) p0 and p1 generate the correlated random r.
        aby3::sbMatrix random_mat(len, bitSize);
        if(pIdx < 2){
            aby3::i64Matrix random_r(len, 1);
            if(pIdx == 0){
                block seed = enc.mShareGen.mNextCommon.getSeed();
                PRNG prng(seed);
                for(size_t i=0; i<len; i++) random_r(i, 0) = (i64)prng.get<int32_t>();
                std::memcpy(random_mat.mShares[0].data(), random_r.data(), len*sizeof(random_r(0, 0)));
                std::fill_n(random_mat.mShares[1].data(), len, 0);
            }
            else{
                block seed = enc.mShareGen.mPrevCommon.getSeed();
                PRNG prng(seed);
                for(size_t i=0; i<len; i++) random_r(i, 0) = (i64)prng.get<int32_t>();
                std::memcpy(random_mat.mShares[1].data(), random_r.data(), len*sizeof(random_r(0, 0)));
                std::fill_n(random_mat.mShares[0].data(), len, 0);
            }
        }
        else{
            std::fill_n(random_mat.mShares[0].data(), len, 0);
            std::fill_n(random_mat.mShares[1].data(), len, 0);
        }

        // 2) compute c = x ADD r.
        aby3::sbMatrix cipher_mat(len, bitSize);
        bool_cipher_add(pIdx, boolInput, random_mat, cipher_mat, enc, eval, runtime);

        // 3) reveal c to P2.
        aby3::si64Matrix arithRes(len, 1);
        if(pIdx == 0){
            for(size_t i=0; i<len; i++) arithRes.mShares[0](i, 0) = -1 * random_mat.mShares[0](i, 0);
            large_data_receiving(pIdx, arithRes.mShares[1], runtime, true);
        }
        if(pIdx == 1){
            aby3::i64Matrix share1(len, 1);
            std::memcpy(share1.data(), cipher_mat.mShares[1].data(), len*sizeof(cipher_mat.mShares[0](0, 0)));
            large_data_sending(pIdx, share1, runtime, true);
            for(size_t i=0; i<len; i++){
                arithRes.mShares[1](i, 0) = -1 * random_mat.mShares[1](i, 0);
            }
            large_data_receiving(pIdx, arithRes.mShares[0], runtime, false);
        }
        if(pIdx == 2){
            aby3::i64Matrix share1(len, 1);
            large_data_receiving(pIdx, share1, runtime, true);
            aby3::i64Matrix plain_c(len, 1);
            for(size_t i=0; i<len; i++){
                plain_c(i, 0) = (cipher_mat.mShares[0](i, 0) ^ cipher_mat.mShares[1](i, 0) ^ share1(i, 0));
            }
            // construct the new shares of the result.
            for(size_t i=0; i<len; i++){
                share1(i, 0) = rand();
                arithRes.mShares[0](i, 0) = plain_c(i, 0) - share1(i, 0);
                arithRes.mShares[1](i, 0) = share1(i, 0);
            }
            large_data_sending(pIdx, arithRes.mShares[0], runtime, true);
            large_data_sending(pIdx, arithRes.mShares[1], runtime, false);
        }

        res = arithRes;
    }

    return;
}

// log(n) rounds with o(nlogn) computation method.
void bool_get_first_zero_mask(int pIdx, std::vector<boolShare> &inputA,
                              aby3::sbMatrix &res, aby3::Sh3Encryptor &enc,
                              aby3::Sh3Evaluator &eval,
                              aby3::Sh3Runtime &runtime) {
    size_t len = inputA.size();
    size_t round = (size_t)floor(log2(len));

    std::vector<boolShare> tmp_notA(len);
    bool_cipher_not(pIdx, inputA, tmp_notA);
    std::vector<boolShare> _init_mask(tmp_notA.begin(), tmp_notA.end());
    std::rotate(_init_mask.rbegin(), _init_mask.rbegin() + 1,
                _init_mask.rend());
    _init_mask[0] = boolShare(0, 0);

    sbMatrix init_mask = vecBoolShares(_init_mask).to_matrix();
    // log-n round prefix summation logic.
    for (size_t i = 0; i < round; i++) {
        size_t stride = (size_t)pow(2, i);
        sbMatrix _init_maskA(len - stride, 1), _init_maskB(len - stride, 1);

        for (size_t j = stride; j < len; j++) {
            _init_maskA.mShares[0](j - stride, 0) = init_mask.mShares[0](j, 0);
            _init_maskA.mShares[1](j - stride, 0) = init_mask.mShares[1](j, 0);
            _init_maskB.mShares[0](j - stride, 0) =
                init_mask.mShares[0](j - stride, 0);
            _init_maskB.mShares[1](j - stride, 0) =
                init_mask.mShares[1](j - stride, 0);
        }
        bool_cipher_or(pIdx, _init_maskA, _init_maskB, _init_maskA, enc, eval,
                       runtime);
        for (size_t j = stride; j < len; j++) {
            init_mask.mShares[0](j, 0) = _init_maskA.mShares[0](j - stride, 0);
            init_mask.mShares[1](j, 0) = _init_maskA.mShares[1](j - stride, 0);
        }
    }
    for (size_t i = 0; i < len - 1; i++) {
        res.mShares[0](i, 0) =
            init_mask.mShares[0](i, 0) ^ init_mask.mShares[0](i + 1, 0);
        res.mShares[1](i, 0) =
            init_mask.mShares[1](i, 0) ^ init_mask.mShares[1](i + 1, 0);
    }
    res.mShares[0](len - 1, 0) = init_mask.mShares[0](len - 1, 0) ^ 1;
    res.mShares[1](len - 1, 0) = init_mask.mShares[1](len - 1, 0) ^ 1;
    return;
}

void bool_init_false(int pIdx, aby3::sbMatrix &res) {
    int bitSize = res.bitCount();
    int i64Size = res.i64Size();

    for (int i = 0; i < i64Size; i++) {
        switch (pIdx) {
            case 0:
                res.mShares[0](i, 0) = 1;
                res.mShares[1](i, 0) = 0;
                break;
            case 1:
                res.mShares[0](i, 0) = 1;
                res.mShares[1](i, 0) = 1;
                break;
            case 2:
                res.mShares[0](i, 0) = 0;
                res.mShares[1](i, 0) = 1;
                break;
            default:
                throw std::runtime_error("bool_init_false: pIdx out of range.");
        }
    }
    return;
}

void bool_init_true(int pIdx, aby3::sbMatrix &res) {
    int bitSize = res.bitCount();
    int i64Size = res.i64Size();

    for (int i = 0; i < i64Size; i++) {
        switch (pIdx) {
            case 0:
                res.mShares[0](i, 0) = 0;
                res.mShares[1](i, 0) = 0;
                break;
            case 1:
                res.mShares[0](i, 0) = 1;
                res.mShares[1](i, 0) = 0;
                break;
            case 2:
                res.mShares[0](i, 0) = 0;
                res.mShares[1](i, 0) = 1;
                break;
            default:
                throw std::runtime_error("bool_init_true: pIdx out of range.");
        }
    }
    return;
}

void bool_init_false(int pIdx, boolShare &res) {
    switch (pIdx) {
        case 0:
            res.bshares[0] = 1;
            res.bshares[1] = 0;
            break;
        case 1:
            res.bshares[0] = 1;
            res.bshares[1] = 1;
            break;
        case 2:
            res.bshares[0] = 0;
            res.bshares[1] = 1;
            break;
        default:
            throw std::runtime_error("bool_init_false: pIdx out of range.");
    }
    return;
}

void bool_init_true(int pIdx, boolShare &res) {
    switch (pIdx) {
        case 0:
            res.bshares[0] = 0;
            res.bshares[1] = 0;
            break;
        case 1:
            res.bshares[0] = 1;
            res.bshares[1] = 0;
            break;
        case 2:
            res.bshares[0] = 0;
            res.bshares[1] = 1;
            break;
        default:
            throw std::runtime_error("bool_init_true: pIdx out of range.");
    }
    return;
}

void bool_shift(int pIdx, boolIndex &sharedA, size_t shift_len,
                boolIndex &res_shift, bool right_flag) {
    int bitSize = 64;  // TODO: support more flexible bitsize.

    if (right_flag) {
        res_shift.indexShares[0] = sharedA.indexShares[0] >> shift_len;
        res_shift.indexShares[1] = sharedA.indexShares[1] >> shift_len;
    } else {
        res_shift.indexShares[0] = sharedA.indexShares[0] << shift_len;
        res_shift.indexShares[1] = sharedA.indexShares[1] << shift_len;
    }

    return;
}

void bool_shift_and_left(int pIdx, aby3::sbMatrix &sharedA, size_t shift_len,
                         aby3::sbMatrix &res_shift, aby3::sbMatrix &res_left) {
    // shift the input.
    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    int left_size = bitSize - shift_len;

    // right shift each input to shift_len bits.
    for (size_t i = 0; i < i64Size; i++) {
        res_shift.mShares[0](i, 0) = sharedA.mShares[0](i, 0) >> shift_len;
        res_shift.mShares[1](i, 0) = sharedA.mShares[1](i, 0) >> shift_len;
    }

    // left shift each input to bitSize - shift_len bits.
    for (size_t i = 0; i < i64Size; i++) {
        res_left.mShares[0](i, 0) =
            sharedA.mShares[0](i, 0) & ((1 << shift_len) - 1);
        res_left.mShares[1](i, 0) =
            sharedA.mShares[1](i, 0) & ((1 << shift_len) - 1);
    }

    return;
}

void bool_shift_and_left(int pIdx, boolIndex &sharedA, size_t shift_len,
                         boolIndex &res_shift, boolIndex &res_left) {
    int bitSize = 64;  // TODO: support more flexible bitsize.
    int left_size = bitSize - shift_len;

    bool_shift(pIdx, sharedA, shift_len, res_shift);

    res_left.indexShares[0] = sharedA.indexShares[0] & ((1 << shift_len) - 1);
    res_left.indexShares[1] = sharedA.indexShares[1] & ((1 << shift_len) - 1);

    return;
}


void bool_aggregation(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &res,
                         aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                         aby3::Sh3Runtime &runtime, const std::string &func) {
    // TODO - log round or invocation.
    size_t len = sharedA.rows();
    size_t bitsize = sharedA.bitCount();

    if (!checkPowerOfTwo(len)) {
        // THROW_RUNTIME_ERROR("The size of sharedA must be power of 2!");
        size_t mix_len = roundUpToPowerOfTwo(len);
        size_t left_size = mix_len - len;

        // padding the input to the power of 2.
        aby3::sbMatrix _sharedMixUp(left_size, bitsize);
        _sharedMixUp.mShares[0].setZero();
        _sharedMixUp.mShares[1].setZero();
        sharedA.resize(mix_len, bitsize);

        std::copy(_sharedMixUp.mShares[0].begin(),
                  _sharedMixUp.mShares[0].end(),
                  sharedA.mShares[0].begin() + len);
        std::copy(_sharedMixUp.mShares[1].begin(),
                  _sharedMixUp.mShares[1].end(),
                  sharedA.mShares[1].begin() + len);

        len = mix_len;
    }

    size_t round = (size_t)floor(log2(len));
    size_t mid_len = len;
    aby3::sbMatrix res_last(mid_len, bitsize);
    std::copy(sharedA.mShares.begin(), sharedA.mShares.end(),
              res_last.mShares.begin());

    for (size_t i = 0; i < round; i++) {
        mid_len /= 2;

        aby3::sbMatrix _left_val(mid_len, bitsize);
        aby3::sbMatrix _right_val(mid_len, bitsize);
        aby3::sbMatrix _res_val(mid_len, bitsize);

        std::copy(res_last.mShares[0].begin(),
                  res_last.mShares[0].begin() + mid_len,
                  _left_val.mShares[0].begin());
        std::copy(res_last.mShares[1].begin(),
                  res_last.mShares[1].begin() + mid_len,
                  _left_val.mShares[1].begin());
        std::copy(res_last.mShares[0].begin() + mid_len,
                  res_last.mShares[0].end(), _right_val.mShares[0].begin());
        std::copy(res_last.mShares[1].begin() + mid_len,
                  res_last.mShares[1].end(), _right_val.mShares[1].begin());

        // using or for aggregation.
        if(func == "OR"){
            bool_cipher_or(pIdx, _left_val, _right_val, _res_val, enc, eval,
                       runtime);
        }
        else if(func == "AND"){
            bool_cipher_and(pIdx, _left_val, _right_val, _res_val, enc, eval,
                       runtime);
        }
        else if(func == "ADD"){
            bool_cipher_add(pIdx, _left_val, _right_val, _res_val, enc, eval,
                       runtime);
        }
        else{
            THROW_RUNTIME_ERROR("The function " + func + " is not supported!");
        }


        res_last.resize(mid_len, bitsize);

        std::copy(_res_val.mShares.begin(), _res_val.mShares.end(),
                  res_last.mShares.begin());
    }
    res.resize(res_last.rows(), bitsize);
    std::copy(res_last.mShares.begin(), res_last.mShares.end(),
              res.mShares.begin());
    return;
}


aby3::i64Matrix back2plain(int pIdx, aby3::sbMatrix &cipher_val,
                           aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                           aby3::Sh3Runtime &runtime) {
    // send the shares to the prevParty.
    runtime.mComm.mPrev.asyncSendCopy(cipher_val.mShares[0]);

    // rece from the nextParty.
    aby3::i64Matrix other_share(cipher_val.i64Size(), 1);
    runtime.mComm.mNext.recv(other_share.data(), other_share.size());

    // get the shares.
    for (size_t i = 0; i < cipher_val.i64Size(); i++) {
        other_share(i, 0) = other_share(i, 0) ^ cipher_val.mShares[1](i, 0) ^
                            cipher_val.mShares[0](i, 0);
    }
    return other_share;
}

bool back2plain(int pIdx, boolShare &cipher_val, aby3::Sh3Encryptor &enc,
                aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime) {
    runtime.mComm.mPrev.asyncSendCopy(cipher_val.bshares[0]);
    bool other_share;
    runtime.mComm.mNext.recv(other_share);
    other_share = other_share ^ cipher_val.bshares[1] ^ cipher_val.bshares[0];
    return other_share;
}

aby3::i64 back2plain(int pIdx, boolIndex &cipher_val, aby3::Sh3Encryptor &enc,
                     aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime) {
    runtime.mComm.mPrev.asyncSendCopy(cipher_val.indexShares[0]);
    aby3::i64 other_share = 0;
    runtime.mComm.mNext.recv(other_share);

    other_share =
        other_share ^ cipher_val.indexShares[1] ^ cipher_val.indexShares[0];
    return other_share;
}

std::vector<bool> back2plain(int pIdx, std::vector<boolShare> &cipher_val,
                             aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                             aby3::Sh3Runtime &runtime) {
    std::vector<char> sending_shares(cipher_val.size());
    for (size_t i = 0; i < cipher_val.size(); i++) {
        sending_shares[i] = (char)cipher_val[i].bshares[0];
    }
    runtime.mComm.mPrev.asyncSendCopy(sending_shares.data(),
                                      sending_shares.size());
    runtime.mComm.mNext.recv(sending_shares.data(), sending_shares.size());

    std::vector<bool> other_shares(cipher_val.size());
    for (size_t i = 0; i < cipher_val.size(); i++) {
        other_shares[i] = ((bool)sending_shares[i]) ^ cipher_val[i].bshares[1] ^
                          cipher_val[i].bshares[0];
    }
    return other_shares;
}

void get_permutation(size_t len, std::vector<size_t> &permutation,
                     block &seed) {
    permutation.resize(len);
    for (size_t i = 0; i < len; i++) {
        permutation[i] = i;
    }
    PRNG prng(seed);
    std::random_shuffle(permutation.begin(), permutation.end(), prng);
    return;
}

void get_inverse_permutation(std::vector<size_t> &permutation,
                             std::vector<size_t> &inverse_permutation) {
    size_t len = permutation.size();
    inverse_permutation.resize(len);
    for (size_t i = 0; i < len; i++) {
        inverse_permutation[permutation[i]] = i;
    }
    return;
}

void combine_permutation(std::vector<std::vector<size_t>> &permutation_list,
                         std::vector<size_t> &final_permutation) {
    size_t len = permutation_list[0].size();
    size_t permutes = permutation_list.size();
    final_permutation.resize(len);
    for (size_t i = 0; i < len; i++) {
        final_permutation[i] = i;
    }
    for (size_t i = 0; i < permutes; i++) {
        std::vector<size_t> tmp_inverse = permutation_list[permutes - i - 1];
        get_inverse_permutation(permutation_list[permutes - i - 1],
                                tmp_inverse);
        plain_permutate(tmp_inverse, final_permutation);
    }
    return;
}

void get_random_mask(int pIdx, i64Matrix &res, block &seed) {
    size_t len = res.rows();
    PRNG prng(seed);
    for (size_t i = 0; i < len; i++) res(i, 0) = (i64)prng.get<int64_t>();
    return;
}

void left_shift_and_fill(int pIdx, aby3::sbMatrix& input, int tag_size, int tag_value){

    assert(input.rows() == 1);
    switch(pIdx) {
        case 0:
            input.mShares[0](0, 0) = input.mShares[0](0, 0) << tag_size;
            input.mShares[1](0, 0) = input.mShares[1](0, 0) << tag_size;
            break;
        case 1:
            input.mShares[0](0, 0) = (input.mShares[0](0, 0) << tag_size) | tag_value;
            input.mShares[1](0, 0) = input.mShares[1](0, 0) << tag_size;
            break;
        case 2:
            input.mShares[0](0, 0) = input.mShares[0](0, 0) << tag_size;
            input.mShares[1](0, 0) = (input.mShares[1](0, 0) << tag_size) | tag_value;
            break;
        default:
            THROW_RUNTIME_ERROR("left_shift_and_fill: pIdx out of range.");
    }
    return;
}

void tag_append(int pIdx, std::vector<aby3::sbMatrix>& inputs){

    size_t len = inputs.size();
    size_t bitsize = inputs[0].bitCount(); // assume that the inputs value and the tag size together is less than 64.
    size_t tag_size = std::ceil(std::log2(len));

    // generate the log(len)-size bit tag shares and append to the inputs.
    for(size_t i=0; i<len; i++){
        left_shift_and_fill(pIdx, inputs[i], tag_size, (int) i);
    }

    return;
}

void tag_remove(int pIdx, size_t tag_len, std::vector<aby3::sbMatrix>& inputs){

    size_t total_bit_size = inputs[0].bitCount();
    size_t bitsize = total_bit_size - tag_len;

    // remove the last tag_len bits from the tail of the inputs.
    for(size_t i=0; i<inputs.size(); i++){
        inputs[i].mShares[0](0, 0) = inputs[i].mShares[0](0, 0) >> tag_len;
        inputs[i].mShares[1](0, 0) = inputs[i].mShares[1](0, 0) >> tag_len;
    }

    return;
}

void plain_permutate(std::vector<size_t> &permutation, aby3::sbMatrix &data){
    size_t len = data.rows();
    size_t bitsize = data.bitCount();

    aby3::sbMatrix res(len, bitsize);
    for(size_t i=0; i<len; i++){
        res.mShares[0](i, 0) = data.mShares[0](permutation[i], 0);
        res.mShares[1](i, 0) = data.mShares[1](permutation[i], 0);
    }
    data = res;
    return;
}

void plain_permutate(std::vector<size_t> &permutation, aby3::i64Matrix &data){
    size_t len = data.rows();

    aby3::i64Matrix res(len, 1);
    for(size_t i=0; i<len; i++){
        res(i, 0) = data(permutation[i], 0);
    }
    data = res;
    return;
}