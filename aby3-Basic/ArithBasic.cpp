#include <aby3/Circuit/CircuitLibrary.h>
#include <aby3/sh3/Sh3BinaryEvaluator.h>

#include <bitset>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "debug.h"
#include "BuildingBlocks.h"
#include "Basics.h"

#define DEBUG_BASIC

using namespace oc;
using namespace aby3;

void arith_aggregation(int pIdx, aby3::si64Matrix &sharedA, aby3::si64Matrix &res,
                         aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
                         aby3::Sh3Runtime &runtime, const std::string& func){
    // TODO - log round or invocation.
    size_t len = sharedA.rows();

    if (!checkPowerOfTwo(len)) {
        // THROW_RUNTIME_ERROR("The size of sharedA must be power of 2!");
        size_t mix_len = roundUpToPowerOfTwo(len);
        size_t left_size = mix_len - len;

        // if(pIdx == 0) debug_info("len = " + std::to_string(len) + ", mix_len = " + std::to_string(mix_len) + ", left_size = " + std::to_string(left_size));

        // for Eigen Matrix data structure.
        sharedA.mShares[0].conservativeResize(mix_len, 1);
        sharedA.mShares[1].conservativeResize(mix_len, 1);

        // for the left size, we need to fill the left size with 0.
        sharedA.mShares[0].block(len, 0, left_size, 1).setZero();
        sharedA.mShares[1].block(len, 0, left_size, 1).setZero();

        len = mix_len;
    }

    size_t round = (size_t)floor(log2(len));
    size_t mid_len = len;
    aby3::si64Matrix res_last(mid_len, 1);
    std::copy(sharedA.mShares.begin(), sharedA.mShares.end(),
              res_last.mShares.begin());
    
    for (size_t i = 0; i < round; i++) {
        mid_len /= 2;

        aby3::si64Matrix _left_val(mid_len, 1);
        aby3::si64Matrix _right_val(mid_len, 1);
        aby3::si64Matrix _res_val(mid_len, 1);

        _left_val.mShares[0] = res_last.mShares[0].block(0, 0, mid_len, 1);
        _left_val.mShares[1] = res_last.mShares[1].block(0, 0, mid_len, 1);
        _right_val.mShares[0] = res_last.mShares[0].block(mid_len, 0, mid_len, 1);
        _right_val.mShares[1] = res_last.mShares[1].block(mid_len, 0, mid_len, 1);

        // using or for aggregation.
        if(func == "ADD"){
            for(size_t j=0; j<_res_val.rows(); j++){
                _res_val.mShares[0](j, 0) = _left_val.mShares[0](j, 0) + _right_val.mShares[0](j, 0);
                _res_val.mShares[1](j, 0) = _left_val.mShares[1](j, 0) + _right_val.mShares[1](j, 0);
            }
        }

        res_last.resize(mid_len, 1);

        std::copy(_res_val.mShares.begin(), _res_val.mShares.end(),
                  res_last.mShares.begin());
    }
    res.resize(res_last.rows(), 1);
    std::copy(res_last.mShares.begin(), res_last.mShares.end(),
              res.mShares.begin());
    return;
}