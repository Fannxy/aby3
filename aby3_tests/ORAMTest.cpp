#include "Test.h"
#include <chrono>
#include <random>
#include <thread>
#include <cmath>

#include "../aby3-Basic/Basics.h"
#include "../aby3-Basic/Shuffle.h"
#include "../aby3-Basic/SqrtOram.h"
#include "../aby3-Basic/timer.h"
#include "../aby3-RTR/BuildingBlocks.h"

using namespace oc;
using namespace aby3;
using namespace std;

const int TEST_SIZE = 16;
const int TEST_UNIT_SIZE = 10;
static size_t MAX_COMM_SIZE = 1 << 25;


int pos_map_test(oc::CLP &cmd) {
    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN SQRT-ORAM Position Map TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // generate the test data.
    size_t LINER_TEST_SIZE = 16;

    std::vector<i64Matrix> input_x(LINER_TEST_SIZE);
    for (size_t i = 0; i < LINER_TEST_SIZE; i++) {
        input_x[i].resize(1, 1);
        input_x[i](0, 0) = i;
    }

    // // TEST1 posMap linear.
    size_t pack = 2, S = 32;

    std::vector<sbMatrix> enc_x(LINER_TEST_SIZE);
    for (size_t i = 0; i < LINER_TEST_SIZE; i++) {
        enc_x[i].resize(1, 1);
        if (role == 0) {
            enc.localBinMatrix(runtime, input_x[i], enc_x[i]).get();
        } else {
            enc.remoteBinMatrix(runtime, enc_x[i]).get();
        }
    }
    std::vector<si64> _shared_permutation(LINER_TEST_SIZE);
    efficient_shuffle_with_random_permutation(
        enc_x, role, enc_x, _shared_permutation, enc, eval, runtime);

    std::vector<boolIndex> shared_permutation(LINER_TEST_SIZE);
    for (size_t i = 0; i < LINER_TEST_SIZE; i++)
        shared_permutation[i] = boolIndex(_shared_permutation[i].mData[0],
                                          _shared_permutation[i].mData[1]);

    ABY3PosMap posMap((size_t)LINER_TEST_SIZE, pack, S, shared_permutation,
                      role, enc, eval, runtime);

    // // access each element.
    boolIndex test_index = boolIndex((int)LINER_TEST_SIZE / 3, role);
    boolShare init_fake = boolShare(false, role);
    aby3::i64 phy_index;
    phy_index = posMap.access(test_index, init_fake);

    if (phy_index >= LINER_TEST_SIZE || phy_index < 0) {
        THROW_RUNTIME_ERROR("ERROR: phy_index >= LINER_TEST_SIZE!");
    }
    aby3::sbMatrix test_data;
    test_data = enc_x[phy_index];

    // check the result.
    i64Matrix test_res(1, 1);
    enc.revealAll(runtime, test_data, test_res).get();

    // posMap - recursive branch test.
    size_t RECURSIVE_TEST_SIZE = 64;
    pack = 8;
    S = 16;

    std::vector<i64Matrix> input_x2(RECURSIVE_TEST_SIZE);
    for (size_t i = 0; i < RECURSIVE_TEST_SIZE; i++) {
        input_x2[i].resize(1, 1);
        input_x2[i](0, 0) = i;
    }
    std::vector<sbMatrix> enc_x2(RECURSIVE_TEST_SIZE);
    for (size_t i = 0; i < RECURSIVE_TEST_SIZE; i++) {
        enc_x2[i].resize(1, 1);
        if (role == 0) {
            enc.localBinMatrix(runtime, input_x2[i], enc_x2[i]).get();
        } else {
            enc.remoteBinMatrix(runtime, enc_x2[i]).get();
        }
    }
    std::vector<si64> _shared_permutation2(RECURSIVE_TEST_SIZE);
    efficient_shuffle_with_random_permutation(
        enc_x2, role, enc_x2, _shared_permutation2, enc, eval, runtime);

    std::vector<boolIndex> shared_permutation2(RECURSIVE_TEST_SIZE);
    for (size_t i = 0; i < RECURSIVE_TEST_SIZE; i++)
        shared_permutation2[i] = boolIndex(_shared_permutation2[i].mData[0],
                                           _shared_permutation2[i].mData[1]);

    ABY3PosMap posMap2((size_t)RECURSIVE_TEST_SIZE, pack, S,
                       shared_permutation2, role, enc, eval, runtime);

    boolIndex test_index2 = boolIndex((int)1, role);
    phy_index = posMap2.access(test_index2, init_fake);

    boolIndex test_index3 = boolIndex((int)5, role);
    aby3::i64 phy_index2 = posMap2.access(test_index3, init_fake);

    boolIndex test_index4 = boolIndex((int)1, role);
    aby3::i64 phy_index3 = posMap2.access(test_index4, init_fake);

    if (phy_index >= RECURSIVE_TEST_SIZE || phy_index < 0) {
        THROW_RUNTIME_ERROR("ERROR: phy_index >= RECURSIVE_TEST_SIZE!");
    }
    test_data = enc_x2[phy_index];

    // check the result.
    i64Matrix test_res2(1, 1);
    enc.revealAll(runtime, test_data, test_res2).get();

    i64Matrix test_res3(1, 1);
    enc.revealAll(runtime, enc_x2[phy_index2], test_res3).get();

    i64Matrix test_res4(1, 1);
    enc.revealAll(runtime, enc_x2[phy_index3], test_res4).get();

    if (role == 0) {
        check_result("PosMap - linear case ", test_res(0, 0),
                     LINER_TEST_SIZE / 3);
        check_result("PosMap - recursive case ", test_res2(0, 0), 1);
        check_result("PosMap - recursive case stash miss", test_res3(0, 0), 5);
        check_result("PosMap - recursive case stash hit", test_res4(0, 0), 1);
    }

    return 0;
}

int sqrt_oram_test(oc::CLP &cmd){
    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN SQRT-ORAM TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // generate the test data.
    size_t TEST_SIZE = 1 << 5;

    // currently we do not consider the stashing fulling case.
    size_t stash_size = TEST_SIZE;
    size_t pack_size = 4;
    size_t block_size = 4;

    std::vector<i64Matrix> input_x(TEST_SIZE);
    std::vector<sbMatrix> enc_x(TEST_SIZE);

    for(size_t i=0; i<TEST_SIZE; i++){
        input_x[i].resize(block_size, 1);
        for(size_t j=0; j<block_size; j++){
            input_x[i](j, 0) = i;
        }

        enc_x[i].resize(block_size, 1);
        if(role == 0){
            enc.localBinMatrix(runtime, input_x[i], enc_x[i]).get();
        }else{
            enc.remoteBinMatrix(runtime, enc_x[i]).get();
        }
    }

    // initialize the oram.
    ABY3SqrtOram oram(TEST_SIZE, stash_size, pack_size, role, enc, eval, runtime);
    oram.initiate(enc_x);

    // access the oram.
    std::vector<sbMatrix> accessing_res(TEST_SIZE);

    for(int i=TEST_SIZE-1; i>-1; i--){

        accessing_res[i].resize(block_size, 1);
        boolIndex logical_index = boolIndex(i, role);
        accessing_res[i] = oram.access(logical_index);
    }

    // check the result.
    std::vector<i64Matrix> test_res(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        test_res[i].resize(block_size, 1);
        enc.revealAll(runtime, accessing_res[i], test_res[i]).get();
    }

    if(role == 0){
        check_result("SQRT-ORAM test", test_res[0], input_x[0]);
    }

    return 0;
}

int splitted_oram_init(oc::CLP &cmd){

    SPLITTED_TEST_INIT
    u64 sizeX = cmd.getMany<int>("dataSize")[0];
    u64 block_size = 16;

    // generate the test data.
    std::vector<sbMatrix> input_x(sizeX);
    for(u64 i=0; i<sizeX; i++){
        input_x[i].resize(block_size, 64);
        get_sb_shares(block_size, input_x[i]);
    }

    // initialize the oram.
    u64 stash_size = std::sqrt((double) sizeX);
    u64 pack_size = 4;
    ABY3SqrtOram oram(sizeX, stash_size, pack_size, role, enc, eval, runtime);
    oram.initiate(input_x);

    return 0;
}