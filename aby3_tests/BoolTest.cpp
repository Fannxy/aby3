#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/debug.h"

using namespace oc;
using namespace aby3;
using namespace std;

// #define DEBUG_TEST

const int TEST_SIZE = 16;
const int TEST_UNIT_SIZE = 10;


int bool_basic_test(CLP &cmd) {
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
        debug_info("RUN BOOL TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    // generate test data.
    i64Matrix input_x(TEST_SIZE, 1);
    i64Matrix input_y(TEST_SIZE, 1);
    i64Matrix input_z(TEST_SIZE, 1);
    i64Matrix input_mask(TEST_SIZE, 1);

    i64Matrix res_xor(TEST_SIZE, 1);
    i64Matrix res_gt(TEST_SIZE, 1);
    i64Matrix res_eq(TEST_SIZE, 1);
    i64Matrix res_add(TEST_SIZE, 1);
    i64Matrix res_or(TEST_SIZE, 1);
    i64Matrix res_and(TEST_SIZE, 1);
    i64Matrix res_dot(1, 1);
    i64Matrix res_not(TEST_SIZE, 1);

    for (int i = 0; i < TEST_SIZE; i++) {
        input_x(i, 0) = i;
        input_y(i, 0) = TEST_SIZE - i;
        input_z(i, 0) = -i;
        input_mask(i, 0) =
            (i == 5) ? -1 : -0;  // note the correct mask-based dot computation.
        res_xor(i, 0) = input_x(i, 0) ^ input_y(i, 0);
        res_or(i, 0) = input_x(i, 0) | input_y(i, 0);
        res_and(i, 0) = input_x(i, 0) & input_y(i, 0);
        res_gt(i, 0) = i > (TEST_SIZE - i);
        res_eq(i, 0) = i == (TEST_SIZE - i);
        res_add(i, 0) = TEST_SIZE;
        res_not(i, 0) = ~i;
    }
    res_dot(0, 0) = 5;

    // encrypt the inputs.
    sbMatrix bsharedX(TEST_SIZE, 1);
    sbMatrix bsharedY(TEST_SIZE, 1);
    sbMatrix bsharedZ(TEST_SIZE, 1);
    sbMatrix bsharedMask(TEST_SIZE, 1);
    bsharedMask.resize(TEST_SIZE, 64);
    bsharedX.resize(TEST_SIZE, 64);
    bsharedY.resize(TEST_SIZE, 64);

    if (role == 0) {
        enc.localBinMatrix(runtime, input_x, bsharedX).get();
        enc.localBinMatrix(runtime, input_y, bsharedY).get();
        enc.localBinMatrix(runtime, input_z, bsharedZ).get();
        enc.localBinMatrix(runtime, input_mask, bsharedMask).get();
    } else {
        enc.remoteBinMatrix(runtime, bsharedX).get();
        enc.remoteBinMatrix(runtime, bsharedY).get();
        enc.remoteBinMatrix(runtime, bsharedZ).get();
        enc.remoteBinMatrix(runtime, bsharedMask).get();
    }

    // xor
    sbMatrix shared_xor(TEST_SIZE, 1);
    i64Matrix test_xor(TEST_SIZE, 1);
    for (int i = 0; i < TEST_SIZE; i++) {
        shared_xor.mShares[0](i) =
            bsharedX.mShares[0](i) ^ bsharedY.mShares[0](i);
        shared_xor.mShares[1](i) =
            bsharedX.mShares[1](i) ^ bsharedY.mShares[1](i);
    }
    enc.revealAll(runtime, shared_xor, test_xor).get();

    // or
    sbMatrix shared_or(TEST_SIZE, 1);
    i64Matrix test_or(TEST_SIZE, 1);
    bool_cipher_or(role, bsharedX, bsharedY, shared_or, enc, eval, runtime);
    enc.revealAll(runtime, shared_or, test_or).get();

    // and
    sbMatrix shared_and(TEST_SIZE, 1);
    i64Matrix test_and(TEST_SIZE, 1);
    bool_cipher_and(role, bsharedX, bsharedY, shared_and, enc, eval, runtime);
    enc.revealAll(runtime, shared_and, test_and).get();

    // gt
    sbMatrix shared_gt(TEST_SIZE, 1);
    i64Matrix test_gt(TEST_SIZE, 1);
    bool_cipher_lt(role, bsharedY, bsharedX, shared_gt, enc, eval, runtime);
    enc.revealAll(runtime, shared_gt, test_gt).get();

    // eq
    sbMatrix shared_eq(TEST_SIZE, 1);
    i64Matrix test_eq(TEST_SIZE, 1);
    bool_cipher_eq(role, bsharedY, bsharedX, shared_eq, enc, eval, runtime);
    enc.revealAll(runtime, shared_eq, test_eq).get();

    // add
    sbMatrix shared_add(TEST_SIZE, 1);
    i64Matrix test_add(TEST_SIZE, 1);
    bool_cipher_add(role, bsharedX, bsharedY, shared_add, enc, eval, runtime);
    enc.revealAll(runtime, shared_add, test_add).get();

    // init false.
    sbMatrix shared_false(TEST_SIZE, 1);
    bool_init_false(role, shared_false);
    i64Matrix res_false(TEST_SIZE, 1);
    i64Matrix test_false(TEST_SIZE, 1);
    for (size_t i = 0; i < TEST_SIZE; i++) res_false(i, 0) = 0;
    enc.revealAll(runtime, shared_false, test_false).get();

    // init true.
    sbMatrix shared_true(TEST_SIZE, 1);
    bool_init_true(role, shared_true);
    i64Matrix res_true(TEST_SIZE, 1);
    i64Matrix test_true(TEST_SIZE, 1);
    for (size_t i = 0; i < TEST_SIZE; i++) res_true(i, 0) = 1;
    enc.revealAll(runtime, shared_true, test_true).get();

    // bool shift and left.
    sbMatrix shared_shift(TEST_SIZE, 1);
    sbMatrix shared_left(TEST_SIZE, 1);
    int shift_val = 5;
    bool_shift_and_left(role, bsharedX, shift_val, shared_shift, shared_left);
    i64Matrix res_shift(TEST_SIZE, 1);
    i64Matrix res_left(TEST_SIZE, 1);

    for (size_t i = 0; i < TEST_SIZE; i++) {
        res_shift(i, 0) = input_x(i, 0) >> shift_val;
        res_left(i, 0) = input_x(i, 0) & ((1 << shift_val) - 1);
    }
    i64Matrix test_shift(TEST_SIZE, 1);
    i64Matrix test_left(TEST_SIZE, 1);
    enc.revealAll(runtime, shared_shift, test_shift).get();
    enc.revealAll(runtime, shared_left, test_left).get();

    // check the back2plain.
    i64Matrix test_back2plain(TEST_SIZE, 1);
    test_back2plain = back2plain(role, shared_xor, enc, eval, runtime);

    boolIndex sharedIndex(shared_xor.mShares[0](0), shared_xor.mShares[1](0));
    i64 test_0 = back2plain(role, sharedIndex, enc, eval, runtime);

    // check the bool dot.
    sbMatrix shared_dot(1, 1);
    bool_cipher_dot(role, bsharedX, bsharedMask, shared_dot, enc, eval,
                    runtime);
    i64Matrix test_dot(1, 1);
    enc.revealAll(runtime, shared_dot, test_dot).get();

    // check not.
    sbMatrix shared_not(TEST_SIZE, 1);
    bool_cipher_not(role, bsharedX, shared_not);
    i64Matrix test_not(TEST_SIZE, 1);
    enc.revealAll(runtime, shared_not, test_not).get();

    // mask 0/1 sbMatrix and boolShare translation.
    sbMatrix shared01(TEST_SIZE, 1);
    i64Matrix mask01(TEST_SIZE, 1);
    std::vector<bool> mask01_res(TEST_SIZE);
    std::vector<boolShare> shared01_vec(TEST_SIZE);
    for (int i = 0; i < TEST_SIZE; i++) {
        mask01(i, 0) = (i % 2 == 0) ? 0 : 1;
        mask01_res[i] = (i % 2 == 0) ? 1 : 0;
    }
    if (role == 0) {
        enc.localBinMatrix(runtime, mask01, shared01).get();
    } else {
        enc.remoteBinMatrix(runtime, shared01).get();
    }
    for (int i = 0; i < TEST_SIZE; i++) {
        shared01_vec[i].from_matrix(shared01.mShares[0](i, 0),
                                    shared01.mShares[1](i, 0));
    }

    std::vector<boolShare> shared_not_vec(TEST_SIZE);
    bool_cipher_not(role, shared01_vec, shared_not_vec);
    std::vector<bool> test_mask01 =
        back2plain(role, shared_not_vec, enc, eval, runtime);

    // check the bool cipher and.
    std::vector<boolShare> input_a(4);
    std::vector<boolShare> input_b(4);
    std::vector<bool> res_and_vec = {false, false, false, true};
    std::vector<bool> res_or_vec = {false, true, true, true};
    // std::vector<bool> plain_input_a = {false, false, true, true};

    input_a[0] = boolShare(false, role);
    input_a[1] = boolShare(false, role);
    input_a[2] = boolShare(true, role);
    input_a[3] = boolShare(true, role);
    input_b[0] = boolShare(false, role);
    input_b[1] = boolShare(true, role);
    input_b[2] = boolShare(false, role);
    input_b[3] = boolShare(true, role);

    std::vector<boolShare> res_shared_and_vec(4);
    std::vector<boolShare> res_shared_or_vec(4);

    for (size_t i = 0; i < 4; i++) {
        bool_cipher_and(role, input_a[i], input_b[i], res_shared_and_vec[i],
                        enc, eval, runtime);
        bool_cipher_or(role, input_a[i], input_b[i], res_shared_or_vec[i], enc,
                       eval, runtime);
    }
    std::vector<bool> test_and_vec =
        back2plain(role, res_shared_and_vec, enc, eval, runtime);
    std::vector<bool> test_or_vec =
        back2plain(role, res_shared_or_vec, enc, eval, runtime);

    // check the bool get_first_zero_mask.
    std::vector<bool> random_bool_array(TEST_SIZE);
    std::vector<bool> res_get_first_zero(TEST_SIZE);
    std::vector<boolShare> initial_mask(TEST_SIZE);

    std::mt19937 gen(16);
    std::bernoulli_distribution dist(0.5);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        random_bool_array[i] = dist(gen);
        initial_mask[i] = boolShare(random_bool_array[i], role);
    }
    bool flag = false;
    for (size_t i = 0; i < TEST_SIZE; i++) {
        if (!flag && (random_bool_array[i] == 0)) {
            res_get_first_zero[i] = 1;
            flag = true;
        } else {
            res_get_first_zero[i] = 0;
        }
    }
    aby3::sbMatrix shared_get_first_zero(TEST_SIZE, 1);
    bool_get_first_zero_mask(role, initial_mask, shared_get_first_zero, enc,
                             eval, runtime);
    aby3::i64Matrix test_get_first_zero_i64(TEST_SIZE, 1);
    enc.revealAll(runtime, shared_get_first_zero, test_get_first_zero_i64)
        .get();

    vecBoolShares bool_shared_get_first_zero(TEST_SIZE);
    bool_shared_get_first_zero.from_matrix(shared_get_first_zero);
    std::vector<bool> test_get_first_zero = back2plain(
        role, bool_shared_get_first_zero.bshares, enc, eval, runtime);

    std::vector<bool> test_boolShare_init =
        back2plain(role, initial_mask, enc, eval, runtime);

    // check the result.
    if (role == 0) {
        check_result("xor", test_xor, res_xor);
        check_result("or", test_or, res_or);
        check_result("gt", test_gt, res_gt);
        check_result("eq", test_eq, res_eq);
        check_result("add", test_add, res_add);
        check_result("and", test_and, res_and);
        check_result("init false", test_false, res_false);
        check_result("init true", test_true, res_true);
        check_result("shift", test_shift, res_shift);
        check_result("left", test_left, res_left);
        check_result("back2plain matrix", test_back2plain, res_xor);
        check_result("back2plain index", test_0, res_xor(0, 0));
        check_result("dot", test_dot(0, 0), res_dot(0, 0));
        check_result("boolShare vector not", test_mask01, mask01_res);
        check_result("boolShare init", test_boolShare_init, random_bool_array);
        check_result("bool and", test_and_vec, res_and_vec);
        check_result("bool or", test_or_vec, res_or_vec);
    }

    return 0;
}

int bool_basic_test2(oc::CLP &cmd) {
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
        debug_info("RUN BOOL TEST2");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    aby3::i64Matrix trueVal(TEST_SIZE, 1);
    aby3::i64Matrix falseVal(TEST_SIZE, 1);
    aby3::i64Matrix max_ref(TEST_SIZE, 1);
    aby3::i64Matrix min_ref(TEST_SIZE, 1);

    for (int i = 0; i < TEST_SIZE; i++) {
        trueVal(i, 0) = 5;
        falseVal(i, 0) = 16;
        max_ref(i, 0) = 16;
        min_ref(i, 0) = 5;
    }

    boolShare tflag(true, role);
    boolShare fflag(false, role);

    aby3::sbMatrix shared_true(TEST_SIZE, 1);
    aby3::sbMatrix shared_false(TEST_SIZE, 1);
    shared_true.resize(TEST_SIZE, 64);
    shared_false.resize(TEST_SIZE, 64);

    if (role == 0) {
        enc.localBinMatrix(runtime, trueVal, shared_true).get();
        enc.localBinMatrix(runtime, falseVal, shared_false).get();
    } else {
        enc.remoteBinMatrix(runtime, shared_true).get();
        enc.remoteBinMatrix(runtime, shared_false).get();
    }

    aby3::sbMatrix test_true(TEST_SIZE, 1);
    aby3::sbMatrix test_false(TEST_SIZE, 1);

    bool_cipher_selector(role, tflag, shared_true, shared_false, test_true, enc,
                         eval, runtime);
    bool_cipher_selector(role, fflag, shared_true, shared_false, test_false,
                         enc, eval, runtime);

    aby3::sbMatrix max_res, min_res;
    bool_cipher_max(role, shared_true, shared_false, max_res, enc, eval,
                    runtime);
    bool_cipher_min(role, shared_true, shared_false, min_res, enc, eval,
                    runtime);

    aby3::i64Matrix res_true(TEST_SIZE, 1);
    aby3::i64Matrix res_false(TEST_SIZE, 1);
    enc.revealAll(runtime, test_true, res_true).get();
    enc.revealAll(runtime, test_false, res_false).get();

    aby3::i64Matrix max_test(TEST_SIZE, 1);
    aby3::i64Matrix min_test(TEST_SIZE, 1);
    enc.revealAll(runtime, max_res, max_test).get();
    enc.revealAll(runtime, min_res, min_test).get();

    if (role == 0) {
        check_result("bool selector true", res_true, trueVal);
        check_result("bool selector false", res_false, falseVal);
        check_result("bool max", max_test, max_ref);
        check_result("bool min", min_test, min_ref);
    }

    bool_cipher_max_min_split(role, shared_true, shared_false, max_res, min_res,
                              enc, eval, runtime);
    enc.revealAll(runtime, max_res, max_test).get();
    enc.revealAll(runtime, min_res, min_test).get();
    if(role == 0){
        check_result("bool max split", max_test, max_ref);
        check_result("bool min split", min_test, min_ref);
    }

    return 0;
}

int get_first_zero_test(CLP &cmd) {
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
        debug_info("RUN GET_FIRST_ZERO TEST!");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    // check the bool get_first_zero_mask.
    std::vector<bool> random_bool_array(TEST_SIZE);
    std::vector<bool> res_get_first_zero(TEST_SIZE);
    std::vector<boolShare> initial_mask(TEST_SIZE);

    std::mt19937 gen(16);
    std::bernoulli_distribution dist(0.5);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        random_bool_array[i] = dist(gen);
        initial_mask[i] = boolShare(random_bool_array[i], role);
    }
    bool flag = false;
    for (size_t i = 0; i < TEST_SIZE; i++) {
        if (!flag && (random_bool_array[i] == 0)) {
            res_get_first_zero[i] = 1;
            flag = true;
        } else {
            res_get_first_zero[i] = 0;
        }
    }
    aby3::sbMatrix shared_get_first_zero(TEST_SIZE, 1);
    bool_get_first_zero_mask(role, initial_mask, shared_get_first_zero, enc,
                             eval, runtime);
    aby3::i64Matrix test_get_first_zero_i64(TEST_SIZE, 1);
    enc.revealAll(runtime, shared_get_first_zero, test_get_first_zero_i64)
        .get();
    vecBoolShares bool_shared_get_first_zero(TEST_SIZE);
    bool_shared_get_first_zero.from_matrix(shared_get_first_zero);
    std::vector<bool> test_get_first_zero = back2plain(
        role, bool_shared_get_first_zero.bshares, enc, eval, runtime);

    // all-zero test case.
    for (size_t i = 0; i < TEST_SIZE; i++) {
        initial_mask[i] = boolShare(false, role);
    }
    std::vector<bool> res_get_first_zero_all_zero(TEST_SIZE, false);
    res_get_first_zero_all_zero[0] = true;
    bool_get_first_zero_mask(role, initial_mask, shared_get_first_zero, enc,
                             eval, runtime);
    bool_shared_get_first_zero.from_matrix(shared_get_first_zero);
    std::vector<bool> test_get_first_zero_all_zero = back2plain(
        role, bool_shared_get_first_zero.bshares, enc, eval, runtime);

    if (role == 0) {
        check_result("get_first_zero_mask random test", test_get_first_zero,
                     res_get_first_zero);
        check_result("get_first_zero_mask all-zero test",
                     test_get_first_zero_all_zero, res_get_first_zero_all_zero);
    }

    return 0;
}

int bool_aggregation_test(oc::CLP& cmd){
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
        debug_info("RUN AGGREGATION TEST!");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    // check the bool aggregation or.
    aby3::i64Matrix input_x(TEST_SIZE, 1);
    aby3::i64 res_agg_or = false;
    aby3::i64 res_agg_add = 0;

    std::mt19937 gen(16);
    std::bernoulli_distribution dist(0.5);
    for(size_t i=0; i<TEST_SIZE; i++){
        bool flag = dist(gen);
        input_x(i, 0) = flag;
        res_agg_or = res_agg_or | flag;

        int iflag = (flag)? 1:0;
        res_agg_add += iflag;
    }

    aby3::sbMatrix bsharedX(TEST_SIZE, 1);
    aby3::sbMatrix bisharedX(TEST_SIZE, BITSIZE);
    aby3::sbMatrix res(1, 1);
    aby3::sbMatrix ires(1, BITSIZE);

    if (role == 0) {
        enc.localBinMatrix(runtime, input_x, bsharedX).get();
        enc.localBinMatrix(runtime, input_x, bisharedX).get();
    } else {
        enc.remoteBinMatrix(runtime, bsharedX).get();
        enc.remoteBinMatrix(runtime, bisharedX).get();
    }

    // bool_aggregation_or(role, bsharedX, res, enc, eval, runtime);
    bool_aggregation(role, bsharedX, res, enc, eval, runtime, "OR");
    bool_aggregation(role, bisharedX, ires, enc, eval, runtime, "ADD");

    aby3::i64Matrix test(1, 1);
    enc.revealAll(runtime, res, test).get();
    aby3::i64Matrix itest(1, 1);
    enc.revealAll(runtime, ires, itest).get();

    if(role == 0){
        check_result("aggregation or", test(0, 0), res_agg_or);
        check_result("aggregation add", itest(0, 0), res_agg_add);
    }

    return 0;
}

int share_conversion_test(oc::CLP& cmd){

    BASIC_TEST_INIT

    if(role == 0){
        debug_info("RUN Share conversion test !");
    }

    // check the bool2arith share conversion.
    aby3::i64Matrix input_x(TEST_SIZE, 1);
    aby3::i64Matrix input_y(TEST_SIZE, 1);
    aby3::i64Matrix lt_res(TEST_SIZE, 1);
    for(size_t i=0; i<TEST_SIZE; i++){
        input_x(i, 0) = i;
        input_y(i, 0) = TEST_SIZE - i;
        lt_res(i, 0) = (i < (TEST_SIZE - i))? 1: 0;
    }

    aby3::sbMatrix bsharedX(TEST_SIZE, 64);
    aby3::sbMatrix bsharedY(TEST_SIZE, 64);
    aby3::sbMatrix bsharedRes(TEST_SIZE, 1);

    if(role == 0){
        enc.localBinMatrix(runtime, input_x, bsharedX).get();
        enc.localBinMatrix(runtime, input_y, bsharedY).get();
    }
    else{
        enc.remoteBinMatrix(runtime, bsharedX).get();
        enc.remoteBinMatrix(runtime, bsharedY).get();
    }

    bool_cipher_lt(role, bsharedX, bsharedY, bsharedRes, enc, eval, runtime);

    aby3::si64Matrix test_lt_res(TEST_SIZE, 1);
    aby3::si64Matrix test_b2a(TEST_SIZE, 1);

    bool2arith(role, bsharedRes, test_lt_res, enc, eval, runtime);
    bool2arith(role, bsharedX, test_b2a, enc, eval, runtime);

    aby3::i64Matrix res_lt(TEST_SIZE, 1);
    enc.revealAll(runtime, test_lt_res, res_lt).get();

    aby3::i64Matrix res_b2a(TEST_SIZE, 1);
    enc.revealAll(runtime, test_b2a, res_b2a).get();

    if(role == 0){
        check_result("bool2arith lt", res_lt, lt_res);
        check_result("bool2arith conversion", res_b2a, input_x);
    }

    return 0;
}