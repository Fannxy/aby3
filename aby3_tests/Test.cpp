#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

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

bool check_result(const std::string &func_name, i64Matrix &test,
                  i64Matrix &res) {
    int size = test.rows();
    bool check_flag = true;
    for (int i = 0; i < size; i++) {
        if (test(i, 0) != res(i, 0)) check_flag = false;
    }

    if (!check_flag) {
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
#ifdef SHOW_TEST_CASE
        debug_info("test case: ");
        debug_output_matrix(test);
        debug_info("expected result: ");
        debug_output_matrix(res);
#endif
    } else {
        debug_info("\033[32m" + func_name + " SUCCESS!" + "\033[0m\n");
    }

    return check_flag;
}

bool check_result(const std::string &func_name, aby3::i64 test, aby3::i64 res) {
    bool check_flag = true;
    if (test != res) check_flag = false;

    if (!check_flag) {
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
        debug_info("test res = " + to_string(test) + "\n");
        debug_info("expected res = " + to_string(res) + "\n");
    } else {
        debug_info("\033[32m" + func_name + " SUCCESS!" + "\033[0m\n");
    }

    return check_flag;
}

bool check_result(const std::string& func_name, std::vector<aby3::i64Matrix> test, std::vector<aby3::i64Matrix> res){
    bool check_flag = true;
    for(size_t i=0; i<test.size(); i++){
        for(size_t j=0; j<test[i].rows(); j++){
            if(test[i](j, 0) != res[i](j, 0)) check_flag = false;
        }
    }
    if(!check_flag){
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
    }else{
        debug_info("\033[32m" + func_name + " SUCCESS!" + "\033[0m\n");
    
    }
}

int arith_basic_test(CLP &cmd) {
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
        debug_info("RUN ARITH TEST");
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
    f64Matrix<D8> finput_x(TEST_SIZE, 1);
    i64Matrix res_gt(TEST_SIZE, 1);
    i64Matrix res_ge(TEST_SIZE, 1);
    i64Matrix res_eq(TEST_SIZE, 1);
    i64Matrix res_mul(TEST_SIZE, 1);
    i64Matrix res_ab_mul(TEST_SIZE, 1);
    i64Matrix res_ib_mul(TEST_SIZE, 1);
    f64Matrix<D8> res_f_mul(TEST_SIZE, 1);
    f64Matrix<D8> res_fb_mul(TEST_SIZE, 1);

    for (int i = 0; i < TEST_SIZE; i++) {
        input_x(i, 0) = i;
        finput_x(i, 0) = i;
        input_y(i, 0) = TEST_SIZE - i;

        res_gt(i, 0) = i > (TEST_SIZE - i);
        res_ge(i, 0) = i >= (TEST_SIZE - i);
        res_mul(i, 0) = i * (TEST_SIZE - i);
        res_eq(i, 0) = i == (TEST_SIZE - i);
        res_ab_mul(i, 0) = res_gt(i, 0) * input_x(i, 0);
        res_ib_mul(i, 0) = i * res_gt(i, 0);
        res_f_mul(i, 0) = i * i;
        res_fb_mul(i, 0) = i * res_gt(i, 0);
    }

    // encrypt the inputs.
    si64Matrix sharedX(TEST_SIZE, 1);
    si64Matrix sharedY(TEST_SIZE, 1);
    sf64Matrix<D8> fsharedX(TEST_SIZE, 1);
    if (role == 0) {
        enc.localIntMatrix(runtime, input_x, sharedX).get();
        enc.localIntMatrix(runtime, input_y, sharedY).get();
        enc.localFixedMatrix(runtime, finput_x, fsharedX).get();
    } else {
        enc.remoteIntMatrix(runtime, sharedX).get();
        enc.remoteIntMatrix(runtime, sharedY).get();
        enc.remoteFixedMatrix(runtime, fsharedX).get();
    }

    // compute for the result.
    sbMatrix shared_gt(TEST_SIZE, 1);
    sbMatrix shared_ge(TEST_SIZE, 1);
    sbMatrix shared_eq(TEST_SIZE, 1);

    si64Matrix shared_mul(TEST_SIZE, 1);
    si64Matrix shared_ab_mul(TEST_SIZE, 1);
    si64Matrix shared_ib_mul(TEST_SIZE, 1);
    sf64Matrix<D8> shared_f_mul(TEST_SIZE, 1);
    sf64Matrix<D8> shared_fb_mul(TEST_SIZE, 1);

    cipher_gt(role, sharedX, sharedY, shared_gt, eval, runtime);
    circuit_cipher_eq(role, sharedX, sharedY, shared_eq, eval, runtime);
    cipher_ge(role, sharedX, sharedY, shared_ge, eval, enc, runtime);

    cipher_mul(role, sharedX, sharedY, shared_mul, eval, enc, runtime);
    cipher_mul(role, sharedX, shared_gt, shared_ab_mul, eval, enc, runtime);
    cipher_mul(role, input_x, shared_gt, shared_ib_mul, eval, enc, runtime);
    // std::cout << "mul3" << std::endl;
    cipher_mul<D8>(role, fsharedX, fsharedX, shared_f_mul, eval, enc, runtime);
    cipher_mul<D8>(role, fsharedX, shared_gt, shared_fb_mul, eval, enc,
                   runtime);

    // check the result.
    i64Matrix test_gt(TEST_SIZE, 1);
    i64Matrix test_ge(TEST_SIZE, 1);
    i64Matrix test_eq(TEST_SIZE, 1);
    i64Matrix test_mul(TEST_SIZE, 1);
    i64Matrix test_ab_mul(TEST_SIZE, 1);
    i64Matrix test_ib_mul(TEST_SIZE, 1);
    f64Matrix<D8> test_f_mul(TEST_SIZE, 1);
    f64Matrix<D8> test_fb_mul(TEST_SIZE, 1);

    enc.revealAll(runtime, shared_gt, test_gt).get();
    enc.revealAll(runtime, shared_ge, test_ge).get();
    enc.revealAll(runtime, shared_eq, test_eq).get();
    enc.revealAll(runtime, shared_mul, test_mul).get();
    enc.revealAll(runtime, shared_ab_mul, test_ab_mul).get();
    enc.revealAll(runtime, shared_ib_mul, test_ib_mul).get();
    enc.revealAll(runtime, shared_f_mul, test_f_mul).get();
    enc.revealAll(runtime, shared_fb_mul, test_fb_mul).get();

    if (role == 0) {
        check_result("gt", test_gt, res_gt);
        check_result("ge", test_ge, res_ge);
        check_result("eq", test_eq, res_eq);
        check_result("mul", test_mul, res_mul);
        check_result("mul_ab", test_ab_mul, res_ab_mul);
        check_result("mul_ib", test_ib_mul, res_ib_mul);
        check_result("mul_f", test_f_mul, res_f_mul);
        check_result("mul_fb", test_fb_mul, res_fb_mul);
    }
    return 0;
}

int arith_basic_test2(CLP& cmd){
    // get the configs.
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
        debug_info("RUN ARITH TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    // generate test data.
    aby3::si64Matrix test_x(TEST_SIZE, 1);
    std::vector<int> test_y(TEST_SIZE);
    aby3::i64Matrix plain_x(TEST_SIZE, 1);
    aby3::i64Matrix res(TEST_SIZE, 1);

    for(size_t i=0; i<TEST_SIZE; i++){
        plain_x(i, 0) = 1;
        test_y[i] = i;
        if(i == 1) res(i, 0) = 1;
        else res(i, 0) = 0;
    }

    if(role == 0){
        enc.localIntMatrix(runtime, plain_x, test_x).get();
    }else{
        enc.remoteIntMatrix(runtime, test_x).get();
    }

    // compute for the result.
    std::vector<si64> test_x_vec(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        test_x_vec[i].mData[0] = test_x.mShares[0](i, 0);
        test_x_vec[i].mData[1] = test_x.mShares[1](i, 0);
    }

    aby3::sbMatrix test_res(TEST_SIZE, 1);
    vector_cipher_eq(role, test_x_vec, test_y, test_res, eval, runtime);

    // check the result.
    aby3::i64Matrix test_res_plain(TEST_SIZE, 1);
    enc.revealAll(runtime, test_res, test_res_plain).get();

    if(role == 0){
        check_result("vector_cipher_eq", test_res_plain, res);
    }


    return 0;
}

int initialization_test(CLP &cmd) {
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
        debug_info("RUN INIT TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // 1. test the correlated randomness.
    oc::block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    oc::block prevSeed = enc.mShareGen.mPrevCommon.getSeed();

    // send the nextSeed to the nextParty.
    runtime.mComm.mNext.asyncSendCopy(nextSeed);

    // send the prevSeed to the prevParty.
    runtime.mComm.mPrev.asyncSendCopy(prevSeed);

    // get the nextParty's prev seed.
    oc::block nextP_seed;
    runtime.mComm.mNext.recv(nextP_seed);

    // get the prevParty's next seed.
    oc::block prevP_seed;
    runtime.mComm.mPrev.recv(prevP_seed);

    // if(role == 0){
    //     debug_info("next seed: ");
    //     debug_info(nextSeed);
    //     debug_info(nextP_seed);
    //     debug_info("prev seed: ");
    //     debug_info(prevSeed);
    //     debug_info(prevP_seed);
    // }

    // check the seeds.
    if (nextSeed != nextP_seed) {
        debug_info("\033[31m P" + to_string(role) + " check: nextSeed ERROR!" +
                   "\033[0m\n");
    } else {
        debug_info("\033[32m P" + to_string(role) +
                   " check: nextSeed SUCCESS!" + "\033[0m\n");
    }

    if (prevSeed != prevP_seed) {
        debug_info("\033[31m P" + to_string(role) + " check: prevSeed ERROR!" +
                   "\033[0m\n");
    } else {
        debug_info("\033[32m P" + to_string(role) +
                   " check: prevSeed SUCCESS!" + "\033[0m\n");
    }

    return 0;
}

int splitted_shuffle_test(oc::CLP &cmd){
    SPLITTED_TEST_INIT

    if(role == 0){
        debug_info("RUN Shuffle TEST");
    }

    // prepare the data.
    size_t data_size;
    if(cmd.isSet("dataSize")){
        data_size = cmd.getMany<int>("dataSize")[0];
    } else {
        throw std::runtime_error(LOCATION);
    }
    size_t bin_size;
    if(cmd.isSet("binSize")){
        bin_size = cmd.getMany<int>("binSize")[0];
    } else {
        throw std::runtime_error(LOCATION);
    }
    aby3::i64Matrix data_plain(data_size, 1);
    aby3::i64Matrix data_res(data_size, 1);
    std::vector<int> data_res_vec(data_size);

    for(size_t i=0; i<data_size; i++){
        data_plain(i, 0) = (data_size - i) % bin_size;
        data_res_vec[i] = (i+1) % bin_size;
    }

    std::sort(data_res_vec.begin(), data_res_vec.end());
    for(size_t i=0; i<data_size; i++){
        data_res(i, 0) = data_res_vec[i];
    }

    // todo - add the deplication tags.
    aby3::sbMatrix enc_data(data_size, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, data_plain, enc_data).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_data).get();
    }

    std::vector<aby3::sbMatrix> enc_data_vec(data_size);
    for(size_t i=0; i<data_size; i++){
        aby3::sbMatrix tmp(1, 64);
        tmp.mShares[0](0, 0) = enc_data.mShares[0](i, 0);
        tmp.mShares[1](0, 0) = enc_data.mShares[1](i, 0);
        enc_data_vec[i] = tmp;
    }

    efficient_shuffle(enc_data_vec, role, enc_data_vec, enc, eval, runtime);

    aby3::sbMatrix enc_test(data_size, 64);
    for(size_t i=0; i<data_size; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    aby3::i64Matrix test(data_size, 1);
    enc.revealAll(runtime, enc_test, test).get();
    return 0;
}

int shuffle_test(oc::CLP &cmd) {
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
        debug_info("RUN SHUFFLE TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // generate the test data.
    std::vector<i64Matrix> input_x(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        input_x[i].resize(TEST_UNIT_SIZE, 1);
        for (size_t j = 0; j < TEST_UNIT_SIZE; j++) {
            input_x[i](j, 0) = i;
        }
    }

    // encrypt the inputs.
    std::vector<sbMatrix> bsharedX(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        bsharedX[i].resize(TEST_UNIT_SIZE, 1);
        if (role == 0) {
            enc.localBinMatrix(runtime, input_x[i], bsharedX[i]).get();
        } else {
            enc.remoteBinMatrix(runtime, bsharedX[i]).get();
        }
    }

    // generate the permutation.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = TEST_SIZE;
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    std::vector<size_t> other_permutation(len);
    runtime.mComm.mPrev.asyncSendCopy(next_permutation.data(),
                                      next_permutation.size());
    runtime.mComm.mNext.recv(other_permutation.data(),
                             other_permutation.size());

    std::vector<size_t> final_permutation;
    std::vector<std::vector<size_t>> permutation_list;
    if (role == 0) {
        permutation_list = {next_permutation, prev_permutation,
                            other_permutation};
    }
    if (role == 1) {
        permutation_list = {prev_permutation, other_permutation,
                            next_permutation};
    }
    if (role == 2) {
        permutation_list = {other_permutation, next_permutation,
                            prev_permutation};
    }
    combine_permutation(permutation_list, final_permutation);

    std::vector<i64Matrix> shuffle_res(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        shuffle_res[i].resize(TEST_UNIT_SIZE, 1);
        shuffle_res[i] = input_x[i];
    }
    plain_permutate(final_permutation, shuffle_res);

    std::vector<sbMatrix> bsharedShuffle(TEST_SIZE);

    efficient_shuffle(bsharedX, role, bsharedShuffle, enc, eval, runtime);

    std::vector<i64Matrix> test_res(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        test_res[i].resize(TEST_UNIT_SIZE, 1);
        enc.revealAll(runtime, bsharedShuffle[i], test_res[i]).get();
    }

    // check the shuffle result.
    if (role == 0) {
        bool check_flag = true;
        for (size_t i = 0; i < TEST_SIZE; i++) {
            for (size_t j = 0; j < TEST_UNIT_SIZE; j++) {
                if (test_res[i](j, 0) != shuffle_res[i](j, 0)) {
                    check_flag = false;
                }
            }
        }
        if (check_flag) {
            debug_info("\033[32m SHUFFLE CHECK SUCCESS ! \033[0m\n");
        } else {
            debug_info("\033[31m SHUFFLE CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(shuffle_res[i]);
            }
            debug_info("Func result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(test_res[i]);
            }
        }
    }

    // shuffle with permutation test.
    std::vector<si64> shared_permutation(TEST_SIZE);
    efficient_shuffle_with_random_permutation(
        bsharedX, role, bsharedShuffle, shared_permutation, enc, eval, runtime);

    std::vector<i64Matrix> test_res2(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        test_res2[i].resize(TEST_UNIT_SIZE, 1);
        enc.revealAll(runtime, bsharedShuffle[i], test_res2[i]).get();
    }
    sbMatrix shared_permutation_matrix(TEST_SIZE, 1);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        shared_permutation_matrix.mShares[0](i) =
            shared_permutation[i].mData[0];
        shared_permutation_matrix.mShares[1](i) =
            shared_permutation[i].mData[1];
    }
    i64Matrix test_permutation(TEST_SIZE, 1);
    enc.revealAll(runtime, shared_permutation_matrix, test_permutation).get();

    bool check_shuffle = true;
    bool check_permutation = true;

    for (size_t i = 0; i < TEST_SIZE; i++) {
        if (test_permutation(i, 0) != final_permutation[i]) {
            check_permutation = false;
        }
        for (size_t j = 0; j < TEST_UNIT_SIZE; j++) {
            if (test_res2[i](j, 0) != shuffle_res[i](j, 0)) {
                check_shuffle = false;
            }
        }
    }

    // check the final result.
    if (role == 0) {
        if (check_shuffle) {
            debug_info(
                "\033[32m SHUFFLE in shuffle and permutation CHECK SUCCESS ! "
                "\033[0m\n");
        } else {
            debug_info(
                "\033[31m SHUFFLE in in shuffle and permutation CHECK ERROR ! "
                "\033[0m\n");
            debug_info("True result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(shuffle_res[i]);
            }
            debug_info("Func result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(test_res2[i]);
            }
        }

        if (check_permutation) {
            debug_info(
                "\033[32m PERMUTATION in shuffle and permutation CHECK SUCCESS "
                "! "
                "\033[0m\n");
        } else {
            debug_info(
                "\033[31m PERMUTATION in in shuffle and permutation CHECK "
                "ERROR ! "
                "\033[0m\n");
            debug_info("True result: \n");
            debug_output_vector(final_permutation);

            debug_info("Func result: \n");
            debug_output_matrix(test_permutation);
        }
    }
    return 0;
}

int large_scale_shuffle_test(oc::CLP &cmd){
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
        debug_info("RUN LARGE-SCALE SHUFFLE TEST \nNote this function will cause 150GB Memory+, BE CAREFUL!");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    size_t LARGE_SCALE_SIZE = 1 << 26;
    size_t LARGE_UNIT_SIZE = 1;

    // generate the test data.
    // std::vector<i64Matrix> input_x(LARGE_SCALE_SIZE);
    // for (size_t i = 0; i < LARGE_SCALE_SIZE; i++) {
    //     input_x[i].resize(LARGE_UNIT_SIZE, 1);
    //     for (size_t j = 0; j < LARGE_UNIT_SIZE; j++) {
    //         input_x[i](j, 0) = i;
    //     }
    // }
    i64Matrix input_x(LARGE_SCALE_SIZE*LARGE_UNIT_SIZE, 1);
    for (size_t i = 0; i < LARGE_SCALE_SIZE; i++) {
        for (size_t j = 0; j < LARGE_UNIT_SIZE; j++) {
            input_x(i*LARGE_UNIT_SIZE + j, 0) = i;
        }
    }
    sbMatrix bsharedXMat(LARGE_SCALE_SIZE*LARGE_UNIT_SIZE, 1);
    size_t round = (size_t)ceil(LARGE_SCALE_SIZE*LARGE_UNIT_SIZE / (double)MAX_COMM_SIZE);
    size_t last_len = LARGE_SCALE_SIZE*LARGE_UNIT_SIZE - (round - 1) * MAX_COMM_SIZE;
    for(size_t i=0; i<round; i++){
        size_t len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
        i64Matrix _data = input_x.block(i*MAX_COMM_SIZE, 0, len, 1);
        sbMatrix _bsharedX(len, 1);
        if(role == 0){
            enc.localBinMatrix(runtime, _data, _bsharedX).get();
        }else{
            enc.remoteBinMatrix(runtime, _bsharedX).get();
        }
        std::memcpy(bsharedXMat.mShares[0].data() + i*MAX_COMM_SIZE, _bsharedX.mShares[0].data(), len*sizeof(_bsharedX.mShares[0](0, 0)));
        std::memcpy(bsharedXMat.mShares[1].data() + i*MAX_COMM_SIZE, _bsharedX.mShares[1].data(), len*sizeof(_bsharedX.mShares[0](0, 0)));
    }


    // encrypt the inputs.
    std::vector<sbMatrix> bsharedX(LARGE_SCALE_SIZE);
    for (size_t i = 0; i < LARGE_SCALE_SIZE; i++) {
        bsharedX[i].resize(LARGE_UNIT_SIZE, 1);
        for(size_t j=0; j<LARGE_UNIT_SIZE; j++){
            bsharedX[i].mShares[0](j, 0) = bsharedXMat.mShares[0](i*LARGE_UNIT_SIZE + j, 0);
            bsharedX[i].mShares[1](j, 0) = bsharedXMat.mShares[1](i*LARGE_UNIT_SIZE + j, 0);
        }
    }

    // if(role == 0) debug_info("after test data construction!");

    // generate the permutation.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = LARGE_SCALE_SIZE;
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    std::vector<size_t> other_permutation(len);

    // runtime.mComm.mPrev.asyncSendCopy(next_permutation.data(),
    //                                   next_permutation.size());
    // runtime.mComm.mNext.recv(other_permutation.data(),
    //                          other_permutation.size());

    for(size_t i=0; i<round; i++){
        size_t len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
        auto sendFu = runtime.mComm.mPrev.asyncSendFuture(next_permutation.data() + i*MAX_COMM_SIZE, len);
        auto recvFu = runtime.mComm.mNext.asyncRecv(other_permutation.data() + i*MAX_COMM_SIZE, len);
        sendFu.get();
        recvFu.get();
    }

    // if(role == 0) debug_info("after get permutation!");

    std::vector<size_t> final_permutation;
    std::vector<std::vector<size_t>> permutation_list;
    if (role == 0) {
        permutation_list = {next_permutation, prev_permutation,
                            other_permutation};
    }
    if (role == 1) {
        permutation_list = {prev_permutation, other_permutation,
                            next_permutation};
    }
    if (role == 2) {
        permutation_list = {other_permutation, next_permutation,
                            prev_permutation};
    }
    combine_permutation(permutation_list, final_permutation);

    std::vector<i64Matrix> shuffle_res(LARGE_SCALE_SIZE);
    for (size_t i = 0; i < LARGE_SCALE_SIZE; i++) {
        shuffle_res[i].resize(LARGE_UNIT_SIZE, 1);
        // shuffle_res[i] = input_x[i];
        for(size_t j=0; j<LARGE_UNIT_SIZE; j++){
            shuffle_res[i](j, 0) = input_x(i*LARGE_UNIT_SIZE + j, 0);
        }
    }
    plain_permutate(final_permutation, shuffle_res);

    // shuffle with permutation test.
    std::vector<aby3::sbMatrix> bsharedShuffle(LARGE_SCALE_SIZE);
    std::vector<si64> shared_permutation(LARGE_SCALE_SIZE);

    // if(role == 0) debug_info("outter, before in permutation!!!");
    efficient_shuffle_with_random_permutation(
        bsharedX, role, bsharedShuffle, shared_permutation, enc, eval, runtime);
    // if(role == 0) debug_info("outter, after in permutation!!!");

    // std::vector<i64Matrix> test_res2(LARGE_SCALE_SIZE);
    // for (size_t i = 0; i < LARGE_SCALE_SIZE; i++) {
    //     test_res2[i].resize(LARGE_UNIT_SIZE, 1);
    //     enc.revealAll(runtime, bsharedShuffle[i], test_res2[i]).get();
    // }

    i64Matrix test_res2_full(LARGE_SCALE_SIZE*LARGE_UNIT_SIZE, 1);
    sbMatrix expand_bshared_full(LARGE_SCALE_SIZE*LARGE_UNIT_SIZE, 1);
    for(size_t i=0; i<LARGE_SCALE_SIZE; i++){
        for(size_t j=0; j<LARGE_UNIT_SIZE; j++){
            expand_bshared_full.mShares[0](i*LARGE_UNIT_SIZE + j, 0) = bsharedShuffle[i].mShares[0](j, 0);
            expand_bshared_full.mShares[1](i*LARGE_UNIT_SIZE + j, 0) = bsharedShuffle[i].mShares[1](j, 0);
        }
    }

    for(size_t i=0; i<round; i++){
        size_t len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
        sbMatrix target_bshared(len, 1);
        std::memcpy(target_bshared.mShares[0].data(), expand_bshared_full.mShares[0].data() + i*MAX_COMM_SIZE, len*sizeof(expand_bshared_full.mShares[0](0, 0)));
        std::memcpy(target_bshared.mShares[1].data(), expand_bshared_full.mShares[1].data() + i*MAX_COMM_SIZE, len*sizeof(expand_bshared_full.mShares[0](0, 0)));

        i64Matrix target_test_res(len, 1);
        enc.revealAll(runtime, target_bshared, target_test_res).get();
        test_res2_full.block(i*MAX_COMM_SIZE, 0, len, 1) = target_test_res;
    }

    // if(role == 0) debug_info("data reveal done!");

    std::vector<i64Matrix> test_res2(LARGE_SCALE_SIZE);
    for (size_t i = 0; i < LARGE_SCALE_SIZE; i++) {
        test_res2[i].resize(LARGE_UNIT_SIZE, 1);
        for(size_t j=0; j<LARGE_UNIT_SIZE; j++){
            test_res2[i](j, 0) = test_res2_full(i*LARGE_UNIT_SIZE + j, 0);
        }
    }

    sbMatrix shared_permutation_matrix(LARGE_SCALE_SIZE, 1);
    for (size_t i = 0; i < LARGE_SCALE_SIZE; i++) {
        shared_permutation_matrix.mShares[0](i) =
            shared_permutation[i].mData[0];
        shared_permutation_matrix.mShares[1](i) =
            shared_permutation[i].mData[1];
    }
    i64Matrix test_permutation(LARGE_SCALE_SIZE, 1);

    round = (size_t)ceil(LARGE_SCALE_SIZE / (double)MAX_COMM_SIZE);
    last_len = LARGE_SCALE_SIZE - (round - 1) * MAX_COMM_SIZE;

    for(size_t i=0; i<round; i++){
        size_t len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
        sbMatrix target_shared_permutation(len, 1);
        std::memcpy(target_shared_permutation.mShares[0].data(), shared_permutation_matrix.mShares[0].data() + i*MAX_COMM_SIZE, len*sizeof(shared_permutation_matrix.mShares[0](0, 0)));
        std::memcpy(target_shared_permutation.mShares[1].data(), shared_permutation_matrix.mShares[1].data() + i*MAX_COMM_SIZE, len*sizeof(shared_permutation_matrix.mShares[0](0, 0)));

        i64Matrix target_test_permutation(len, 1);
        enc.revealAll(runtime, target_shared_permutation, target_test_permutation).get();
        test_permutation.block(i*MAX_COMM_SIZE, 0, len, 1) = target_test_permutation;
    }

    bool check_shuffle = true;
    bool check_permutation = true;

    for (size_t i = 0; i < LARGE_SCALE_SIZE; i++) {
        if (test_permutation(i, 0) != final_permutation[i]) {
            check_permutation = false;
        }
        for (size_t j = 0; j < LARGE_UNIT_SIZE; j++) {
            if (test_res2[i](j, 0) != shuffle_res[i](j, 0)) {
                check_shuffle = false;
            }
        }
    }

    // check the final result.
    if (role == 0) {
        if (check_shuffle) {
            debug_info(
                "\033[32m SHUFFLE in shuffle and permutation CHECK SUCCESS ! "
                "\033[0m\n");
        } else {
            debug_info(
                "\033[31m SHUFFLE in in shuffle and permutation CHECK ERROR ! "
                "\033[0m\n");
            debug_info("True result: \n");
            int print_length = (LARGE_SCALE_SIZE > 25) ? 25 : LARGE_SCALE_SIZE;
            for (size_t i = 0; i < print_length; i++) {
                debug_output_matrix(shuffle_res[i]);
            }
            debug_info("Func result: \n");
            for (size_t i = 0; i < print_length; i++) {
                debug_output_matrix(test_res2[i]);
            }
        }

        if (check_permutation) {
            debug_info(
                "\033[32m PERMUTATION in shuffle and permutation CHECK SUCCESS "
                "! "
                "\033[0m\n");
        } else {
            debug_info(
                "\033[31m PERMUTATION in in shuffle and permutation CHECK "
                "ERROR ! "
                "\033[0m\n");
            debug_info("True result: \n");
            debug_output_vector(final_permutation);

            debug_info("Func result: \n");
            debug_output_matrix(test_permutation);
        }
    }
    return 0;
}

int correlation_test(oc::CLP &cmd) {
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
        debug_info("RUN CORRELATION TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    i64Matrix prevMask(TEST_SIZE, 1);
    i64Matrix nextMask(TEST_SIZE, 1);

    // generate the permutation.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();

    get_random_mask(role, prevMask, prevSeed);
    get_random_mask(role, nextMask, nextSeed);

    // check.
    i64Matrix nextP_mask(TEST_SIZE, 1);
    i64Matrix prevP_mask(TEST_SIZE, 1);

    runtime.mComm.mNext.asyncSendCopy(nextMask.data(), nextMask.size());
    runtime.mComm.mPrev.asyncSendCopy(prevMask.data(), prevMask.size());

    runtime.mComm.mNext.recv(nextP_mask.data(), nextP_mask.size());
    runtime.mComm.mPrev.recv(prevP_mask.data(), prevP_mask.size());

    bool next_check = true;
    bool prev_check = true;

    for (size_t i = 0; i < TEST_SIZE; i++) {
        if (nextMask(i, 0) != nextP_mask(i, 0)) {
            debug_info(to_string(nextMask(i, 0)) + " " +
                       to_string(nextP_mask(i, 0)) + "\n");
            next_check = false;
        }
        if (prevMask(i, 0) != prevP_mask(i, 0)) {
            debug_info(to_string(prevMask(i, 0)) + " " +
                       to_string(prevP_mask(i, 0)) + "\n");
            prev_check = false;
        }
    }

    if (!next_check) {
        debug_info("\033[31m P" + to_string(role) +
                   " check: next randomness ERROR!" + "\033[0m\n");
    } else {
        debug_info("\033[32m P" + to_string(role) +
                   " check: next randomness SUCCESS!" + "\033[0m\n");
    }

    if (!prev_check) {
        debug_info("\033[31m P" + to_string(role) +
                   " check: prev randomness ERROR!" + "\033[0m\n");
    } else {
        debug_info("\033[32m P" + to_string(role) +
                   " check: prev randomness SUCCESS!" + "\033[0m\n");
    }

    return 0;
}

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

int communication_test(oc::CLP &cmd){
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
        debug_info("RUN Communication TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // test the communication.
    size_t test_size = 1 << 30;
    aby3::i64Matrix test_a(test_size, 1);
    for(size_t i=0; i<test_size; i++){
        test_a(i, 0) = role;
    }
    aby3::i64Matrix buffer(test_size, 1);


    // block-wise communicate
    size_t block_size = 1 << 28;
    size_t round = (size_t)ceil(test_size / (double)block_size);
    size_t last_len = test_size - (round - 1) * block_size;

    for(size_t i=0; i<round; i++){
        size_t len = (i == round - 1) ? last_len : block_size;
        aby3::i64Matrix sub_test_a = test_a.block(i * block_size, 0, len, 1);
        aby3::i64Matrix sub_buffer(len, 1);
        auto send_fu = runtime.mComm.mNext.asyncSendFuture(sub_test_a.data(), sub_test_a.size());
        auto recv_fu = runtime.mComm.mPrev.asyncRecv(sub_buffer.data(), sub_buffer.size());
        send_fu.get();
        recv_fu.get();
        buffer.block(i * block_size, 0, len, 1) = sub_buffer;
    }

    aby3::i64Matrix checking_res(test_size, 1);
    for(size_t i=0; i<test_size; i++){
        if(role == 0) checking_res(i, 0) = buffer(i, 0) - 2;
        else checking_res(i, 0) = buffer(i, 0) + 1;
    }

    check_result("Ring Communication test ROLE " + std::to_string(role), checking_res, test_a);

    // only 0 -> 1 and 1 -> 2.
    if(role == 0){
        for(size_t i=0; i<round; i++){
            size_t len = (i == round - 1) ? last_len : block_size;
            aby3::i64Matrix sub_test_a = test_a.block(i * block_size, 0, len, 1);
            aby3::i64Matrix sub_buffer(len, 1);
            auto sendFu = runtime.mComm.mNext.asyncSendFuture(sub_test_a.data(), sub_test_a.size());
            sendFu.get();
        }
    }
    if(role == 1){
        aby3::i64Matrix new_buf(test_size, 1);
        // receiving from 0.
        for(size_t i=0; i<round; i++){
            size_t len = (i == round - 1) ? last_len : block_size;

            // recv from p0.
            aby3::i64Matrix sub_buffer(len, 1);
            auto recvFu = runtime.mComm.mPrev.asyncRecv(sub_buffer.data(), sub_buffer.size());

            // send to p2.
            aby3::i64Matrix send_buf = test_a.block(i * block_size, 0, len, 1);
            auto sendFu = runtime.mComm.mNext.asyncSendFuture(send_buf.data(), send_buf.size());

            recvFu.get();
            new_buf.block(i * block_size, 0, len, 1) = sub_buffer;
            sendFu.get();
        }
        check_result("Single Communication test ROLE " + std::to_string(role), new_buf, buffer);
    }
    if(role == 2){
        // receiving from 1.
        aby3::i64Matrix new_buf(test_size, 1);
        for(size_t i=0; i<round; i++){
            size_t len = (i == round - 1) ? last_len : block_size;
            aby3::i64Matrix sub_buffer(len, 1);
            runtime.mComm.mPrev.recv(sub_buffer.data(), sub_buffer.size());
            new_buf.block(i * block_size, 0, len, 1) = sub_buffer;
        }
        check_result("Single Communication test ROLE " + std::to_string(role), new_buf, buffer);
    }

    // large-scale communication test, between p1 and p0.
    if(role == 1){
        large_data_sending(role, test_a, runtime, false);
    }
    if(role == 0){
        aby3::i64Matrix from1_test(test_size, 1);
        large_data_receiving(role, from1_test, runtime, false);
        aby3::i64Matrix res(test_size, 1);
        std::fill_n(res.data(), test_size, 1);
        check_result("Large-scale Communication test ROLE " + std::to_string(role), from1_test, res);
    }
    double testd = 0;
    synchronized_time(role, testd, runtime);

    return 0;
}