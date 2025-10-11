#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Feddb-Core/feddb.h"
#include "../aby3-Feddb-Core/utils.h"
#include "../aby3-Feddb-Core/genperm.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/debug.h"

using namespace oc;
using namespace aby3;
using namespace std;


int oblivious_idx_select_test(CLP &cmd) {
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
        debug_info("RUN OblIdx TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    size_t data_size = 5;
    size_t idx_size = 4;
    i64Matrix data_plain(data_size, 1);
    i64Matrix idx_plain(idx_size, 1);
    i64Matrix res_plain(idx_size, 1);

    data_plain(0,0)=1;
    data_plain(1,0)=4;
    data_plain(2,0)=9;
    data_plain(3,0)=16;
    data_plain(4,0)=25;

    idx_plain(0,0)=2;
    idx_plain(1,0)=1;
    idx_plain(2,0)=3;
    idx_plain(3,0)=2;

    res_plain(0,0)=9;
    res_plain(1,0)=4;
    res_plain(2,0)=16;
    res_plain(3,0)=9;

    // encrypt the inputs.
    si64Matrix data_shared(data_size, 1);
    si64Matrix idx_shared(idx_size, 1);
    si64Matrix res_shared(idx_size, 1);

    if (role == 0) {
        enc.localIntMatrix(runtime, data_plain, data_shared).get();
        enc.localIntMatrix(runtime, idx_plain, idx_shared).get();
    } else {
        enc.remoteIntMatrix(runtime, data_shared).get();
        enc.remoteIntMatrix(runtime, idx_shared).get();
    }

    oblivious_idx_select(role, data_shared, idx_shared, res_shared, enc, eval, runtime);

    i64Matrix res_test(idx_size, 1);
    enc.revealAll(runtime, res_shared, res_test).get();

    if(role == 0){
        check_result("OblIdx Test", res_test, res_plain);
    }
    
    return 0;
}

int genperm_test(CLP &cmd){
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
        debug_info("RUN GenPerm TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    size_t data_size = 9;
    size_t perm_size = 9;

    i64Matrix data_plain(data_size, 1);
    i64Matrix perm_plain(perm_size, 1);

    
    for(size_t i=0;i<5;i++){
        data_plain(i,0)=i;
    }

    data_plain(5,0)=2;
    data_plain(6,0)=1;
    data_plain(7,0)=3;
    data_plain(8,0)=2;

    perm_plain(0,0)=0;
    perm_plain(1,0)=1;
    perm_plain(2,0)=3;
    perm_plain(3,0)=6;
    perm_plain(4,0)=8;
    perm_plain(5,0)=4;
    perm_plain(6,0)=2;
    perm_plain(7,0)=7;
    perm_plain(8,0)=5;

    si64Matrix data_shared(data_size, 1);

    if (role == 0) {
        enc.localIntMatrix(runtime, data_plain, data_shared).get();
    } else {
        enc.remoteIntMatrix(runtime, data_shared).get();
    }

    si64Matrix perm_shared(perm_size, 1);
    genPerm(role, data_shared, perm_shared, enc, eval, runtime);

    i64Matrix perm_test(perm_size, 1);
    enc.revealAll(runtime, perm_shared, perm_test).get();

    //DEBUG
    // if(role == 0){
    //     std::cout << "perm_test: " << std::endl;
    //     for(size_t i=0;i<perm_size;i++){
    //         std::cout << perm_test(i,0) << " ";
    //     }
    //     std::cout << std::endl ;
    // }

    if(role == 0){
        check_result("GenPerm Test", perm_test, perm_plain);
    }

    return 0;

}

int feddb_shuffle_test(CLP &cmd){
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

    size_t TEST_SIZE = 16;
    // generate the test data.
    i64Matrix input_x(TEST_SIZE, 1);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        input_x(i, 0) = i;
    }
 
     // encrypt the inputs.
    sbMatrix bsharedX(TEST_SIZE, 1);
    if (role == 0) {
        enc.localBinMatrix(runtime, input_x, bsharedX).get();
    } else {
        enc.remoteBinMatrix(runtime, bsharedX).get();
    }
     
    si64Matrix sisharedX(TEST_SIZE, 1);
    if (role == 0) {
        enc.localIntMatrix(runtime, input_x, sisharedX).get();
    } else {
        enc.remoteIntMatrix(runtime, sisharedX).get();
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

    //sbmatrix test :pi2*pi0*pi1
    std::vector<size_t> final_permutation_sb;
    std::vector<std::vector<size_t>> permutation_list_sb;
    if (role == 0) {
        permutation_list_sb = {next_permutation, prev_permutation,
                            other_permutation};
    }
    if (role == 1) {
        permutation_list_sb = {prev_permutation, other_permutation,
                            next_permutation};
    }
    if (role == 2) {
        permutation_list_sb = {other_permutation, next_permutation,
                            prev_permutation};
    }
    combine_permutation(permutation_list_sb, final_permutation_sb);
 
    //明文shuffle的结果
    i64Matrix shuffle_sb_res(TEST_SIZE, 1);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        shuffle_sb_res(i, 0) = input_x(i, 0);
    }
    plain_permutate(final_permutation_sb, shuffle_sb_res);
 
    sbMatrix bsharedShuffle(TEST_SIZE, 1);
    shuffle(role, bsharedX, bsharedShuffle, enc, eval, runtime);

    i64Matrix test_res_1(TEST_SIZE, 1);
    enc.revealAll(runtime, bsharedShuffle, test_res_1).get();
    

    // check the sbMatrix shuffle result.
    if (role == 0) {
        // bool check_flag = true;
        // for (size_t i = 0; i < TEST_SIZE; i++) {
        //     if (test_res_1(i, 0) != shuffle_sb_res(i, 0)) {
        //         check_flag = false;
        //     }  
        // }
        // if (check_flag) {
        //     debug_info("\033[32m sbMatrix SHUFFLE CHECK SUCCESS ! \033[0m\n");
        // } else {
        //     debug_info("\033[31m sbMatrix SHUFFLE CHECK ERROR ! \033[0m\n");
        //     debug_info("True result: \n");
        //     debug_output_matrix(shuffle_sb_res);
        //     debug_info("Func result: \n");
        //     debug_output_matrix(test_res_1);
        // }
        check_result("sbMatrix Shuffle Test", test_res_1, shuffle_sb_res);
    }


    //si64Matrix test : pi2*pi1*pi0
    std::vector<size_t> final_permutation_si;
    std::vector<std::vector<size_t>> permutation_list_si;
    if (role == 0) {
        permutation_list_si = {prev_permutation, next_permutation,
                            other_permutation};
    }
    if (role == 1) {
        permutation_list_si = {other_permutation, prev_permutation,
                            next_permutation};
    }
    if (role == 2) {
        permutation_list_si = {next_permutation, other_permutation,
                            prev_permutation};
    }
    combine_permutation(permutation_list_si, final_permutation_si);

    //明文shuffle的结果
    i64Matrix shuffle_si_res(TEST_SIZE, 1);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        shuffle_si_res(i, 0) = input_x(i, 0);
    }
    plain_permutate(final_permutation_si, shuffle_si_res);

    si64Matrix sisharedShuffle(TEST_SIZE, 1);
    shuffle(role, sisharedX, sisharedShuffle, enc, eval, runtime);

    i64Matrix test_res_2(TEST_SIZE, 1);
    enc.revealAll(runtime, sisharedShuffle, test_res_2).get();

     // check the si64Matrix shuffle result.
     if (role == 0) {
        check_result("si64Matrix Shuffle Test", test_res_2, shuffle_si_res);
    }

    //vector<si64Matrix> test
    vector<i64Matrix> input_y(TEST_SIZE);
    size_t UNIT_SIZE = TEST_SIZE;
    for (size_t i = 0; i < TEST_SIZE; i++) {
        input_y[i].resize(UNIT_SIZE, 1);
        for (size_t j = 0; j < UNIT_SIZE; j++) {
            input_y[i](j, 0) = j;
        }
    }

    vector<si64Matrix> sisharedY(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        sisharedY[i].resize(UNIT_SIZE, 1);
        if (role == 0) {
            enc.localIntMatrix(runtime, input_y[i], sisharedY[i]).get();
        } else {
            enc.remoteIntMatrix(runtime, sisharedY[i]).get();
        }
    }

    //vector明文shuffle的结果
    vector<i64Matrix> shuffle_vec_res(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        shuffle_vec_res[i].resize(UNIT_SIZE, 1);
        for (size_t j = 0; j < UNIT_SIZE; j++) {
            shuffle_vec_res[i](j, 0) = input_y[i](j, 0);
        }
    }
    for(size_t i=0; i<TEST_SIZE; i++){
        plain_permutate(final_permutation_si, shuffle_vec_res[i]);
    }

    vector<si64Matrix> vecSharedShuffle(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        vecSharedShuffle[i].resize(UNIT_SIZE, 1);
    }
    shuffle(role, sisharedY, vecSharedShuffle, enc, eval, runtime);

    vector<i64Matrix> test_res_3(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        test_res_3[i].resize(UNIT_SIZE, 1);
        enc.revealAll(runtime, vecSharedShuffle[i], test_res_3[i]).get();
    }

     // check the si64Matrix shuffle result.
     if (role == 0) {
        bool check_flag = true;
        for (size_t i = 0; i < TEST_SIZE; i++) {
            for (size_t j = 0; j < UNIT_SIZE; j++) {
                if (test_res_3[i](j,0) != shuffle_vec_res[i](j, 0)) {
                    check_flag = false; 
                }
            }


        }
        if (check_flag) {
            debug_info("\033[32m Vector<si64Matrix> SHUFFLE CHECK SUCCESS ! \033[0m\n");
        } else {
            debug_info("\033[31m Vector<si64Matrix> SHUFFLE CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(shuffle_vec_res[i]);
            }
            debug_info("Func result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(test_res_3[i]);
            }
        }
    }

    return 0;

}

int persist_test(CLP &cmd){
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
        debug_info("RUN PERSIST TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    size_t TEST_SIZE = 10;
    // generate the test data.
    i64Matrix input_x(TEST_SIZE, 1);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        input_x(i, 0) = i;
    }

    // persist_plain test
    persist_plain(role, "test_plain", input_x);

    // encrypt the inputs.
    sbMatrix bsharedX(TEST_SIZE, 1);
    if (role == 0) {
        enc.localBinMatrix(runtime, input_x, bsharedX).get();
    } else {
        enc.remoteBinMatrix(runtime, bsharedX).get();
    }
    // persist_sbMatrix test
    persist_cipher(role, "test_sbMatrix", bsharedX);

    si64Matrix sisharedX(TEST_SIZE, 1);
    if (role == 0) {
        enc.localIntMatrix(runtime, input_x, sisharedX).get();
    } else {
        enc.remoteIntMatrix(runtime, sisharedX).get();
    }
    // persist_si64Matrix test
    persist_cipher(role, "test_si64Matrix", sisharedX);

    // read the plain text
    i64Matrix read_plain_x(TEST_SIZE, 1);
    read_plain(role, "test_plain", read_plain_x);
    if (role == 0) {
        check_result("Persist&Read Plain Test", read_plain_x, input_x);
    }
    
    // read the sbMatrix
    sbMatrix read_sbMatrix_x(TEST_SIZE, 1);
    read_cipher(role, "test_sbMatrix", read_sbMatrix_x);
    i64Matrix sbMatrix_test(TEST_SIZE, 1);
    enc.revealAll(runtime, read_sbMatrix_x, sbMatrix_test).get();
  
    if (role == 0) {
        check_result("Persist&Read sbMatrix Test", sbMatrix_test, input_x);
    }

    // read the si64Matrix
    si64Matrix read_si64Matrix_x(TEST_SIZE, 1);
    read_cipher(role, "test_si64Matrix", read_si64Matrix_x);
    i64Matrix si64Matrix_test(TEST_SIZE, 1);
    enc.revealAll(runtime, read_si64Matrix_x, si64Matrix_test).get();
    if (role == 0) {
        check_result("Persist&Read si64Matrix Test", si64Matrix_test, input_x);
    }
    
    return 0;

}

int index_agg_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN INDEX_AGG TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    size_t TEST_SIZE = 8;
    std::vector<i64Matrix> input_data(2);
    input_data[0].resize(TEST_SIZE,1);
    input_data[1].resize(TEST_SIZE,1);

    input_data[0](0,0)=1;
    input_data[0](1,0)=3;
    input_data[0](2,0)=5;
    input_data[0](3,0)=4;
    input_data[0](4,0)=3;
    input_data[0](5,0)=2;
    input_data[0](6,0)=2;
    input_data[0](7,0)=3;

    for(int i=0; i< TEST_SIZE; i++){
        input_data[1](i,0)=i;
    }

    std::vector<si64Matrix> dataShared(2);
    dataShared[0].resize(TEST_SIZE,1);
    dataShared[1].resize(TEST_SIZE,1);
    if (role == 0) {
        enc.localIntMatrix(runtime, input_data[0], dataShared[0]).get();
        enc.localIntMatrix(runtime, input_data[1], dataShared[1]).get();
    } else {
        enc.remoteIntMatrix(runtime, dataShared[0]).get();
        enc.remoteIntMatrix(runtime, dataShared[1]).get();
    }

    i64Matrix idx(TEST_SIZE,1);
    idx(0,0)=0;
    idx(1,0)=5;
    idx(2,0)=6;
    idx(3,0)=1;
    idx(4,0)=4;
    idx(5,0)=7;
    idx(6,0)=3;
    idx(7,0)=2;

    i64Matrix equalFlag(TEST_SIZE,1);
    equalFlag(0,0)=0;
    equalFlag(1,0)=0;
    equalFlag(2,0)=1;
    equalFlag(3,0)=0;
    equalFlag(4,0)=1;
    equalFlag(5,0)=1;
    equalFlag(6,0)=0;
    equalFlag(7,0)=0;

    si64Matrix equalFlagShared(TEST_SIZE,1);
    if (role == 0) {
        enc.localIntMatrix(runtime, equalFlag, equalFlagShared).get();
    } else {
        enc.remoteIntMatrix(runtime, equalFlagShared).get();
    }

    std::vector<si64Matrix> test_res(2),test_res_sorted(2);
    index_agg(role, equalFlagShared, idx, dataShared, test_res, enc, eval, runtime);
    int resRows=test_res[0].rows();
    test_res_sorted[0].resize(resRows,1);
    test_res_sorted[1].resize(resRows,1);

    si64Matrix perm(resRows,1);
    genPerm(role, test_res[0], perm, enc, eval, runtime);
    i64Matrix perm_plain(resRows,1);
    enc.revealAll(runtime, perm, perm_plain).get();
    permutate(role, test_res[0], test_res_sorted[0], perm_plain);
    permutate(role, test_res[1], test_res_sorted[1], perm_plain);

    i64Matrix res_key(resRows,1),res_val(resRows,1);
    res_key(0,0)=2;
    res_key(1,0)=3;
    res_key(2,0)=4;
    res_key(3,0)=5;

    res_val(0,0)=11;
    res_val(1,0)=12;
    res_val(2,0)=3;
    res_val(3,0)=2;

    i64Matrix test_res_key(resRows,1),test_res_val(resRows,1);
    enc.revealAll(runtime, test_res_sorted[0], test_res_key).get();
    enc.revealAll(runtime, test_res_sorted[1], test_res_val).get();
    if (role == 0) {
        check_result("Index Agg Test-key",test_res_key, res_key);
        check_result("Index Agg Test-val",test_res_val, res_val);
    }

    return 0;

}