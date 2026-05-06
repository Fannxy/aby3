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

#define EMPTY_VALUE std::numeric_limits<i64>::max()


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

    // ========== vector<si64Matrix> version test ==========
    if (role == 0) {
        debug_info("RUN OblIdx Vector TEST");
    }

    // 构造3列数据: data[0]={1,4,9,16,25}, data[1]={10,20,30,40,50}, data[2]={100,200,300,400,500}
    size_t v_data_size = 5;
    size_t v_idx_size = 4;
    size_t v_num_cols = 3;

    std::vector<i64Matrix> v_data_plain(v_num_cols);
    v_data_plain[0].resize(v_data_size, 1);
    v_data_plain[0](0,0)=1;  v_data_plain[0](1,0)=4;  v_data_plain[0](2,0)=9;  v_data_plain[0](3,0)=16; v_data_plain[0](4,0)=25;
    v_data_plain[1].resize(v_data_size, 1);
    v_data_plain[1](0,0)=10; v_data_plain[1](1,0)=20; v_data_plain[1](2,0)=30; v_data_plain[1](3,0)=40; v_data_plain[1](4,0)=50;
    v_data_plain[2].resize(v_data_size, 1);
    v_data_plain[2](0,0)=100;v_data_plain[2](1,0)=200;v_data_plain[2](2,0)=300;v_data_plain[2](3,0)=400;v_data_plain[2](4,0)=500;

    // idx={2,1,3,2}, 期望结果: col0={9,4,16,9}, col1={30,20,40,30}, col2={300,200,400,300}
    i64Matrix v_idx_plain(v_idx_size, 1);
    v_idx_plain(0,0)=2; v_idx_plain(1,0)=1; v_idx_plain(2,0)=3; v_idx_plain(3,0)=2;

    std::vector<i64Matrix> v_res_expected(v_num_cols);
    v_res_expected[0].resize(v_idx_size, 1);
    v_res_expected[0](0,0)=9;  v_res_expected[0](1,0)=4;  v_res_expected[0](2,0)=16; v_res_expected[0](3,0)=9;
    v_res_expected[1].resize(v_idx_size, 1);
    v_res_expected[1](0,0)=30; v_res_expected[1](1,0)=20; v_res_expected[1](2,0)=40; v_res_expected[1](3,0)=30;
    v_res_expected[2].resize(v_idx_size, 1);
    v_res_expected[2](0,0)=300;v_res_expected[2](1,0)=200;v_res_expected[2](2,0)=400;v_res_expected[2](3,0)=300;

    // encrypt
    std::vector<si64Matrix> v_data_shared(v_num_cols);
    si64Matrix v_idx_shared(v_idx_size, 1);
    for(size_t c = 0; c < v_num_cols; c++){
        v_data_shared[c].resize(v_data_size, 1);
        if (role == 0) {
            enc.localIntMatrix(runtime, v_data_plain[c], v_data_shared[c]).get();
        } else {
            enc.remoteIntMatrix(runtime, v_data_shared[c]).get();
        }
    }
    if (role == 0) {
        enc.localIntMatrix(runtime, v_idx_plain, v_idx_shared).get();
    } else {
        enc.remoteIntMatrix(runtime, v_idx_shared).get();
    }

    // call vector version
    std::vector<si64Matrix> v_res_shared;
    oblivious_idx_select(role, v_data_shared, v_idx_shared, v_res_shared, enc, eval, runtime);

    // verify
    if(role == 0){
        for(size_t c = 0; c < v_num_cols; c++){
            i64Matrix v_res_test(v_idx_size, 1);
            enc.revealAll(runtime, v_res_shared[c], v_res_test).get();
            check_result("OblIdx Vector Test col" + std::to_string(c), v_res_test, v_res_expected[c]);
        }
    } else {
        for(size_t c = 0; c < v_num_cols; c++){
            i64Matrix v_res_test(v_idx_size, 1);
            enc.revealAll(runtime, v_res_shared[c], v_res_test).get();
        }
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


    //genPerm_kv test
    i64Matrix key_plain(data_size, 1);
    i64Matrix value_plain(data_size, 1);
    i64Matrix perm_kv_plain(perm_size, 1);
    for(size_t i=0;i<5;i++){
        key_plain(i,0)=i;
    }

    key_plain(5,0)=2;
    key_plain(6,0)=1;
    key_plain(7,0)=3;
    key_plain(8,0)=2;

    value_plain(0,0)=0;
    value_plain(1,0)=6;
    value_plain(2,0)=5;
    value_plain(3,0)=7;
    value_plain(4,0)=4;
    value_plain(5,0)=2;
    value_plain(6,0)=1;
    value_plain(7,0)=3;
    value_plain(8,0)=8;

    perm_kv_plain(0,0)=0;
    perm_kv_plain(1,0)=2;
    perm_kv_plain(2,0)=4;
    perm_kv_plain(3,0)=7;
    perm_kv_plain(4,0)=8;
    perm_kv_plain(5,0)=3;
    perm_kv_plain(6,0)=1;
    perm_kv_plain(7,0)=6;
    perm_kv_plain(8,0)=5;

    si64Matrix key_shared(data_size, 1);
    si64Matrix value_shared(data_size, 1);
    si64Matrix k_v_shared(data_size, 2);
    si64Matrix perm_kv_shared(perm_size, 1);
    if (role == 0) {
        enc.localIntMatrix(runtime, key_plain, key_shared).get();
        enc.localIntMatrix(runtime, value_plain, value_shared).get();
    } else {
        enc.remoteIntMatrix(runtime, key_shared).get();
        enc.remoteIntMatrix(runtime, value_shared).get();
    }
    for(size_t i=0; i<data_size; i++){
        k_v_shared.mShares[0](i, 1) = key_shared.mShares[0](i, 0);
        k_v_shared.mShares[1](i, 1) = key_shared.mShares[1](i, 0);
        k_v_shared.mShares[0](i, 0) = value_shared.mShares[0](i, 0);
        k_v_shared.mShares[1](i, 0) = value_shared.mShares[1](i, 0);
    }
    //genPerm_kv(role, key_shared, value_shared, perm_kv_shared, enc, eval, runtime);
    genPerm(role, k_v_shared, perm_kv_shared, enc, eval, runtime);
    i64Matrix perm_kv_test(perm_size, 1);
    enc.revealAll(runtime, perm_kv_shared, perm_kv_test).get();
    if(role == 0){
        check_result("GenPerm_kv Test", perm_kv_test, perm_kv_plain);
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
    size_t UNIT_SIZE = TEST_SIZE;
    // generate the test data.
    std::vector<i64Matrix> input_x(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        input_x[i].resize(UNIT_SIZE, 1);
        for (size_t j = 0; j < UNIT_SIZE; j++) {
            input_x[i](j, 0) = j;
        }
    }

    // persist_plain test
    persist_plain(role, "test_plain", input_x);

    // encrypt the inputs.
    std::vector<sbMatrix> bsharedX(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        bsharedX[i].resize(UNIT_SIZE, 1);
        if (role == 0) {
            enc.localBinMatrix(runtime, input_x[i], bsharedX[i]).get();
        } else {
            enc.remoteBinMatrix(runtime, bsharedX[i]).get();
        }
    }
    // persist_sbMatrix test
    persist_cipher(role, "test_sbMatrix", bsharedX);

    std::vector<si64Matrix> sisharedX(TEST_SIZE);
    for (size_t i = 0; i < TEST_SIZE; i++) {
        sisharedX[i].resize(UNIT_SIZE, 1);
        if (role == 0) {
            enc.localIntMatrix(runtime, input_x[i], sisharedX[i]).get();
        } else {
            enc.remoteIntMatrix(runtime, sisharedX[i]).get();
        }
    }
    // persist_si64Matrix test
    persist_cipher(role, "test_si64Matrix", sisharedX);

    // read the plain text
    std::vector<i64Matrix> read_plain_x(TEST_SIZE);
    read_plain(role, "test_plain", read_plain_x);

    if (role == 0) {
        bool check_flag = true;
        for (size_t i = 0; i < TEST_SIZE; i++) {
            for (size_t j = 0; j < UNIT_SIZE; j++) {
                if (read_plain_x[i](j,0) != input_x[i](j, 0)) {
                    check_flag = false; 
                }
            }


        }
        if (check_flag) {
            debug_info("\033[32m Vector<i64Matrix> PRESIST&READ CHECK SUCCESS ! \033[0m\n");
        } else {
            debug_info("\033[31m Vector<i64Matrix> PRESIST&READ CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(input_x[i]);
            }
            debug_info("Func result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(read_plain_x[i]);
            }
        }
    }
    
    // read the sbMatrix
    std::vector<sbMatrix> read_sbMatrix_x(TEST_SIZE);
    read_cipher(role, "test_sbMatrix", read_sbMatrix_x);
    std::vector<i64Matrix> sbMatrix_test(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        sbMatrix_test[i].resize(UNIT_SIZE, 1);
        enc.revealAll(runtime, read_sbMatrix_x[i], sbMatrix_test[i]).get();
    }

    if (role == 0) {
        bool check_flag = true;
        for (size_t i = 0; i < TEST_SIZE; i++) {
            for (size_t j = 0; j < UNIT_SIZE; j++) {
                if (sbMatrix_test[i](j,0) != input_x[i](j, 0)) {
                    check_flag = false; 
                }
            }


        }
        if (check_flag) {
            debug_info("\033[32m Vector<sb64Matrix> PRESIST&READ CHECK SUCCESS ! \033[0m\n");
        } else {
            debug_info("\033[31m Vector<sb64Matrix> PRESIST&READ CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(input_x[i]);
            }
            debug_info("Func result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(sbMatrix_test[i]);
            }
        }
    }

    // read the si64Matrix
    std::vector<si64Matrix> read_si64Matrix_x(TEST_SIZE);
    read_cipher(role, "test_si64Matrix", read_si64Matrix_x);
    std::vector<i64Matrix> si64Matrix_test(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        si64Matrix_test[i].resize(UNIT_SIZE, 1);
        enc.revealAll(runtime, read_si64Matrix_x[i], si64Matrix_test[i]).get();
    }

    if (role == 0) {
        bool check_flag = true;
        for (size_t i = 0; i < TEST_SIZE; i++) {
            for (size_t j = 0; j < UNIT_SIZE; j++) {
                if (si64Matrix_test[i](j,0) != input_x[i](j, 0)) {
                    check_flag = false; 
                }
            }


        }
        if (check_flag) {
            debug_info("\033[32m Vector<si64Matrix> PRESIST&READ CHECK SUCCESS ! \033[0m\n");
        } else {
            debug_info("\033[31m Vector<si64Matrix> PRESIST&READ CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(input_x[i]);
            }
            debug_info("Func result: \n");
            for (size_t i = 0; i < TEST_SIZE; i++) {
                debug_output_matrix(si64Matrix_test[i]);
            }
        }
    }
    
    return 0;

}

//TODO key_col expand test
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

    std::vector<si64Matrix> dataShared_key(1);
    si64Matrix dataShared_val(TEST_SIZE,1);
    dataShared_key[0].resize(TEST_SIZE,1);

    if (role == 0) {
        enc.localIntMatrix(runtime, input_data[0], dataShared_key[0]).get();
        enc.localIntMatrix(runtime, input_data[1], dataShared_val).get();
    } else {
        enc.remoteIntMatrix(runtime, dataShared_key[0]).get();
        enc.remoteIntMatrix(runtime, dataShared_val).get();
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
    index_agg(role, equalFlagShared, idx, dataShared_key, dataShared_val, test_res, enc, eval, runtime);
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

int index_agg_maxmin_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN INDEX_AGG_MAXMIN TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // Same fixture style as index_agg_test.
    size_t TEST_SIZE = 8;
    std::vector<i64Matrix> input_data(2);
    input_data[0].resize(TEST_SIZE,1); // key
    input_data[1].resize(TEST_SIZE,1); // val

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

    std::vector<si64Matrix> dataShared_key(1);
    si64Matrix dataShared_val(TEST_SIZE,1);
    dataShared_key[0].resize(TEST_SIZE,1);

    if (role == 0) {
        enc.localIntMatrix(runtime, input_data[0], dataShared_key[0]).get();
        enc.localIntMatrix(runtime, input_data[1], dataShared_val).get();
    } else {
        enc.remoteIntMatrix(runtime, dataShared_key[0]).get();
        enc.remoteIntMatrix(runtime, dataShared_val).get();
    }

    // Public sorted order by key and segment flags.
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

    // Run MAX and MIN.
    std::vector<si64Matrix> max_res(2), min_res(2);
    index_agg_maxmin(role, equalFlagShared, idx, dataShared_key, dataShared_val, true,  max_res, enc, eval, runtime);
    index_agg_maxmin(role, equalFlagShared, idx, dataShared_key, dataShared_val, false, min_res, enc, eval, runtime);

    auto sort_by_key = [&](std::vector<si64Matrix>& in_res, std::vector<si64Matrix>& out_res){
        int resRows = in_res[0].rows();
        out_res.resize(2);
        out_res[0].resize(resRows,1);
        out_res[1].resize(resRows,1);
        si64Matrix perm(resRows,1);
        genPerm(role, in_res[0], perm, enc, eval, runtime);
        i64Matrix perm_plain(resRows,1);
        enc.revealAll(runtime, perm, perm_plain).get();
        permutate(role, in_res[0], out_res[0], perm_plain);
        permutate(role, in_res[1], out_res[1], perm_plain);
    };

    std::vector<si64Matrix> max_sorted, min_sorted;
    sort_by_key(max_res, max_sorted);
    sort_by_key(min_res, min_sorted);

    int outRows = max_sorted[0].rows();
    i64Matrix max_key(outRows,1), max_val(outRows,1);
    i64Matrix min_key(outRows,1), min_val(outRows,1);
    enc.revealAll(runtime, max_sorted[0], max_key).get();
    enc.revealAll(runtime, max_sorted[1], max_val).get();
    enc.revealAll(runtime, min_sorted[0], min_key).get();
    enc.revealAll(runtime, min_sorted[1], min_val).get();

    // Expected:
    // sorted keys: [1,2,3,4,5]
    // vals by group: key1->{0}, key2->{5,6}, key3->{1,4,7}, key4->{3}, key5->{2}
    // max: [0,6,7,3,2], min: [0,5,1,3,2]
    i64Matrix expected_key(5,1), expected_max(5,1), expected_min(5,1);
    expected_key(0,0)=1; expected_key(1,0)=2; expected_key(2,0)=3; expected_key(3,0)=4; expected_key(4,0)=5;
    expected_max(0,0)=0; expected_max(1,0)=6; expected_max(2,0)=7; expected_max(3,0)=3; expected_max(4,0)=2;
    expected_min(0,0)=0; expected_min(1,0)=5; expected_min(2,0)=1; expected_min(3,0)=3; expected_min(4,0)=2;

    if (role == 0) {
        check_result("Index Agg MAX Test-key", max_key, expected_key);
        check_result("Index Agg MAX Test-val", max_val, expected_max);
        check_result("Index Agg MIN Test-key", min_key, expected_key);
        check_result("Index Agg MIN Test-val", min_val, expected_min);
    }

    return 0;
}

int group_by_common_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN GROUP_BY_COMMON TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    
    size_t TEST_SIZE = 5;
    i64Matrix key(TEST_SIZE,2);
    i64Matrix val(TEST_SIZE,1);
    key(0,0)=1;
    key(1,0)=3;
    key(2,0)=1;
    key(3,0)=3;
    key(4,0)=2;
    key(0,1)=1;
    key(1,1)=2;
    key(2,1)=1;
    key(3,1)=3;
    key(4,1)=1;

    val(0,0)=2;
    val(1,0)=4;
    val(2,0)=3;
    val(3,0)=5;
    val(4,0)=1;

    i64Matrix key_G(TEST_SIZE,2);
    key_G(0,0)=1;
    key_G(1,0)=1;
    key_G(2,0)=2;
    key_G(3,0)=3;
    key_G(4,0)=3;
    key_G(0,1)=1;
    key_G(1,1)=1;
    key_G(2,1)=1;
    key_G(3,1)=2;
    key_G(4,1)=3;

    i64Matrix val_G(TEST_SIZE,1);
    val_G(0,0)=2;
    val_G(1,0)=3;
    val_G(2,0)=1;
    val_G(3,0)=4;
    val_G(4,0)=5;

    i64Matrix e(TEST_SIZE,1);
    e(0,0)=0;
    e(1,0)=1;
    e(2,0)=1;
    e(3,0)=1;
    e(4,0)=1;

    i64Matrix key_GN(TEST_SIZE,2);
    key_GN(0,0)=EMPTY_VALUE;
    key_GN(1,0)=1;
    key_GN(2,0)=2;
    key_GN(3,0)=3;
    key_GN(4,0)=3;
    key_GN(0,1)=EMPTY_VALUE;
    key_GN(1,1)=1;
    key_GN(2,1)=1;
    key_GN(3,1)=2;
    key_GN(4,1)=3; 

    i64Matrix perm_GN(TEST_SIZE,1);
    perm_GN(0,0)=4;
    perm_GN(1,0)=0;
    perm_GN(2,0)=1;
    perm_GN(3,0)=2;
    perm_GN(4,0)=3;

    i64Matrix key_out(TEST_SIZE,2);
    key_out(0,0)=1;
    key_out(1,0)=2;
    key_out(2,0)=3;
    key_out(3,0)=3;
    key_out(4,0)=EMPTY_VALUE;
    key_out(0,1)=1;
    key_out(1,1)=1;
    key_out(2,1)=2;
    key_out(3,1)=3;
    key_out(4,1)=EMPTY_VALUE;


    si64Matrix keyShared(TEST_SIZE,2);
    si64Matrix valShared(TEST_SIZE,1);
    if (role == 0) {
        enc.localIntMatrix(runtime, key, keyShared).get();
        enc.localIntMatrix(runtime, val, valShared).get();
    } else {
        enc.remoteIntMatrix(runtime, keyShared).get();
        enc.remoteIntMatrix(runtime, valShared).get();
    }

    si64Matrix keyGShared(TEST_SIZE,2);
    si64Matrix valGShared(TEST_SIZE,1);
    si64Matrix eShared(TEST_SIZE,1);
    si64Matrix perm_GNShared(TEST_SIZE,1);
    si64Matrix key_GNShared(TEST_SIZE,2);
    si64Matrix key_outShared(TEST_SIZE,2);
    group_by_common(role, keyShared, valShared, keyGShared, valGShared, eShared, perm_GNShared, key_GNShared, key_outShared, enc, eval, runtime);

    i64Matrix keyGPlain(TEST_SIZE,2);
    enc.revealAll(runtime, keyGShared, keyGPlain).get();
    i64Matrix valGPlain(TEST_SIZE,1);
    enc.revealAll(runtime, valGShared, valGPlain).get();
    i64Matrix ePlain(TEST_SIZE,1);
    enc.revealAll(runtime, eShared, ePlain).get();
    i64Matrix perm_GNPlain(TEST_SIZE,1);
    enc.revealAll(runtime, perm_GNShared, perm_GNPlain).get();
    i64Matrix key_GNPlain(TEST_SIZE,2);
    enc.revealAll(runtime, key_GNShared, key_GNPlain).get();
    i64Matrix key_outPlain(TEST_SIZE,2);
    enc.revealAll(runtime, key_outShared, key_outPlain).get();

    if (role == 0) {
        check_result("Group By Common Test-key",keyGPlain, key_G);
        check_result("Group By Common Test-val",valGPlain, val_G);
        check_result("Group By Common Test-e",ePlain, e);
        check_result("Group By Common Test-perm_GN",perm_GNPlain, perm_GN);
        check_result("Group By Common Test-key_GN",key_GNPlain, key_GN);
        check_result("Group By Common Test-key_out",key_outPlain, key_out);
    }

    return 0;

}

int group_by_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN GROUP_BY TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    
    size_t TEST_SIZE = 5;
    std::vector<i64Matrix> key(2);
    key[0].resize(TEST_SIZE,1);
    key[1].resize(TEST_SIZE,1);
    i64Matrix val(TEST_SIZE,1);
    key[0](0,0)=1;
    key[0](1,0)=3;
    key[0](2,0)=1;
    key[0](3,0)=3;
    key[0](4,0)=2;
    key[1](0,0)=1;
    key[1](1,0)=2;
    key[1](2,0)=1;
    key[1](3,0)=3;
    key[1](4,0)=1;


    val(0,0)=2;
    val(1,0)=4;
    val(2,0)=3;
    val(3,0)=5;
    val(4,0)=1;
    
    std::vector<i64Matrix> key_out(2);
    key_out[0].resize(TEST_SIZE,1);
    key_out[1].resize(TEST_SIZE,1);
    key_out[0](0,0)=1;
    key_out[0](1,0)=2;
    key_out[0](2,0)=3;
    key_out[0](3,0)=3;
    key_out[0](4,0)=EMPTY_VALUE;
    key_out[1](0,0)=1;
    key_out[1](1,0)=1;
    key_out[1](2,0)=2;
    key_out[1](3,0)=3;
    key_out[1](4,0)=EMPTY_VALUE;


    i64Matrix count(TEST_SIZE, 1);
    count(0,0)=2;
    count(1,0)=1;
    count(2,0)=1;
    count(3,0)=1;
    count(4,0)=0;

    i64Matrix sum(TEST_SIZE, 1);
    sum(0,0)=5;
    sum(1,0)=1;
    sum(2,0)=4;
    sum(3,0)=5;
    sum(4,0)=0;

    i64Matrix max(TEST_SIZE, 1);
    max(0,0)=3;
    max(1,0)=1;
    max(2,0)=4;
    max(3,0)=5;
    max(4,0)=0;

    i64Matrix min(TEST_SIZE, 1);
    min(0,0)=2;
    min(1,0)=1;
    min(2,0)=4;
    min(3,0)=5;
    min(4,0)=0;

    std::vector<si64Matrix> keyShared(2);
    keyShared[0].resize(TEST_SIZE,1);
    keyShared[1].resize(TEST_SIZE,1);
    si64Matrix valShared(TEST_SIZE,1);
    if (role == 0) {
        enc.localIntMatrix(runtime, key[0], keyShared[0]).get();
        enc.localIntMatrix(runtime, key[1], keyShared[1]).get();
        enc.localIntMatrix(runtime, val, valShared).get();
    } else {
        enc.remoteIntMatrix(runtime, keyShared[0]).get();
        enc.remoteIntMatrix(runtime, keyShared[1]).get();
        enc.remoteIntMatrix(runtime, valShared).get();
    }

    //group_count
    si64Matrix countShared(TEST_SIZE,1);
    std::vector<si64Matrix> key_outShared(2);
    key_outShared[0].resize(TEST_SIZE,1);
    key_outShared[1].resize(TEST_SIZE,1);
    group_count(role, keyShared, key_outShared, countShared, enc, eval, runtime);

    i64Matrix countPlain(TEST_SIZE,1);
    enc.revealAll(runtime, countShared, countPlain).get();
    std::vector<i64Matrix> key_outPlain(2);
    key_outPlain[0].resize(TEST_SIZE,1);
    key_outPlain[1].resize(TEST_SIZE,1);
    enc.revealAll(runtime, key_outShared[0], key_outPlain[0]).get();
    enc.revealAll(runtime, key_outShared[1], key_outPlain[1]).get();

    if (role == 0) {
        check_result("Group By Count Test-key-0",key_outPlain[0], key_out[0]);
        check_result("Group By Count Test-key-1",key_outPlain[1], key_out[1]);
        check_result("Group By Count Test-count",countPlain, count);

    }

    //group_sum
    si64Matrix sumShared(TEST_SIZE,1);
    group_sum(role, keyShared, valShared, key_outShared, sumShared, enc, eval, runtime);

    i64Matrix sumPlain(TEST_SIZE,1);
    enc.revealAll(runtime, sumShared, sumPlain).get();
    if (role == 0) {
        check_result("Group By Sum Test-sum",sumPlain, sum);
    }

    //group_max
    si64Matrix maxShared(TEST_SIZE,1);
    group_max(role, keyShared, valShared, key_outShared, maxShared, enc, eval, runtime);

    i64Matrix maxPlain(TEST_SIZE,1);
    enc.revealAll(runtime, maxShared, maxPlain).get();
    if (role == 0) {
        check_result("Group By Max Test-max",maxPlain, max);
    }

    //group_min
    si64Matrix minShared(TEST_SIZE,1);
    group_min(role, keyShared, valShared, key_outShared, minShared, enc, eval, runtime);

    i64Matrix minPlain(TEST_SIZE,1);
    enc.revealAll(runtime, minShared, minPlain).get();
    if (role == 0) {
        check_result("Group By Min Test-min",minPlain, min);
    }

    return 0;

}

int odd_even_merge_sort_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN ODD_EVEN_MERGE_SORT TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    
    size_t TEST_SIZE = 5;
    i64Matrix data(TEST_SIZE,1);
    data(0,0)=15;
    data(1,0)=14;
    data(2,0)=13;
    data(3,0)=10;
    data(4,0)=4;

    i64Matrix res(TEST_SIZE,1);
    res(0,0)=4;
    res(1,0)=10;
    res(2,0)=13;
    res(3,0)=14;
    res(4,0)=15;

    si64Matrix dataSi64Shared(TEST_SIZE,1);
    if (role == 0) {
        enc.localIntMatrix(runtime, data, dataSi64Shared).get();
    } else {
        enc.remoteIntMatrix(runtime, dataSi64Shared).get();
    }

    sbMatrix dataSbShared(TEST_SIZE,64);
    if (role == 0) {
        enc.localBinMatrix(runtime, data, dataSbShared).get();
    } else {
        enc.remoteBinMatrix(runtime, dataSbShared).get();
    }
    
    si64Matrix resSi64Shared(TEST_SIZE,1);
    odd_even_merge_sort(dataSi64Shared, resSi64Shared, role, enc, eval, runtime);  

    sbMatrix resSbShared(TEST_SIZE,64);
    odd_even_merge_sort(dataSbShared, resSbShared, role, enc, eval, runtime);  
    
    i64Matrix resPlain_1(TEST_SIZE,1);
    enc.revealAll(runtime, resSi64Shared, resPlain_1).get();
    if (role == 0) {
        check_result("Odd Even Merge Sort Test(si64)",resPlain_1, res);
    }

    i64Matrix resPlain_2(TEST_SIZE,1);
    enc.revealAll(runtime, resSbShared, resPlain_2).get();
    if (role == 0) {
        check_result("Odd Even Merge Sort Test(sb)",resPlain_2, res);
    }

    return 0;

}

int join_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }
    
    if (role == 0) {
        debug_info("RUN JOIN TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    //distribute_setup((u64)role, ios, enc, eval, runtime);

    std::vector<i64Matrix> T_1_key(2);
    std::vector<i64Matrix> T_1_other(2);
    std::vector<i64Matrix> T_2_key(2);
    std::vector<i64Matrix> T_2_other(1);
    std::vector<i64Matrix> T_join(5);
    // std::vector<i64Matrix> T_1_auged(4);
    // std::vector<i64Matrix> T_2_auged(4);
    size_t len1 = 5;
    size_t len2 = 6;
    size_t m = 5 ;
    T_1_key[0].resize(len1, 1);
    T_1_key[1].resize(len1, 1);
    T_1_other[0].resize(len1, 1);
    T_1_other[1].resize(len1, 1);
    T_2_key[0].resize(len2, 1);
    T_2_key[1].resize(len2, 1);
    T_2_other[0].resize(len2, 1);
    for(size_t i = 0; i < 5; i++){
        T_join[i].resize(m, 1);
    }

    //set T_1
    T_1_key[0](0,0)=1;
    T_1_key[0](1,0)=1;
    T_1_key[0](2,0)=2;
    T_1_key[0](3,0)=3;
    T_1_key[0](4,0)=4;
    T_1_key[1](0,0)=2;
    T_1_key[1](1,0)=2;
    T_1_key[1](2,0)=2;
    T_1_key[1](3,0)=3;
    T_1_key[1](4,0)=1;

    T_1_other[0](0,0)=14;
    T_1_other[0](1,0)=2;
    T_1_other[0](2,0)=4;
    T_1_other[0](3,0)=4;
    T_1_other[0](4,0)=2;
    T_1_other[1](0,0)=15;
    T_1_other[1](1,0)=4;
    T_1_other[1](2,0)=6;
    T_1_other[1](3,0)=6;
    T_1_other[1](4,0)=0;


    //set T_2
    T_2_key[0](0,0)=1;
    T_2_key[0](1,0)=2;
    T_2_key[0](2,0)=3;
    T_2_key[0](3,0)=3;
    T_2_key[0](4,0)=4;
    T_2_key[0](5,0)=4;
    T_2_key[1](0,0)=2;
    T_2_key[1](1,0)=3;
    T_2_key[1](2,0)=3;
    T_2_key[1](3,0)=3;
    T_2_key[1](4,0)=1;
    T_2_key[1](5,0)=2;

    T_2_other[0](0,0)=0;
    T_2_other[0](1,0)=1;
    T_2_other[0](2,0)=4;
    T_2_other[0](3,0)=5;
    T_2_other[0](4,0)=2;
    T_2_other[0](5,0)=0;


    // //set T_join
    T_join[0](0,0)=1;
    T_join[0](1,0)=1;
    T_join[0](2,0)=3;
    T_join[0](3,0)=3;
    T_join[0](4,0)=4;
    T_join[1](0,0)=2;
    T_join[1](1,0)=2;
    T_join[1](2,0)=3;
    T_join[1](3,0)=3;
    T_join[1](4,0)=1;
    T_join[2](0,0)=2;
    T_join[2](1,0)=14;
    T_join[2](2,0)=4;
    T_join[2](3,0)=4;
    T_join[2](4,0)=2;
    T_join[3](0,0)=4;
    T_join[3](1,0)=15;
    T_join[3](2,0)=6;
    T_join[3](3,0)=6;
    T_join[3](4,0)=0;
    T_join[4](0,0)=0;
    T_join[4](1,0)=0;
    T_join[4](2,0)=4;
    T_join[4](3,0)=5;
    T_join[4](4,0)=2;

    std::vector<si64Matrix> T_1_keyShared(2);
    std::vector<si64Matrix> T_1_otherShared(2);
    std::vector<si64Matrix> T_2_keyShared(2);
    std::vector<si64Matrix> T_2_otherShared(1);
    T_1_keyShared[0].resize(len1, 1);
    T_1_keyShared[1].resize(len1, 1);
    T_1_otherShared[0].resize(len1, 1);
    T_1_otherShared[1].resize(len1, 1);
    T_2_keyShared[0].resize(len2, 1);
    T_2_keyShared[1].resize(len2, 1);
    T_2_otherShared[0].resize(len2, 1);
  
    if (role == 0) {
        enc.localIntMatrix(runtime, T_1_key[0], T_1_keyShared[0]).get();
        enc.localIntMatrix(runtime, T_1_key[1], T_1_keyShared[1]).get();
        enc.localIntMatrix(runtime, T_1_other[0], T_1_otherShared[0]).get();
        enc.localIntMatrix(runtime, T_1_other[1], T_1_otherShared[1]).get();
        enc.localIntMatrix(runtime, T_2_key[0], T_2_keyShared[0]).get();
        enc.localIntMatrix(runtime, T_2_key[1], T_2_keyShared[1]).get();
        enc.localIntMatrix(runtime, T_2_other[0], T_2_otherShared[0]).get();
    } else {
        enc.remoteIntMatrix(runtime, T_1_keyShared[0]).get();
        enc.remoteIntMatrix(runtime, T_1_keyShared[1]).get();
        enc.remoteIntMatrix(runtime, T_1_otherShared[0]).get();
        enc.remoteIntMatrix(runtime, T_1_otherShared[1]).get();
        enc.remoteIntMatrix(runtime, T_2_keyShared[0]).get();
        enc.remoteIntMatrix(runtime, T_2_keyShared[1]).get();
        enc.remoteIntMatrix(runtime, T_2_otherShared[0]).get();
    }

    std::vector<si64Matrix> T_join_Shared(5);
    join(role, T_1_keyShared, T_1_otherShared, T_2_keyShared, T_2_otherShared, T_join_Shared, enc, eval, runtime);

    std::vector<i64Matrix> T_join_Plain(5);
    for(size_t i = 0; i < 5; i++){
        T_join_Plain[i].resize(m, 1);
        enc.revealAll(runtime, T_join_Shared[i], T_join_Plain[i]).get();
    }

    if (role == 0) {
        bool check_flag = true;
        for (size_t i = 0; i < 5; i++) {
            for (size_t j = 0; j < m; j++) {
                if (T_join_Plain[i](j,0) != T_join[i](j, 0)) {
                    check_flag = false; 
                }
            }

        }
        if (check_flag) {
            debug_info("\033[32m JOIN CHECK SUCCESS ! \033[0m\n");
        } else {
            debug_info("\033[31m JOIN CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            for (size_t i = 0; i < 5; i++) {
                debug_output_matrix(T_join[i]);
            }
            debug_info("Func result: \n");
            for (size_t i = 0; i < 5; i++) {
                debug_output_matrix(T_join_Plain[i]);
            }
        }
        //DEBUG
            // for (size_t i = 0; i < 5; i++) {
            //     debug_output_matrix(T_join_Plain[i]);
            // }
    }

    
    // //DEBUG
    // int m = output_size(0, 0);
    // std::vector<i64Matrix> T_2_aligned_Plain(2);
    // for(size_t i = 0; i < 2; i++){
    //     T_2_aligned_Plain[i].resize(m, 1);
    //     enc.revealAll(runtime, T_2_aligned[i], T_2_aligned_Plain[i]).get();
    // }
    // if(role==1){
    //     std::cout<<"j_plain: "<<std::endl;
    //     for(size_t i=0; i<m; i++){
    //         std::cout<<T_2_aligned_Plain[0](i, 0)<<" ";
    //     }
    //     std::cout<<std::endl;
    //     std::cout<<"d_plain: "<<std::endl;
    //     for(size_t i=0; i<m; i++){
    //         std::cout<<T_2_aligned_Plain[1](i, 0)<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    
    // //----


    return 0;

}

int filter_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }
    if (role == 0) {
        debug_info("RUN FILTER TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    std::vector<i64Matrix> T(3);
    size_t len = 6;
    T[0].resize(len, 1);
    T[1].resize(len, 1);
    T[2].resize(len, 1);
    
    T[0](0,0)=5;
    T[0](1,0)=5;
    T[0](2,0)=6;
    T[0](3,0)=6;
    T[0](4,0)=6;
    T[0](5,0)=6;
    
    
    for(size_t i=0; i<len; i++){
        T[1](i,0)=i;
        T[2](i,0)=i;
    }

    std::vector<si64Matrix> T_Shared(3);
    for(size_t i=0; i<3; i++){
        T_Shared[i].resize(len, 1);
        if (role == 0) {
            enc.localIntMatrix(runtime, T[i], T_Shared[i]).get();
        } else {
            enc.remoteIntMatrix(runtime, T_Shared[i]).get();
        }
    }

    std::vector<si64Matrix> T_Filtered(3);
    filter(role, T_Shared, 0, 6, true, "<", T_Filtered, enc, eval, runtime);
    //filter(role, T_Shared, 0, 1, false, "<=", T_Filtered, enc, eval, runtime);

    std::vector<i64Matrix> T_Filtered_Plain(3);
    for(size_t i = 0; i < 3; i++){
        T_Filtered_Plain[i].resize(len, 1);
        enc.revealAll(runtime, T_Filtered[i], T_Filtered_Plain[i]).get();
    }

    if (role == 0) {
        for(size_t i=0; i<3; i++){
            debug_output_matrix(T_Filtered_Plain[i]);
        }
    }

    return 0;
}

int secret_rshift64_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }
    if (role == 0) {
        debug_info("RUN SECRET_RSHIFT64 TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    i64Matrix q(2, 1);
    i64Matrix k(2, 1);
    i64Matrix y(2, 1);
    
    q(0,0)=15;
    q(1,0)=24;

    k(0,0)=1;
    k(1,0)=4;

    y(0,0)=15;
    y(1,0)=6;

    sbMatrix q_Shared(2, 64);
    sbMatrix k_Shared(2, 6);
    sbMatrix y_Shared(2, 64);

    if (role == 0) {
        enc.localBinMatrix(runtime, q, q_Shared).get();
        enc.localBinMatrix(runtime, k, k_Shared).get();
      
    } else {
        enc.remoteBinMatrix(runtime, q_Shared).get();
        enc.remoteBinMatrix(runtime, k_Shared).get();
    }

    //bool_cipher_secret_rshift64(role, q_Shared, k_Shared, y_Shared, eval, runtime);
    bool_cipher_div_pow2(role, q_Shared, k_Shared, y_Shared, eval, runtime);
    i64Matrix y_test(2, 1);
    enc.revealAll(runtime, y_Shared, y_test).get();
    if (role == 0) {
        check_result("Secrect shift Test",y_test, y);
    }
    return 0;

}

int semi_join_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN SEMI JOIN TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // T_1: 5 rows, 2 key cols, 1 other col
    // T_2: 6 rows, 2 key cols, 1 other col (other not used in semi-join result)
    size_t len1 = 5, len2 = 6;

    std::vector<i64Matrix> T_1_key(2), T_1_other(1);
    std::vector<i64Matrix> T_2_key(2), T_2_other(1);

    T_1_key[0].resize(len1, 1); T_1_key[1].resize(len1, 1);
    T_1_other[0].resize(len1, 1);
    T_2_key[0].resize(len2, 1); T_2_key[1].resize(len2, 1);
    T_2_other[0].resize(len2, 1);

    // T_1: (key0, key1, other0)
    // row0: (1, 2, 14)  -> key (1,2) in T_2 -> KEEP
    // row1: (1, 2,  2)  -> key (1,2) in T_2 -> KEEP
    // row2: (2, 2,  4)  -> key (2,2) NOT in T_2 -> ZERO
    // row3: (3, 3,  4)  -> key (3,3) in T_2 -> KEEP
    // row4: (4, 1,  2)  -> key (4,1) in T_2 -> KEEP
    T_1_key[0](0,0)=1; T_1_key[0](1,0)=1; T_1_key[0](2,0)=2; T_1_key[0](3,0)=3; T_1_key[0](4,0)=4;
    T_1_key[1](0,0)=2; T_1_key[1](1,0)=2; T_1_key[1](2,0)=2; T_1_key[1](3,0)=3; T_1_key[1](4,0)=1;
    T_1_other[0](0,0)=2; T_1_other[0](1,0)=14; T_1_other[0](2,0)=4; T_1_other[0](3,0)=4; T_1_other[0](4,0)=2;

    // T_2: keys that exist
    // (1,2), (2,3), (3,3), (3,3), (4,1), (4,2)
    T_2_key[0](0,0)=1; T_2_key[0](1,0)=2; T_2_key[0](2,0)=3; T_2_key[0](3,0)=3; T_2_key[0](4,0)=4; T_2_key[0](5,0)=4;
    T_2_key[1](0,0)=2; T_2_key[1](1,0)=3; T_2_key[1](2,0)=3; T_2_key[1](3,0)=3; T_2_key[1](4,0)=1; T_2_key[1](5,0)=2;
    T_2_other[0](0,0)=0; T_2_other[0](1,0)=1; T_2_other[0](2,0)=4; T_2_other[0](3,0)=5; T_2_other[0](4,0)=2; T_2_other[0](5,0)=0;

    // Expected semi-join result: T_1 rows where key exists in T_2, else zeroed
    // result cols: key0, key1, other0
    std::vector<i64Matrix> T_expected(3);
    for(int i = 0; i < 3; i++) T_expected[i].resize(len1, 1);
    // row0: kept (1,2,14)
    T_expected[0](0,0)=1; T_expected[1](0,0)=2; T_expected[2](0,0)=2;
    // row1: kept (1,2,2)
    T_expected[0](1,0)=1; T_expected[1](1,0)=2; T_expected[2](1,0)=14;
    // row2: zeroed (2,2) not in T_2
    T_expected[0](2,0)=0; T_expected[1](2,0)=0; T_expected[2](2,0)=0;
    // row3: kept (3,3,4)
    T_expected[0](3,0)=3; T_expected[1](3,0)=3; T_expected[2](3,0)=4;
    // row4: kept (4,1,2)
    T_expected[0](4,0)=4; T_expected[1](4,0)=1; T_expected[2](4,0)=2;

    // Encrypt
    std::vector<si64Matrix> T_1_keyShared(2), T_1_otherShared(1);
    std::vector<si64Matrix> T_2_keyShared(2), T_2_otherShared(1);
    T_1_keyShared[0].resize(len1, 1); T_1_keyShared[1].resize(len1, 1);
    T_1_otherShared[0].resize(len1, 1);
    T_2_keyShared[0].resize(len2, 1); T_2_keyShared[1].resize(len2, 1);
    T_2_otherShared[0].resize(len2, 1);

    if (role == 0) {
        enc.localIntMatrix(runtime, T_1_key[0], T_1_keyShared[0]).get();
        enc.localIntMatrix(runtime, T_1_key[1], T_1_keyShared[1]).get();
        enc.localIntMatrix(runtime, T_1_other[0], T_1_otherShared[0]).get();
        enc.localIntMatrix(runtime, T_2_key[0], T_2_keyShared[0]).get();
        enc.localIntMatrix(runtime, T_2_key[1], T_2_keyShared[1]).get();
        enc.localIntMatrix(runtime, T_2_other[0], T_2_otherShared[0]).get();
    } else {
        enc.remoteIntMatrix(runtime, T_1_keyShared[0]).get();
        enc.remoteIntMatrix(runtime, T_1_keyShared[1]).get();
        enc.remoteIntMatrix(runtime, T_1_otherShared[0]).get();
        enc.remoteIntMatrix(runtime, T_2_keyShared[0]).get();
        enc.remoteIntMatrix(runtime, T_2_keyShared[1]).get();
        enc.remoteIntMatrix(runtime, T_2_otherShared[0]).get();
    }

    // Run semi_join
    std::vector<si64Matrix> T_semi_joined;
    semi_join(role, T_1_keyShared, T_1_otherShared, T_2_keyShared, T_2_otherShared, T_semi_joined, enc, eval, runtime);

    // Reveal results
    std::vector<i64Matrix> T_result(3);
    for(int i = 0; i < 3; i++){
        T_result[i].resize(len1, 1);
        enc.revealAll(runtime, T_semi_joined[i], T_result[i]).get();
    }

    // Check
    if (role == 0) {
        bool check_flag = true;
        for (int i = 0; i < 3; i++) {
            for (size_t j = 0; j < len1; j++) {
                if (T_result[i](j, 0) != T_expected[i](j, 0)) {
                    check_flag = false;
                }
            }
        }
        if (check_flag) {
            debug_info("\033[32m SEMI JOIN CHECK SUCCESS ! \033[0m\n");
        } else {
            debug_info("\033[31m SEMI JOIN CHECK ERROR ! \033[0m\n");
            debug_info("Expected result: ");
            for (int i = 0; i < 3; i++) {
                debug_output_matrix(T_expected[i]);
            }
            debug_info("Actual result: ");
            for (int i = 0; i < 3; i++) {
                debug_output_matrix(T_result[i]);
            }
        }
    }

    return 0;
}

int mul_and_sum_test(CLP &cmd){
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN MUL_AND_SUM TEST");
    }

    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // ---------- 子用例 1: 固定数据手算验证 ----------
    // a = [1, 2, 3, 4, 5], b = [10, 20, 30, 40, 50]
    // 期望 inner product = 1*10 + 2*20 + 3*30 + 4*40 + 5*50 = 10+40+90+160+250 = 550
    {
        size_t TEST_SIZE = 5;
        i64Matrix a_plain(TEST_SIZE, 1), b_plain(TEST_SIZE, 1);
        a_plain(0,0)=1; a_plain(1,0)=2; a_plain(2,0)=3; a_plain(3,0)=4; a_plain(4,0)=5;
        b_plain(0,0)=10; b_plain(1,0)=20; b_plain(2,0)=30; b_plain(3,0)=40; b_plain(4,0)=50;

        i64 expected = 0;
        for (size_t i = 0; i < TEST_SIZE; i++) {
            expected += a_plain(i, 0) * b_plain(i, 0);
        }

        si64Matrix a_shared(TEST_SIZE, 1), b_shared(TEST_SIZE, 1);
        if (role == 0) {
            enc.localIntMatrix(runtime, a_plain, a_shared).get();
            enc.localIntMatrix(runtime, b_plain, b_shared).get();
        } else {
            enc.remoteIntMatrix(runtime, a_shared).get();
            enc.remoteIntMatrix(runtime, b_shared).get();
        }

        si64Matrix sum_shared;
        mul_and_sum(role, a_shared, b_shared, sum_shared, enc, eval, runtime);

        i64Matrix sum_plain(1, 1);
        enc.revealAll(runtime, sum_shared, sum_plain).get();

        if (role == 0) {
            i64Matrix expected_mat(1, 1);
            expected_mat(0, 0) = expected;
            check_result("MulAndSum Test - small fixed", sum_plain, expected_mat);
        }
    }

    // ---------- 子用例 2: 含负数 / 跨零的随机数据 ----------
    {
        size_t TEST_SIZE = 64;
        i64Matrix a_plain(TEST_SIZE, 1), b_plain(TEST_SIZE, 1);
        std::mt19937_64 rng(12345);  // 固定种子，便于复现
        for (size_t i = 0; i < TEST_SIZE; i++) {
            a_plain(i, 0) = (i64)(rng() % 2001) - 1000;  // [-1000, 1000]
            b_plain(i, 0) = (i64)(rng() % 2001) - 1000;
        }
        i64 expected = 0;
        for (size_t i = 0; i < TEST_SIZE; i++) {
            expected += a_plain(i, 0) * b_plain(i, 0);
        }

        si64Matrix a_shared(TEST_SIZE, 1), b_shared(TEST_SIZE, 1);
        if (role == 0) {
            enc.localIntMatrix(runtime, a_plain, a_shared).get();
            enc.localIntMatrix(runtime, b_plain, b_shared).get();
        } else {
            enc.remoteIntMatrix(runtime, a_shared).get();
            enc.remoteIntMatrix(runtime, b_shared).get();
        }

        si64Matrix sum_shared;
        mul_and_sum(role, a_shared, b_shared, sum_shared, enc, eval, runtime);

        i64Matrix sum_plain(1, 1);
        enc.revealAll(runtime, sum_shared, sum_plain).get();

        if (role == 0) {
            i64Matrix expected_mat(1, 1);
            expected_mat(0, 0) = expected;
            check_result("MulAndSum Test - random signed (n=64)", sum_plain, expected_mat);
        }
    }

    // ---------- 子用例 3: 与 cipher_mul + 累加 一致性对比 ----------
    {
        size_t TEST_SIZE = 32;
        i64Matrix a_plain(TEST_SIZE, 1), b_plain(TEST_SIZE, 1);
        std::mt19937_64 rng(0xC0FFEEULL);
        for (size_t i = 0; i < TEST_SIZE; i++) {
            a_plain(i, 0) = (i64)(rng() % 1000);
            b_plain(i, 0) = (i64)(rng() % 1000);
        }

        si64Matrix a_shared(TEST_SIZE, 1), b_shared(TEST_SIZE, 1);
        if (role == 0) {
            enc.localIntMatrix(runtime, a_plain, a_shared).get();
            enc.localIntMatrix(runtime, b_plain, b_shared).get();
        } else {
            enc.remoteIntMatrix(runtime, a_shared).get();
            enc.remoteIntMatrix(runtime, b_shared).get();
        }

        // 路径 A: 优化版 mul_and_sum
        si64Matrix sum_fast;
        mul_and_sum(role, a_shared, b_shared, sum_fast, enc, eval, runtime);
        i64Matrix sum_fast_plain(1, 1);
        enc.revealAll(runtime, sum_fast, sum_fast_plain).get();

        // 路径 B: cipher_mul 得到逐元素积，再在密文上累加成 1×1
        si64Matrix prod_shared(TEST_SIZE, 1);
        cipher_mul(role, a_shared, b_shared, prod_shared, eval, enc, runtime);
        si64Matrix sum_ref(1, 1);
        sum_ref.mShares[0](0, 0) = 0;
        sum_ref.mShares[1](0, 0) = 0;
        for (size_t i = 0; i < TEST_SIZE; i++) {
            sum_ref.mShares[0](0, 0) += prod_shared.mShares[0](i, 0);
            sum_ref.mShares[1](0, 0) += prod_shared.mShares[1](i, 0);
        }
        i64Matrix sum_ref_plain(1, 1);
        enc.revealAll(runtime, sum_ref, sum_ref_plain).get();

        if (role == 0) {
            check_result("MulAndSum Test - matches cipher_mul+sum", sum_fast_plain, sum_ref_plain);
        }
    }

    // ---------- 子用例 4: 边界 —— 长度为 1 的情形 ----------
    {
        size_t TEST_SIZE = 1;
        i64Matrix a_plain(TEST_SIZE, 1), b_plain(TEST_SIZE, 1);
        a_plain(0, 0) = -7;
        b_plain(0, 0) = 11;

        si64Matrix a_shared(TEST_SIZE, 1), b_shared(TEST_SIZE, 1);
        if (role == 0) {
            enc.localIntMatrix(runtime, a_plain, a_shared).get();
            enc.localIntMatrix(runtime, b_plain, b_shared).get();
        } else {
            enc.remoteIntMatrix(runtime, a_shared).get();
            enc.remoteIntMatrix(runtime, b_shared).get();
        }

        si64Matrix sum_shared;
        mul_and_sum(role, a_shared, b_shared, sum_shared, enc, eval, runtime);

        i64Matrix sum_plain(1, 1);
        enc.revealAll(runtime, sum_shared, sum_plain).get();

        if (role == 0) {
            i64Matrix expected_mat(1, 1);
            expected_mat(0, 0) = -77;
            check_result("MulAndSum Test - n=1 boundary", sum_plain, expected_mat);
        }
    }

    return 0;
}