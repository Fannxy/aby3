#include "Test.h"

#include "../aby3-Basic/Search.h"
#include "../aby3-Benchmark/benchmark.h"
#include "../aby3-RTR/PtATests.h"
#include "../aby3-Basic/Basics.h"
#include "../aby3-RTR/BuildingBlocks.h"

using namespace oc;
using namespace aby3;
using namespace std;

int constant_dot_test(oc::CLP& cmd){

    BASIC_TEST_INIT

    if (role == 0) {
        debug_info("RUN Constant dot TEST");
    }

    int n = 5, m = 10;

    // compute the dot product of two n*m matrices.
    i64Matrix a(n*m, 1), b(n*m, 1);
    i64Matrix ref_res(n, 1);
    i64Matrix ref_res_bool(n, 1);
    for(int i = 0; i < n; ++i){
        ref_res(i, 0) = 0;
        ref_res_bool(i, 0) = 0;
        for(int j = 0; j < m; ++j){
            a(i*m + j, 0) = i + j;
            b(i*m + j, 0) = i + j;
            ref_res(i, 0) += (i + j) * (i + j);
            ref_res_bool(i, 0) ^= (i + j) & (i + j);
        }
    }
    
    // generate the ciphertext data.
    si64Matrix enc_a(n*m, 1), enc_b(n*m, 1);
    if(role == 0){
        enc.localIntMatrix(runtime, a, enc_a).get();
        enc.localIntMatrix(runtime, b, enc_b).get();
    }
    else{
        enc.remoteIntMatrix(runtime, enc_a).get();
        enc.remoteIntMatrix(runtime, enc_b).get();
    }
    std::vector<si64Matrix> enc_a_vec(n);
    std::vector<si64Matrix> enc_b_vec(n);
    for(int i = 0; i < n; ++i){
        enc_a_vec[i].resize(m, 1); enc_b_vec[i].resize(m, 1);
        for(int j=0; j<m; j++){
            enc_a_vec[i].mShares[0](j, 0) = enc_a.mShares[0](i*m + j, 0);
            enc_a_vec[i].mShares[1](j, 0) = enc_a.mShares[1](i*m + j, 0);
            enc_b_vec[i].mShares[0](j, 0) = enc_b.mShares[0](i*m + j, 0);
            enc_b_vec[i].mShares[1](j, 0) = enc_b.mShares[1](i*m + j, 0);
        }
    }

    // generate the ciphertext data for bool vase.
    sbMatrix enc_a_bool(n*m, 64), enc_b_bool(n*m, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, a, enc_a_bool).get();
        enc.localBinMatrix(runtime, b, enc_b_bool).get();
    }
    else{
        enc.remoteBinMatrix(runtime, enc_a_bool).get();
        enc.remoteBinMatrix(runtime, enc_b_bool).get();
    }
    std::vector<sbMatrix> enc_a_vec_bool(n);
    std::vector<sbMatrix> enc_b_vec_bool(n);
    for(int i = 0; i < n; ++i){
        enc_a_vec_bool[i].resize(m, 16); enc_b_vec_bool[i].resize(m, 16);
        for(int j=0; j<m; j++){
            enc_a_vec_bool[i].mShares[0](j, 0) = enc_a_bool.mShares[0](i*m + j, 0);
            enc_a_vec_bool[i].mShares[1](j, 0) = enc_a_bool.mShares[1](i*m + j, 0);
            enc_b_vec_bool[i].mShares[0](j, 0) = enc_b_bool.mShares[0](i*m + j, 0);
            enc_b_vec_bool[i].mShares[1](j, 0) = enc_b_bool.mShares[1](i*m + j, 0);
        }
    }


    // compute the dot product.
    si64Matrix enc_res(n, 1);
    constant_sint_dot(role, enc_a_vec, enc_b_vec, enc_res, enc, eval, runtime);

    // compute the dot product for the bool case.
    sbMatrix enc_res_bool(n, 64);
    constant_bool_dot(role, enc_a_vec_bool, enc_b_vec_bool, enc_res_bool, enc, eval, runtime);

    // decrypt the result.
    i64Matrix res;
    enc.revealAll(runtime, enc_res, res).get();
    i64Matrix res_bool;
    enc.revealAll(runtime, enc_res_bool, res_bool).get();

    // check the result.
    if(role == 0){
        bool check_flag = check_result("Sint dot", res, ref_res);
        check_flag = check_result("Bool dot", res_bool, ref_res_bool);
    }

    return 0;
}

int binary_search_test(oc::CLP& cmd){
    BASIC_TEST_INIT

    if (role == 0) {
        debug_info("RUN Binary Search TEST");
    }

    // prepare the test data, keyset, data, input key.
    int n = 1 << 8;
    std::vector<sbMatrix> keyset(n);
    std::vector<si64Matrix> data(n);

    i64Matrix key_and_data(n, 1);
    for(int i=0; i<n; i++){
        key_and_data(i, 0) = i;
    }
    sbMatrix keyset_mat(n, 16);
    si64Matrix data_mat(n, 1);

    if(role == 0){
        enc.localBinMatrix(runtime, key_and_data, keyset_mat).get();
        enc.localIntMatrix(runtime, key_and_data, data_mat).get();
    }
    else{
        enc.remoteBinMatrix(runtime, keyset_mat).get();
        enc.remoteIntMatrix(runtime, data_mat).get();
    }

    for(int i=0; i<n; i++){
        keyset[i].resize(1, 16);
        data[i].resize(1, 1);
        keyset[i].mShares[0](0, 0) = keyset_mat.mShares[0](i, 0);
        keyset[i].mShares[1](0, 0) = keyset_mat.mShares[1](i, 0);
        data[i].mShares[0](0, 0) = data_mat.mShares[0](i, 0);
        data[i].mShares[1](0, 0) = data_mat.mShares[1](i, 0);
    }

    // prepare the input key and the reference result.
    i64Matrix key(1, 1);
    key(0, 0) = n-2;
    sbMatrix key_mat(1, 16);
    i64Matrix ref_res(1, 1); ref_res(0, 0) = n-2;
    i64Matrix ref_mres(n, 1); 
    for(int i=0; i<n; i++){
        ref_mres(i, 0) = (i == n-2);
    }
    if(role == 0){
        enc.localBinMatrix(runtime, key, key_mat).get();
    }
    else{
        enc.remoteBinMatrix(runtime, key_mat).get();
    }

    sbMatrix res_mat(n, 1);
    i64Matrix res;
    // test the mBSs, mcompBS.
    mcompBS(keyset, key_mat, res_mat, role, enc, eval, runtime);
    enc.revealAll(runtime, res_mat, res).get();
    if(role == 0){
        bool check_flag = check_result("mcompBS", res, ref_mres);
    }

    // test the mBSs, mtagBS.
    mtagBS(keyset, key_mat, res_mat, role, enc, eval, runtime);
    enc.revealAll(runtime, res_mat, res).get();
    if(role == 0){
        bool check_flag = check_result("mtagBS", res, ref_mres);
    }

    // // test the BSs, compBS.
    si64Matrix res_mat2;
    i64Matrix res2;
    compBS(data, keyset, key_mat, res_mat2, role, enc, eval, runtime);
    enc.revealAll(runtime, res_mat2, res2).get();
    if(role == 0){
        bool check_flag = check_result("compBS", res2, ref_res);
    }

    // test the BSs, tagBS.
    tagBS(data, keyset, key_mat, res_mat2, role, enc, eval, runtime);
    enc.revealAll(runtime, res_mat2, res2).get();
    if(role == 0){
        bool check_flag = check_result("tagBS", res2, ref_res);
    }

    int alpha = sqrtToPowerOfTwo(n);
    subHBS(data, keyset, key_mat, res_mat2, role, enc, eval, runtime, alpha);
    enc.revealAll(runtime, res_mat2, res2).get();
    if(role == 0){
        bool check_flag = check_result("subHBS", res2, ref_res);    
    }

    return 0;
}

int data_generation_test(oc::CLP& cmd){
    SETUP_PROCESS

    if (role == 0 && rank == 0) {   
        debug_info("RUN Data Generation TEST");
    }

    int n = 1 << 8;
    std::vector<si64Matrix> data;
    std::vector<sbMatrix> data_bool;

    if(rank == 0){
        generate_data(data, n, 1, role, enc, eval, runtime);
        generate_data(data_bool, n, 1, 32, role, enc, eval, runtime);
    }

    std::vector<si64Matrix> p_data;
    std::vector<sbMatrix> p_data_bool;
    generate_data_parallel(p_data, n, 1, role, enc, eval, runtime);
    generate_data_parallel(p_data_bool, n, 1, 32, role, enc, eval, runtime);

    // test whether datas are the same.
    if(rank == 0){
        i64Matrix test_data(n, 1);
        i64Matrix test_data_bool(n, 1);
        i64Matrix test_p_data(n, 1);
        i64Matrix test_p_data_bool(n, 1);

        si64Matrix data_mat(n, 1);
        sbMatrix data_bool_mat(n, 32);
        si64Matrix p_data_mat(n, 1);
        sbMatrix p_data_bool_mat(n, 32);

        for(int i=0; i<n; i++){
            data_mat.mShares[0](i, 0) = data[i].mShares[0](0, 0);
            data_mat.mShares[1](i, 0) = data[i].mShares[1](0, 0);
            data_bool_mat.mShares[0](i, 0) = data_bool[i].mShares[0](0, 0);
            data_bool_mat.mShares[1](i, 0) = data_bool[i].mShares[1](0, 0);
            p_data_mat.mShares[0](i, 0) = p_data[i].mShares[0](0, 0);
            p_data_mat.mShares[1](i, 0) = p_data[i].mShares[1](0, 0);
            p_data_bool_mat.mShares[0](i, 0) = p_data_bool[i].mShares[0](0, 0);
            p_data_bool_mat.mShares[1](i, 0) = p_data_bool[i].mShares[1](0, 0);
        }

        enc.revealAll(runtime, data_mat, test_data).get();
        enc.revealAll(runtime, data_bool_mat, test_data_bool).get();
        enc.revealAll(runtime, p_data_mat, test_p_data).get();
        enc.revealAll(runtime, p_data_bool_mat, test_p_data_bool).get();

        check_result("generate data parallel", test_data, test_p_data);
        check_result("generate data bool parallel", test_data_bool, test_p_data_bool);
    }

    return 0;
}

int mpi_binary_search_test(oc::CLP& cmd){

    SETUP_PROCESS

    if (role == 0 && rank == 0) {
        debug_info("RUN MPI Binary Search TEST");
    }

    // prepare the test data, keyset, data, input key.
    std::vector<sbMatrix> keyset;
    std::vector<si64Matrix> data;
    sbMatrix starget(1, 32);

    size_t n=1, m=64, optB=16, k=1;

    // prepare the pta tasks.
    auto ptaBSTask = new ABY3MPITask<sb64, sb64, si64, si64, PtABS>(size, optB, role, enc, runtime, eval);
    auto ptaMBSTask = new ABY3MPIPairOnlyTask<sb64, sb64, sb64, sb64, PtAMBS>(size, optB, role, enc, runtime, eval);

    i64Matrix target(1, 1); target(0, 0) = m - 2;
    i64Matrix target_m(m, 1);
    for(size_t i=0; i<m; i++){
        target_m(i, 0) = (i == m-2);
    }

    if(rank == 0){
        // data preparation.
        i64Matrix plain_data(m, 1);
        for(size_t i=0; i<m; i++){
            plain_data(i, 0) = i;
        }
        si64Matrix sdata(m, 1);
        sbMatrix skeyset(m, 32);

        if(role == 0){
            enc.localIntMatrix(runtime, plain_data, sdata).get();
            enc.localBinMatrix(runtime, plain_data, skeyset).get();
            enc.localBinMatrix(runtime, target, starget).get();
        }
        else{
            enc.remoteIntMatrix(runtime, sdata).get();
            enc.remoteBinMatrix(runtime, skeyset).get();
            enc.remoteBinMatrix(runtime, starget).get();
        }
        // get the vector data.
        for(size_t i=0; i<m; i++){
            aby3::si64Matrix _data(1, 1);
            aby3::sbMatrix _keyset(1, 32);
            _data.mShares[0](0, 0) = sdata.mShares[0](i, 0);
            _data.mShares[1](0, 0) = sdata.mShares[1](i, 0);
            _keyset.mShares[0](0, 0) = skeyset.mShares[0](i, 0);
            _keyset.mShares[1](0, 0) = skeyset.mShares[1](i, 0);
            data.push_back(_data);
            keyset.push_back(_keyset);
        }
    }

    // run the ptaBS.
    si64Matrix bs_res;
    sbMatrix bs_res_m;
    i64Matrix test_bs_res;
    i64Matrix test_bs_res_m;

    ptaBS(data, keyset, starget, bs_res, ptaBSTask);

    if(rank == 0){
        enc.revealAll(runtime, bs_res, test_bs_res).get();
        if(role == 0){
            bool check_flag = check_result("ptaBS", test_bs_res, target);
        }
    }

    ptaMBS(keyset, starget, bs_res_m, ptaMBSTask);

    if(rank == 0){
        enc.revealAll(runtime, bs_res_m, test_bs_res_m).get();
        if(role == 0){
            bool check_flag = check_result("ptaMBS", test_bs_res_m, target_m);
        }
    }

    subHBS_with_PtA(data, keyset, starget, bs_res, role, enc, eval, runtime, 16, 16, ptaBSTask, ptaMBSTask);

    if(rank == 0){
        enc.revealAll(runtime, bs_res, test_bs_res).get();
        if(role == 0){
            bool check_flag = check_result("subHBS_with_PtA", test_bs_res, target);
        }
    }

    return 0;
}