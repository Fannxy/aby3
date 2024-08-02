#include "Test.h"

#include "../aby3-Basic/Search.h"
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
    for(int i = 0; i < n; ++i){
        ref_res(i, 0) = 0;
        for(int j = 0; j < m; ++j){
            a(i*m + j, 0) = i + j;
            b(i*m + j, 0) = i + j;
            ref_res(i, 0) += (i + j) * (i + j);
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


    // compute the dot product.
    si64Matrix enc_res(n, 1);
    constant_sint_dot(role, enc_a_vec, enc_b_vec, enc_res, enc, eval, runtime);

    // decrypt the result.
    i64Matrix res;
    enc.revealAll(runtime, enc_res, res).get();

    // check the result.
    if(role == 0){
        bool check_flag = check_result("Sint dot", res, ref_res);
    }

    return 0;
}

int binary_search_test(oc::CLP& cmd){
    BASIC_TEST_INIT

    if (role == 0) {
        debug_info("RUN Binary Search TEST");
    }

    // prepare the test data, keyset, data, input key.
    int n = 1 << 3;
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
    mtagBS(data, keyset, key_mat, res_mat, role, enc, eval, runtime);
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

    return 0;
}