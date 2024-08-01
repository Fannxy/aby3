#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-Basic/Sort.h"
#include "../aby3-RTR/BuildingBlocks.h"


using namespace oc;
using namespace aby3;

int bc_sort_test(oc::CLP &cmd){

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
        debug_info("RUN BC Sort TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    
    // prepare the data.
    aby3::i64Matrix data(10, 1);
    aby3::i64Matrix res(10, 1);
    for(size_t i=0; i<10; i++) {
        data(i, 0) = 9 - i;
        if(i < 5){
            res(i, 0) = 5 + i;
        }
        else{
            res(i, 0) = i - 5;
        }
    }
    
    aby3::sbMatrix enc_data(10, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, data, enc_data).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_data).get();
    }

    std::vector<aby3::sbMatrix> enc_data_vec(10);
    for(size_t i=0; i<10; i++){
        aby3::sbMatrix tmp(1, 64);
        tmp.mShares[0](0, 0) = enc_data.mShares[0](i, 0);
        tmp.mShares[1](0, 0) = enc_data.mShares[1](i, 0);
        enc_data_vec[i] = tmp;
    }

    // sort the data.
    std::vector<size_t> lows = {0, 5};
    std::vector<size_t> highs = {5, 10};

    bc_sort_different(enc_data_vec, lows, highs, role, enc, eval, runtime, 1024);
    aby3::sbMatrix enc_test(10, 64);
    for(size_t i=0; i<10; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    // check the result.
    aby3::i64Matrix test(10, 1);
    enc.revealAll(runtime, enc_test, test).get();


    if(role == 0){
        bool check_flag = check_result("BC Sort Test", test, res);
        if(!check_flag){
            debug_info("test res = ");
            debug_output_matrix(test);
            debug_info("reference res = ");
            debug_output_matrix(res);
        }
    }

    return 0;
}

int bc_sort_corner_test(oc::CLP &cmd){

    BASIC_TEST_INIT

    if(role == 0){
        debug_info("RUN BC Sort Corner TEST");
    }

    // prepare the data.
    aby3::i64Matrix data(10, 1);
    std::vector<int> res_vec = {9, 7, 8, 5, 6, 2, 3, 4, 0, 1};
    aby3::i64Matrix res(10, 1);
    for(size_t i=0; i<10; i++) {
        data(i, 0) = 9 - i;
        res(i, 0) = res_vec[i];
    }

    // sort the data.
    std::vector<size_t> lows = {1, 3, 5, 8};
    std::vector<size_t> highs = {3, 5, 8, 10};

    aby3::sbMatrix enc_data(10, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, data, enc_data).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_data).get();
    }

    std::vector<aby3::sbMatrix> enc_data_vec(10);
    for(size_t i=0; i<10; i++){
        aby3::sbMatrix tmp(1, 64);
        tmp.mShares[0](0, 0) = enc_data.mShares[0](i, 0);
        tmp.mShares[1](0, 0) = enc_data.mShares[1](i, 0);
        enc_data_vec[i] = tmp;
    }


    bc_sort_different(enc_data_vec, lows, highs, role, enc, eval, runtime, 1024);
    aby3::sbMatrix enc_test(10, 64);
    for(size_t i=0; i<10; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    // check the result.
    aby3::i64Matrix test(10, 1);
    enc.revealAll(runtime, enc_test, test).get();
    if(role == 0){
        bool flag = check_result("BC Sort Corner Test", test, res);
        if(!flag){
            debug_info("test res = ");
            debug_output_matrix(test);
            debug_info("reference res = ");
            debug_output_matrix(res);
        }
    }
    

    return 0;
}

int bc_sort_multiple_times(oc::CLP& cmd){

    BASIC_TEST_INIT

    // this test can not pass
    if(role == 0){
        debug_info("RUN BC Sort Multiple Times TEST");
    }

    size_t test_size = 50;
    size_t half_size = test_size / 2;

    // prepare the data.
    aby3::i64Matrix data(test_size, 1);
    aby3::i64Matrix res(test_size, 1);
    for(size_t i=0; i<test_size; i++) {
        data(i, 0) = test_size - 1 - i;
        // res(i, 0) = i;
        if(i < half_size){
            res(i, 0) = half_size + i;
        }
        else{
            res(i, 0) = i - half_size;
        }
    }

    // sort the data.
    std::vector<size_t> lows = {0, half_size};
    std::vector<size_t> highs = {half_size, test_size};

    aby3::sbMatrix enc_data(test_size, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, data, enc_data).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_data).get();
    }

    std::vector<aby3::sbMatrix> enc_data_vec(test_size);
    for(size_t i=0; i<test_size; i++){
        aby3::sbMatrix tmp(1, 64);
        tmp.mShares[0](0, 0) = enc_data.mShares[0](i, 0);
        tmp.mShares[1](0, 0) = enc_data.mShares[1](i, 0);
        enc_data_vec[i] = tmp;
    }

    bc_sort_different(enc_data_vec, lows, highs, role, enc, eval, runtime, 625);

    aby3::sbMatrix enc_test(test_size, 64);
    for(size_t i=0; i<test_size; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    // check the result.
    aby3::i64Matrix test(test_size, 1);
    enc.revealAll(runtime, enc_test, test).get();
    if(role == 0){
        bool flag = check_result("BC Sort Multiple Times Test", test, res);
        if(!flag){
            debug_info("test res = ");
            debug_output_matrix(test);
            debug_info("reference res = ");
            debug_output_matrix(res);
        }
    }

    return 0;
}

int quick_sort_test(oc::CLP &cmd){
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
        debug_info("RUN Quick Sort TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // prepare the data.
    size_t data_size = 1 << 15;
    aby3::i64Matrix data_plain(data_size, 1);
    aby3::i64Matrix data_res(data_size, 1);

    for(size_t i=0; i<data_size; i++){
        data_plain(i, 0) = data_size - i;
        data_res(i, 0) = i + 1;
    }

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

    // sort the data.
    size_t min_size = (1 << 5);
    quick_sort_different(enc_data_vec, role, enc, eval, runtime, min_size);

    aby3::sbMatrix enc_test(data_size, 64);
    for(size_t i=0; i<data_size; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    aby3::i64Matrix test(data_size, 1);
    enc.revealAll(runtime, enc_test, test).get();

    if(role == 0){
        check_result("Quick Sort Test", test, data_res);
    }

    return 0;
}

int quick_sort_with_duplicate_elements_test(oc::CLP& cmd){

    BASIC_TEST_INIT

    if(role == 0){
        debug_info("RUN Quick Sort With Duplicate Elements TEST");
    }

    // prepare the data.
    size_t data_size = 1 << 20;
    size_t bin_size = 1 << 5;
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
    size_t min_size = (1 << 3);

    // tag append, assume that max(data) + log(n) << 64bits, otherwise need to change the tag size.
    quick_sort(enc_data_vec, role, enc, eval, runtime, min_size);
    // tag_append(role, enc_data_vec);


    aby3::sbMatrix enc_test(data_size, 64);
    for(size_t i=0; i<data_size; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    aby3::i64Matrix test(data_size, 1);
    enc.revealAll(runtime, enc_test, test).get();

    if(role == 0){
        check_result("Quick Sort with Duplicated Elements Test", test, data_res);
    }
    return 0;
}

int odd_even_merge_test(oc::CLP& cmd){

    BASIC_TEST_INIT

    if(role == 0){
        debug_info("RUN ODD EVEN MERGE SORT!");
    }

    // prepare the data.
    size_t arr1_len = 50, arr2_len = 98, arr3_len = 50, arr4_len = 98;

    aby3::i64Matrix arr1(arr1_len, 1);
    aby3::i64Matrix arr2(arr2_len, 1);
    aby3::i64Matrix arr3(arr3_len, 1);
    aby3::i64Matrix arr4(arr4_len, 1);

    std::vector<int> ref_res(arr1_len + arr2_len);
    aby3::i64Matrix res_res_(arr1_len + arr2_len, 1);

    std::vector<int> multi_ref_res(arr1_len + arr2_len + arr3_len + arr4_len);
    aby3::i64Matrix multi_res_res_(arr1_len + arr2_len + arr3_len + arr4_len, 1);

    for(size_t i=0; i<arr1_len; i++){
        arr1(i, 0) = i + 10;
        ref_res[i] = i + 10;
        multi_ref_res[i] = i + 10;
    }
    for(size_t i=0; i<arr2_len; i++){
        arr2(i, 0) = i + 20;
        ref_res[i + arr1_len] = i + 20;
        multi_ref_res[i + arr1_len] = i + 20;
    }
    for(size_t i=0; i<arr3_len; i++){
        arr3(i, 0) = i + 30;
        multi_ref_res[i + arr1_len + arr2_len] = i + 30;
    }
    for(size_t i=0; i<arr4_len; i++){
        arr4(i, 0) = i + 40;
        multi_ref_res[i + arr1_len + arr2_len + arr3_len] = i + 40;
    }

    std::sort(ref_res.begin(), ref_res.end());
    for(size_t i=0; i<arr1_len + arr2_len; i++) res_res_(i, 0) = ref_res[i];

    std::sort(multi_ref_res.begin(), multi_ref_res.end());
    for(size_t i=0; i<arr1_len + arr2_len + arr3_len + arr4_len; i++) multi_res_res_(i, 0) = multi_ref_res[i];
    
    // enc the data.
    aby3::sbMatrix enc_arr1(arr1_len, 64);
    aby3::sbMatrix enc_arr2(arr2_len, 64);
    aby3::sbMatrix enc_arr3(arr3_len, 64);
    aby3::sbMatrix enc_arr4(arr4_len, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, arr1, enc_arr1).get();
        enc.localBinMatrix(runtime, arr2, enc_arr2).get();
        enc.localBinMatrix(runtime, arr3, enc_arr3).get();
        enc.localBinMatrix(runtime, arr4, enc_arr4).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_arr1).get();
        enc.remoteBinMatrix(runtime, enc_arr2).get();
        enc.remoteBinMatrix(runtime, enc_arr3).get();
        enc.remoteBinMatrix(runtime, enc_arr4).get();
    }

    aby3::sbMatrix sort_test;
    odd_even_merge(enc_arr1, enc_arr2, sort_test, role, enc, eval, runtime);

    aby3::i64Matrix test(arr1_len + arr2_len, 1);
    enc.revealAll(runtime, sort_test, test).get();

    aby3::sbMatrix multi_sort_test;
    std::vector<aby3::sbMatrix> enc_data_vec = {enc_arr1, enc_arr2, enc_arr3, enc_arr4};
    odd_even_multi_merge(enc_data_vec, multi_sort_test, role, enc, eval, runtime);

    aby3::i64Matrix multi_test(arr1_len + arr2_len + arr3_len + arr4_len, 1);
    enc.revealAll(runtime, multi_sort_test, multi_test).get();

    if(role == 0){
        check_result("Odd Even Merge Sort Test", test, res_res_);
        check_result("Odd Even Multi Merge Sort Test", multi_test, multi_res_res_);
    }

    std::vector<sbMatrix> enc_data_vec1 = {enc_arr1, enc_arr3};
    std::vector<sbMatrix> enc_data_vec2 = {enc_arr3, enc_arr1};
    std::vector<sbMatrix> high_dim_res(2);
    high_dimensional_odd_even_merge(enc_data_vec1, enc_data_vec2, high_dim_res, role, enc, eval, runtime);

    aby3::i64Matrix high_dim_test1(arr1_len + arr3_len, 1);
    enc.revealAll(runtime, high_dim_res[0], high_dim_test1).get();
    aby3::i64Matrix high_dim_test2(arr1_len + arr3_len, 1);
    enc.revealAll(runtime, high_dim_res[1], high_dim_test2).get();

    std::vector<int> high_dim_res1(arr1_len + arr3_len);
    for(size_t i=0; i<arr1_len; i++) high_dim_res1[i] = arr1(i, 0);
    for(size_t i=0; i<arr3_len; i++) high_dim_res1[i + arr1_len] = arr3(i, 0);
    std::sort(high_dim_res1.begin(), high_dim_res1.end());
    aby3::i64Matrix high_dim_res1_(arr1_len + arr3_len, 1);
    for(size_t i=0; i<arr1_len + arr3_len; i++) high_dim_res1_(i, 0) = high_dim_res1[i];

    if(role == 0){
        check_result("High Dimensional Odd Even Merge Sort Test 1", high_dim_test1, high_dim_res1_);
        check_result("High Dimensional Odd Even Merge Sort Test 2", high_dim_test2, high_dim_res1_);
    }

    std::vector<sbMatrix> multi_enc_data_vec_1 = {enc_arr1, enc_arr2};
    std::vector<sbMatrix> multi_enc_data_vec_2 = {enc_arr3, enc_arr4};
    std::vector<std::vector<sbMatrix>> hd_multi_enc_data_vec(2);
    for(size_t i=0; i<2; i++){
        hd_multi_enc_data_vec[i] = (i == 0) ? multi_enc_data_vec_1 : multi_enc_data_vec_2;
    }

    std::vector<sbMatrix> hd_multi_res(2);
    high_dimensional_odd_even_multi_merge(hd_multi_enc_data_vec, hd_multi_res, role, enc, eval, runtime);

    aby3::i64Matrix hd_multi_test1(arr1_len + arr2_len, 1);
    enc.revealAll(runtime, hd_multi_res[0], hd_multi_test1).get();
    aby3::i64Matrix hd_multi_test2(arr3_len + arr4_len, 1);
    enc.revealAll(runtime, hd_multi_res[1], hd_multi_test2).get();

    std::vector<int> hd_multi_res2(arr3_len + arr4_len);
    for(size_t i=0; i<arr3_len; i++) hd_multi_res2[i] = arr3(i, 0);
    for(size_t i=0; i<arr4_len; i++) hd_multi_res2[i + arr3_len] = arr4(i, 0);
    std::sort(hd_multi_res2.begin(), hd_multi_res2.end());
    aby3::i64Matrix hd_multi_res2_(arr3_len + arr4_len, 1);
    for(size_t i=0; i<arr3_len + arr4_len; i++) hd_multi_res2_(i, 0) = hd_multi_res2[i];

    if(role == 0){
        check_result("High Dimensional Odd Even Multi Merge Sort Test 1", hd_multi_test1, res_res_);
        check_result("High Dimensional Odd Even Multi Merge Sort Test 2", hd_multi_test2, hd_multi_res2_);
    }


    return 0;
}