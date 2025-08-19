#include "utils.h"

using namespace oc;
using namespace aby3;

void prefixsum(int pIdx, si64Matrix &v, si64Matrix &result){

    size_t v_len = v.rows();
    size_t result_len = v_len;

    result.resize(result_len, v.cols());
    i64 sum_0=0,sum_1=0;
    for(size_t i=0; i<v_len; i++){
        sum_0 += v.mShares[0](i, 0);
        sum_1 += v.mShares[1](i, 0);
        result.mShares[0](i, 0) = sum_0;
        result.mShares[1](i, 0) = sum_1;
    }

    return;
}

void prefixsum_inv(int pIdx, si64Matrix &v, si64Matrix &result){
    size_t v_len = v.rows();
    size_t result_len = v_len;

    result.resize(result_len, v.cols());
    for(size_t i=v_len-1; i>=1; i--){
        result.mShares[0](i, 0) = v.mShares[0](i, 0)-v.mShares[0](i-1, 0);
        result.mShares[1](i, 0) = v.mShares[1](i, 0)-v.mShares[1](i-1, 0);
    }
    result.mShares[0](0, 0) = v.mShares[0](0, 0);
    result.mShares[1](0, 0) = v.mShares[1](0, 0);

    return;
}

void permutation_inverse(i64Matrix& rsigma_plain, i64Matrix& rsigma_inv_plain) {
    int n = rsigma_plain.rows();
    rsigma_inv_plain.resize(n, rsigma_plain.cols());
    
    // hash
    std::unordered_map<int, int> value_to_position;
    
    // step0：建立值到位置的映射
    for (int i = 0; i < n; i++) {
        int value = rsigma_plain(i, 0);
        value_to_position[value] = i ;  
    }
    
    // step1：构建逆置换
    for (int i = 0; i < n; i++) {
        int value = rsigma_plain(i, 0);
        rsigma_inv_plain(value , 0) = value_to_position[value];
    }
    
    return ;
}

void fed_argsort(int pIdx, si64Matrix &v, si64Matrix &result, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    size_t v_len = v.rows();
    size_t result_len = v_len;

    std::vector<i64Matrix> value_plain(v_len);
    std::vector<si64Matrix> value_enc(v_len);

    for(size_t i=0; i<v_len; i++){
        value_plain[i].resize(1, 1);
        value_enc[i].resize(1, 1);
        value_plain[i](0,0)=i;
        if(pIdx==0){
            enc.localIntMatrix(runtime, value_plain[i], value_enc[i]).get();
        }
        else{
            enc.remoteIntMatrix(runtime, value_enc[i]).get();
        }
    }

    //DEBUG
    if(pIdx==0){
        std::cout << "value_enc_plain: " << std::endl;
    }
    i64Matrix temp(v_len, 1);
    for(size_t i=0; i<v_len; i++){
        enc.revealAll(runtime, value_enc[i], temp).get();  
        if(pIdx==0){
            std::cout << temp(0, 0) << " ";    
        }
    }
    if(pIdx==0){
        std::cout << std::endl;  
        std::cout.flush(); 
    }
 
        
    
    //----value_enc_plain correct


    //DEBUG
    if(pIdx==0){
        std::cout << "v_key_plain1: " << std::endl;
    }
    i64Matrix temp3(v_len, 1);
    for(size_t i=0; i<v_len; i++){
        enc.revealAll(runtime, v, temp3).get();
        if(pIdx==0){
            std::cout << temp3(i, 0) << " ";  
        }
    }
    if(pIdx==0){
        std::cout << std::endl;  
        std::cout.flush(); 
    }   
    //TODO： minsize
    quick_sort_with_other_elements(v, value_enc, pIdx, enc, eval, runtime,5);
    if(pIdx==0){
        std::cout << "v_key_plain2: " << std::endl;
    }
    for(size_t i=0; i<v_len; i++){
        enc.revealAll(runtime, v, temp3).get();
        if(pIdx==0){
            std::cout << temp3(i, 0) << " ";  
        }
    }
    if(pIdx==0){
        std::cout << std::endl;  
        std::cout.flush(); 
    } 
    //DEBUG
    if(pIdx==0){
        std::cout << "arg_plain: " << std::endl;
    }
    i64Matrix temp2(v_len, 1);
    for(size_t i=0; i<v_len; i++){
        enc.revealAll(runtime, value_enc[i], temp2).get();
        if(pIdx==0){
            std::cout << temp2(0, 0) << " ";  
        }
    }
    if(pIdx==0){
        std::cout << std::endl;  
        std::cout.flush(); 
    }   
    
    //----arg_plain correct

    for(size_t i=0; i<v_len; i++){
        result.mShares[0](i, 0) = value_enc[i].mShares[0](0, 0);
        result.mShares[1](i, 0) = value_enc[i].mShares[1](0, 0);
    }

    return ;
}