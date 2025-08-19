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