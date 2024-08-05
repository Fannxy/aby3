#include "Search.h"

using namespace oc;
using namespace aby3;


int mcompBS(std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    // get the length.
    int n = keyset.size();
    if(n == 1) {
        boolShare le_res(1, pIdx);
        res.mShares[0](0, 0) = le_res.bshares[0];
        res.mShares[1](0, 0) = le_res.bshares[1];
        return 0;
    }
    int bitsize = keyset[0].bitCount();

    // comp key with keyset - can be compute in parallel.
    aby3::sbMatrix expand_key_set(n, bitsize);
    aby3::sbMatrix expand_key(n, bitsize);
    for(int i = 0; i < n; i++){
        expand_key_set.mShares[0](i, 0) = keyset[i].mShares[0](0, 0);
        expand_key_set.mShares[1](i, 0) = keyset[i].mShares[1](0, 0);
        expand_key.mShares[0](i, 0) = key.mShares[0](0, 0);
        expand_key.mShares[1](i, 0) = key.mShares[1](0, 0);
    }

    // compute comparison.
    aby3::sbMatrix le_res(n, 1);
    bool_cipher_lt(pIdx, expand_key_set, expand_key, le_res, enc, eval, runtime);
    bool_cipher_not(pIdx, le_res, le_res);

    // differential subtraction.
    res.resize(n, 1);
    res.mShares[0](0, 0) = le_res.mShares[0](0, 0);
    res.mShares[1](0, 0) = le_res.mShares[1](0, 0);

    aby3::sbMatrix diff_res(n-1, 1), left_le(n-1, 1), right_le(n-1, 1);

    for(int i = 0; i < n-1; i++){
        left_le.mShares[0](i, 0) = le_res.mShares[0](i, 0);
        left_le.mShares[1](i, 0) = le_res.mShares[1](i, 0);
        right_le.mShares[0](i, 0) = le_res.mShares[0](i+1, 0);
        right_le.mShares[1](i, 0) = le_res.mShares[1](i+1, 0);
    }

    bool_cipher_not(pIdx, left_le, left_le);
    bool_cipher_and(pIdx, left_le, right_le, diff_res, enc, eval, runtime);

    for(int i = 0; i < n-1; i++){
        res.mShares[0](i+1, 0) = diff_res.mShares[0](i, 0);
        res.mShares[1](i+1, 0) = diff_res.mShares[1](i, 0);
    }

    return 0;
}

int compBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    int n = data.size();
    aby3::sbMatrix comp_flag(n, 1);
    mcompBS(keyset, key, comp_flag, pIdx, enc, eval, runtime);

    // get the result.
    aby3::si64Matrix data_mat(n, 1);
    for(int i = 0; i < n; i++){
        data_mat.mShares[0](i, 0) = data[i].mShares[0](0, 0);
        data_mat.mShares[1](i, 0) = data[i].mShares[1](0, 0);
    }
    // compute the result.
    res.resize(1, 1);
    aby3::si64Matrix res_tmp(n, 1);
    // eval.asyncMul(runtime, data_mat, comp_flag, res_tmp).get();
    arith_bool_mul(pIdx, data_mat, comp_flag, res_tmp, enc, eval, runtime);

    res.mShares[0](0, 0) = res_tmp.mShares[0](0, 0);
    res.mShares[1](0, 0) = res_tmp.mShares[1](0, 0);
    for(int i=1; i<n; i++){
        res.mShares[0](0, 0) += res_tmp.mShares[0](i, 0);
        res.mShares[1](0, 0) += res_tmp.mShares[1](i, 0);
    }
    return 0;
}

int mtagBS(std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    int n = keyset.size();
    int bitsize = keyset[0].bitCount();

    if(checkPowerOfTwo(n) == false){
        THROW_RUNTIME_ERROR("Currently, the size of data should be power of 2!");
    }
    int logn = log2(n);

    // get the mid point.
    int p = n/2 - 1;
    sbMatrix c(1, 1);
    boolShare q(1, pIdx);

    bool_cipher_lt(pIdx, keyset[p], key, c, enc, eval, runtime);
    bool_cipher_not(pIdx, c, c);

    aby3::sbMatrix tag(1, 1);
    aby3::sbMatrix tag_next;
    tag.mShares[0](0, 0) = q.bshares[0]; tag.mShares[1](0, 0) = q.bshares[1];

    // log-round binary search.
    for(size_t i=0; i<logn; i++){
        int tag_size = (1 << i);

        aby3::sbMatrix keymat(tag_size, bitsize);
        aby3::sbMatrix keytag(tag_size, bitsize);
        aby3::sbMatrix res_key_tmp(tag_size, bitsize);
        for(size_t j=0; j<tag_size; j++){
            int tag_ind = (((2*(j+1) - 1) * n) / (2*tag_size)) - 1;
            keymat.mShares[0](j, 0) = keyset[tag_ind].mShares[0](0, 0);
            keymat.mShares[1](j, 0) = keyset[tag_ind].mShares[1](0, 0);
            keytag.mShares[0](j, 0) = (tag.mShares[0](j, 0) == 1) ? -1 : 0;
            keytag.mShares[1](j, 0) = (tag.mShares[1](j, 0) == 1) ? -1 : 0;
        }
        // bool_cipher_and(pIdx, keymat, keytag, res_key_tmp, enc, eval, runtime);

        aby3::sbMatrix akey(1, bitsize);
        // akey.mShares[0](0, 0) = res_key_tmp.mShares[0](0, 0);
        // akey.mShares[1](0, 0) = res_key_tmp.mShares[1](0, 0);
        // for(size_t j=1; j<tag_size; j++){
        //     akey.mShares[0](0, 0) ^= res_key_tmp.mShares[0](j, 0);
        //     akey.mShares[1](0, 0) ^= res_key_tmp.mShares[1](j, 0);
        // }

        std::vector<sbMatrix> enc_key_mat = {keymat};
        std::vector<sbMatrix> enc_key_tag = {keytag};
        constant_bool_dot(pIdx, enc_key_mat, enc_key_tag, akey, enc, eval, runtime);

        bool_cipher_lt(pIdx, akey, key, c, enc, eval, runtime);
        bool_cipher_not(pIdx, c, c);

        aby3::sbMatrix notC(1, 1);
        bool_cipher_not(pIdx, c, notC);

        // compute the tag for the next layer.
        int next_size = (1 << (i+1));
        tag_next.resize(next_size, 1);
        aby3::sbMatrix tag_expand(next_size, 1);
        aby3::sbMatrix comp_and_inv_comp(next_size, 1);

        for(size_t j=0; j<next_size; j++){
            tag_expand.mShares[0](j, 0) = tag.mShares[0](j/2, 0);
            tag_expand.mShares[1](j, 0) = tag.mShares[1](j/2, 0);
            if(j%2 == 0){
                comp_and_inv_comp.mShares[0](j, 0) = c.mShares[0](0, 0);
                comp_and_inv_comp.mShares[1](j, 0) = c.mShares[1](0, 0);
            }
            else{
                comp_and_inv_comp.mShares[0](j, 0) = notC.mShares[0](0, 0);
                comp_and_inv_comp.mShares[1](j, 0) = notC.mShares[1](0, 0);
            }
        }

        bool_cipher_and(pIdx, tag_expand, comp_and_inv_comp, tag_next, enc, eval, runtime);
        tag = tag_next;
    }
    
    res = tag;

    return 0;
}

int tagBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    int n = data.size();
    int bitsize = keyset[0].bitCount();

    if(checkPowerOfTwo(n) == false){
        THROW_RUNTIME_ERROR("Currently, the size of data should be power of 2!");
    }
    int logn = log2(n);
    // get the mid point.
    int p = n/2 - 1;
    sbMatrix c(1, 1);
    boolShare q(1, pIdx);
    bool_cipher_lt(pIdx, keyset[p], key, c, enc, eval, runtime);
    bool_cipher_not(pIdx, c, c);

    si64Matrix pivot(1, 1);
    pivot.mShares[0](0, 0) = data[p].mShares[0](0, 0) - data[n-1].mShares[0](0, 0);
    pivot.mShares[1](0, 0) = data[p].mShares[1](0, 0) - data[n-1].mShares[1](0, 0);
    si64Matrix _pivot(1, 1);
    arith_bool_mul(pIdx, pivot, c, _pivot, enc, eval, runtime);
    pivot.mShares[0](0, 0) = _pivot.mShares[0](0, 0) + data[n-1].mShares[0](0, 0);
    pivot.mShares[1](0, 0) = _pivot.mShares[1](0, 0) + data[n-1].mShares[1](0, 0);

    aby3::sbMatrix prev_tag(1, 1);
    prev_tag.mShares[0](0, 0) = q.bshares[0];
    prev_tag.mShares[1](0, 0) = q.bshares[1];

    // log-round binary search.
    for(size_t i=1; i<logn; i++){
        int tag_size = (1 << i);
        aby3::sbMatrix tag(tag_size, 1);
        aby3::sbMatrix comp_inv(1, 1); bool_cipher_not(pIdx, c, comp_inv);

        aby3::sbMatrix prev_tag_expand(tag_size, 1);
        aby3::sbMatrix prev_comp_and_inv_comp(tag_size, 1);
    
        for(size_t j=0; j<tag_size; j++){
            prev_tag_expand.mShares[0](j, 0) = prev_tag.mShares[0](j/2, 0);
            prev_tag_expand.mShares[1](j, 0) = prev_tag.mShares[1](j/2, 0);
            if(j%2 == 0){
                prev_comp_and_inv_comp.mShares[0](j, 0) = c.mShares[0](0, 0);
                prev_comp_and_inv_comp.mShares[1](j, 0) = c.mShares[1](0, 0);
            }
            else{
                prev_comp_and_inv_comp.mShares[0](j, 0) = comp_inv.mShares[0](0, 0);
                prev_comp_and_inv_comp.mShares[1](j, 0) = comp_inv.mShares[1](0, 0);
            }
        }

        // compute the tag.
        bool_cipher_and(pIdx, prev_tag_expand, prev_comp_and_inv_comp, tag, enc, eval, runtime);

        prev_tag = tag;

        // aggregate this level.
        aby3::si64Matrix data_mat(tag_size, 1);
        aby3::sbMatrix key_mat(tag_size, bitsize);
        aby3::sbMatrix keytag(tag_size, bitsize);
        for(size_t j=0; j<tag_size; j++){
            int tag_ind = (((2*(j+1) - 1) * n) / (2*tag_size)) - 1;
            data_mat.mShares[0](j, 0) = data[tag_ind].mShares[0](0, 0);
            data_mat.mShares[1](j, 0) = data[tag_ind].mShares[1](0, 0);
            key_mat.mShares[0](j, 0) = keyset[tag_ind].mShares[0](0, 0);
            key_mat.mShares[1](j, 0) = keyset[tag_ind].mShares[1](0, 0);
            keytag.mShares[0](j, 0) = (tag.mShares[0](j, 0) == 1) ? -1 : 0;
            keytag.mShares[1](j, 0) = (tag.mShares[1](j, 0) == 1) ? -1 : 0;
        }

        // compute the result.
        aby3::si64Matrix res_tmp(tag_size, 1);
        aby3::sbMatrix res_key_tmp(tag_size, bitsize);
        arith_bool_mul(pIdx, data_mat, tag, res_tmp, enc, eval, runtime);
        // bool_cipher_and(pIdx, key_mat, keytag, res_key_tmp, enc, eval, runtime);
        
        aby3::si64Matrix a(1, 1);
        aby3::sbMatrix akey(1, bitsize);
        a.mShares[0](0, 0) = res_tmp.mShares[0](0, 0);
        a.mShares[1](0, 0) = res_tmp.mShares[1](0, 0);
        // akey.mShares[0](0, 0) = res_key_tmp.mShares[0](0, 0);
        // akey.mShares[1](0, 0) = res_key_tmp.mShares[1](0, 0);
        for(size_t j=1; j<tag_size; j++){
            a.mShares[0](0, 0) += res_tmp.mShares[0](j, 0);
            a.mShares[1](0, 0) += res_tmp.mShares[1](j, 0);
            // akey.mShares[0](0, 0) ^= res_key_tmp.mShares[0](j, 0);
            // akey.mShares[1](0, 0) ^= res_key_tmp.mShares[1](j, 0);
        }

        std::vector<sbMatrix> enc_key_mat = {key_mat};
        std::vector<sbMatrix> enc_key_tag = {keytag};
        constant_bool_dot(pIdx, enc_key_mat, enc_key_tag, akey, enc, eval, runtime);

        // update the pivot and the comp to the next level.
        bool_cipher_lt(pIdx, akey, key, c, enc, eval, runtime);
        bool_cipher_not(pIdx, c, c);

        aby3::si64Matrix pivot_tmp(1, 1);
        pivot_tmp.mShares[0](0, 0) = a.mShares[0](0, 0) - pivot.mShares[0](0, 0);
        pivot_tmp.mShares[1](0, 0) = a.mShares[1](0, 0) - pivot.mShares[1](0, 0);
        arith_bool_mul(pIdx, pivot_tmp, c, a, enc, eval, runtime);
        pivot.mShares[0](0, 0) += a.mShares[0](0, 0);
        pivot.mShares[1](0, 0) += a.mShares[1](0, 0);
    }

    // get the result.
    res = pivot;

    return 0;
}

int subHBS(std::vector<aby3::si64Matrix> &data, std::vector<aby3::sbMatrix> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, int alpha, int threshold){
    int n = data.size();
    int bitsize = keyset[0].bitCount();
    int beta = n / alpha;
    if(alpha * beta != n) THROW_RUNTIME_ERROR("The size of data should be alpha * beta!");

    if(n < threshold){
        return compBS(data, keyset, key, res, pIdx, enc, eval, runtime);
    }

    // construct the upper-level tree.
    std::vector<sbMatrix> keyset_upper;
    for(size_t i=0; i<alpha; i++){
        keyset_upper.push_back(keyset[i*beta + (beta-1)]);
    }

    // call MBS.
    aby3::sbMatrix res_upper(alpha, 1); // alpha * 1
    mcompBS(keyset_upper, key, res_upper, pIdx, enc, eval, runtime);
    // b2a
    aby3::si64Matrix all_ones(alpha, 1);
    init_ones(pIdx, enc, runtime, all_ones, alpha);
    aby3::si64Matrix res_upper_a(alpha, 1);
    // eval.asyncMul(runtime, all_ones, res_upper, res_upper_a).get();
    arith_bool_mul(pIdx, all_ones, res_upper, res_upper_a, enc, eval, runtime);

    // constant inner product for data.
    std::vector<si64Matrix> data_lowers(beta);
    std::vector<si64Matrix> tag_lowers(beta);
    si64Matrix data_res(beta, 1);
    for(size_t i=0; i<beta; i++){
        tag_lowers[i] = res_upper_a;
        data_lowers[i].resize(alpha, 1);
        for(size_t j=0; j<alpha; j++){
            data_lowers[i].mShares[0](j, 0) = data[j*beta + i].mShares[0](0, 0);
            data_lowers[i].mShares[1](j, 0) = data[j*beta + i].mShares[1](0, 0);
        }
    }
    constant_sint_dot(pIdx, data_lowers, tag_lowers, data_res, enc, eval, runtime);


    // constant inner product for key.
    std::vector<sbMatrix> key_lowers(beta);
    std::vector<sbMatrix> keytag_lowers(beta);
    sbMatrix key_res(beta, bitsize);
    for(size_t i=0; i<beta; i++){
        // keytag_lowers[i] = res_upper;
        key_lowers[i].resize(alpha, bitsize);
        keytag_lowers[i].resize(alpha, bitsize);
        for(size_t j=0; j<alpha; j++){
            key_lowers[i].mShares[0](j, 0) = keyset[j*beta + i].mShares[0](0, 0);
            key_lowers[i].mShares[1](j, 0) = keyset[j*beta + i].mShares[1](0, 0);
            keytag_lowers[i].mShares[0](j, 0) = (res_upper.mShares[0](j, 0) == 1) ? -1 : 0;
            keytag_lowers[i].mShares[1](j, 0) = (res_upper.mShares[1](j, 0) == 1) ? -1 : 0;
        }
    }
    constant_bool_dot(pIdx, key_lowers, keytag_lowers, key_res, enc, eval, runtime);

    // call BS.
    std::vector<si64Matrix> data_to_be_search(beta);
    std::vector<sbMatrix> keyset_to_be_search(beta);

    for(size_t i=0; i<beta; i++){
        data_to_be_search[i].resize(1, 1);
        keyset_to_be_search[i].resize(1, bitsize);
        data_to_be_search[i].mShares[0](0, 0) = data_res.mShares[0](i, 0);
        data_to_be_search[i].mShares[1](0, 0) = data_res.mShares[1](i, 0);
        keyset_to_be_search[i].mShares[0](0, 0) = key_res.mShares[0](i, 0);
        keyset_to_be_search[i].mShares[1](0, 0) = key_res.mShares[1](i, 0);
    }

    // tagBS(data_to_be_search, keyset_to_be_search, key, res, pIdx, enc, eval, runtime);
    int sub_alpha = sqrtToPowerOfTwo(beta);
    subHBS(data_to_be_search, keyset_to_be_search, key, res, pIdx, enc, eval, runtime, sub_alpha);

    return 0;
}