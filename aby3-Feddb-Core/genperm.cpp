#include "genperm.h"
#include "utils.h"

using namespace oc;
using namespace aby3;


si64Matrix reshare_matrix(int pIdx, si64Matrix& x, int targetPartyIdx, Sh3Encryptor& enc, Sh3Runtime& runtime){
    
    int prevPartyIdx = (targetPartyIdx + 2) % 3;  // P_{i-1}
    int nextPartyIdx = (targetPartyIdx + 1) % 3;  // P_{i+1}
    
    auto& comm = runtime.mComm;
    
    int n = x.rows();
    int cols = x.cols();
    //generate random shares r_1, r_2, r_3 such that r_1 + r_2 + r_3 = 0
    i64Matrix r(n,cols);   
    std::fill_n(r.data(), n*cols, 0);
    si64Matrix r_share(n,cols);
    // for(size_t i=0;i<n;i++){
    //     for(size_t j=0;j<cols;j++){
    //         r(i,j)=0;
    //     }
    // }
    if (pIdx == prevPartyIdx) {
        enc.localIntMatrix(runtime, r, r_share).get();
    } else {
        enc.remoteIntMatrix(runtime, r_share).get();
    }

    si64Matrix x_re(n, cols);
    //i64Matrix x_reshare1(n, cols*64), x_reshare2(n, cols*64);
    if (pIdx == prevPartyIdx) {
        for(size_t i=0;i<n;i++){
            for(size_t j=0;j<cols;j++){
                x_re.mShares[0](i,j) = x.mShares[0](i,j) + r_share.mShares[0](i,j);
                x_re.mShares[1](i,j) = x.mShares[1](i,j) + r_share.mShares[1](i,j);    
            }
        }
        large_data_sending(pIdx, x_re.mShares[0], runtime, true);
        return x_re;

    } else if (pIdx == nextPartyIdx) {//p2
        for(size_t i=0;i<n;i++){
            for(size_t j=0;j<cols;j++){
                x_re.mShares[0](i,j) = x.mShares[0](i,j) + r_share.mShares[0](i,j);
                x_re.mShares[1](i,j) = x.mShares[1](i,j) + r_share.mShares[1](i,j);
            }
        }
        large_data_sending(pIdx, x_re.mShares[1], runtime, false);
        return x_re;
    }
    else if(pIdx == targetPartyIdx){
        large_data_receiving(pIdx, x_re.mShares[1], runtime, true);
        large_data_receiving(pIdx, x_re.mShares[0], runtime, false);
        return x_re;
    }

}

//取出key的第d位
void getBitKey(int pIdx, sbMatrix &k_bool, si64Matrix &k_j_arith, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    // 保留旧实现（按单个 d 位提取）
    // sbMatrix k_j(k_bool.rows(),64);
    // for (size_t i = 0; i < k_bool.rows(); i++)
    // {
    //     auto word = d / 64;
    //     auto offset = d % 64;
    //
    //     k_j.mShares[0](i, 0) = (k_bool.mShares[0](i, word) >> offset) & 1;
    //     k_j.mShares[1](i, 0) = (k_bool.mShares[1](i, word) >> offset) & 1;
    // }
    // bool2arith(pIdx, k_j, k_j_arith, enc, eval, runtime);
    // return;

    // 新实现：把所有 bit 按 [bit * rows, (bit + 1) * rows) 拼接到一个大 sbMatrix，
    // 再一次性 bool2arith，输出尺寸为 (bitcount * rows, 1)。

    const size_t rows = k_bool.rows();
    const size_t bitcount = k_bool.bitCount();
    sbMatrix k_j(rows * bitcount, 64);

    for (size_t bit = 0; bit < bitcount; ++bit) {
        const size_t word = bit / 64;
        const size_t offset = bit % 64;
        const size_t base = bit * rows;

        for (size_t i = 0; i < rows; ++i) {
            const size_t dst = base + i;
            k_j.mShares[0](dst, 0) = (k_bool.mShares[0](i, word) >> offset) & 1;
            k_j.mShares[1](dst, 0) = (k_bool.mShares[1](i, word) >> offset) & 1;
        }
    }

    k_j_arith.resize(rows * bitcount, 1);
    bool2arith(pIdx, k_j, k_j_arith, enc, eval, runtime);
    return;
}

//k_j只有一位
void genBitPerm(int pIdx, si64Matrix &k_j, si64Matrix &perm, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    int n = k_j.rows();

    si64Matrix f0(n,1),f1(n,1),one(n,1);
    set_const_share(pIdx, 1, one, enc, eval, runtime);

    f1 = k_j;
    f0 = one - k_j;
    si64Matrix s0(n,1),s1(n,1),s_sub(n,1);
    prefixsum(pIdx, f0, s0);
    prefixsum_with_initial_elements(pIdx, f1, s1, s0);
    s_sub = s1 - s0;
    //DEBUG
    // i64Matrix s0_plain(n,1),s1_plain(n,1);
    // enc.revealAll(runtime, s0, s0_plain).get();
    // enc.revealAll(runtime, s1, s1_plain).get();
    // std::cout << "s0: " << std::endl;
    // for(size_t i=0;i<n;i++){
    //     std::cout << s0_plain(i,0) << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "s1: " << std::endl;
    // for(size_t i=0;i<n;i++){
    //     std::cout << s1_plain(i,0) << " ";
    // }
    // std::cout << std::endl;
    //----s0,s1 correct

    // for(size_t i=0;i<n;i++){
    //     s_sub.mShares[0](i,0) = s1.mShares[0](i,0) - s0.mShares[0](i,0);
    //     s_sub.mShares[1](i,0) = s1.mShares[1](i,0) - s0.mShares[1](i,0);
    // }

    si64Matrix t(n,1);
    cipher_mul_seq(pIdx, k_j, s_sub, t, eval, enc, runtime);

    // for(size_t i=0;i<n;i++){
    //     perm.mShares[0](i,0) = s0.mShares[0](i,0) + t.mShares[0](i,0)-one.mShares[0](i,0);
    //     perm.mShares[1](i,0) = s0.mShares[1](i,0) + t.mShares[1](i,0)-one.mShares[1](i,0);
    // }
    perm = s0 + t - one;
    
    return ;
}

void applyPerm(int pIdx, si64Matrix &perm, si64Matrix &k_j, si64Matrix &k_j_prime, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = k_j.rows();
    size_t cols = k_j.cols();

    // generate the prev, next - correlated randomness.
    // 1 - generate the permutations. pi
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    //2- shuffle:pi(perm)=perm_prime
    si64Matrix perm_prime(len,1);
    si64Matrix next_perm(len,1),prev_perm(len,1);
    prev_perm=perm;
    for(size_t i=0;i<3;i++){
        if(pIdx==i){
            permutate(pIdx, prev_perm, next_perm, prev_permutation);
            
        }
        else if(pIdx==((i+2)%3)){
            permutate(pIdx, prev_perm, next_perm, next_permutation);
            
        }
        else if(pIdx==(i+1)%3){
            std::fill_n(next_perm.mShares[0].data(),len,0);
            std::fill_n(next_perm.mShares[1].data(),len,0);
 
        }
        prev_perm=reshare_matrix(pIdx, next_perm, (i+1)%3, enc, runtime);

    }
    perm_prime=prev_perm;


    //3-shuffle:pi(k_j)=k_j_perm
    si64Matrix k_j_perm(len,cols);
    next_perm.resize(len,cols);
    prev_perm.resize(len,cols);
    prev_perm=k_j;
    for(int id=0;id<3;id++){
        if(pIdx==id){
            permutate(pIdx, prev_perm, next_perm, prev_permutation);
            
        }
        else if(pIdx==((id+2)%3)){
            permutate(pIdx, prev_perm, next_perm, next_permutation);
           
        }
        
        else if(pIdx==(id+1)%3){
            // for(size_t i=0;i<len;i++){
            //     for(size_t j=0;j<cols;j++){
            //         next_perm.mShares[0](i, j) = 0;
            //         next_perm.mShares[1](i, j) = 0;
            //     }
            // }
            const size_t total = len * cols;
            std::fill_n(next_perm.mShares[0].data(), total, 0);
            std::fill_n(next_perm.mShares[1].data(), total, 0);
        }
        prev_perm=reshare_matrix(pIdx, next_perm, (id+1)%3, enc, runtime);
        //DEBUG
        // i64Matrix tmp(len,1);
        // enc.revealAll(runtime, prev_perm, tmp).get();
        // std::cout << " prev_perm: " << std::endl;
        // for(size_t i=0;i<len;i++){
        //     std::cout << "pidx: " << pIdx << " " << tmp(i,0) << " ";
        // }
        // std::cout << std::endl;
    }
    k_j_perm=prev_perm;

    //4-revealall perm_prime
    i64Matrix perm_prime_plain(len,1);
    enc.revealAll(runtime, perm_prime, perm_prime_plain).get();
    permutate(pIdx, k_j_perm, k_j_prime, perm_prime_plain);



    return ;
}

void composePerm(int pIdx, si64Matrix &perm_a, si64Matrix &perm_b, si64Matrix &perm, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = perm_a.rows();

    // generate the prev, next - correlated randomness.
    // 1 - generate the permutations. pi
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);
    std::vector<size_t> prev_inverse_permutation;
    std::vector<size_t> next_inverse_permutation;
    get_inverse_permutation(prev_permutation, prev_inverse_permutation);
    get_inverse_permutation(next_permutation, next_inverse_permutation);

    //2 - shuffle:pi(perm_a)=perm_a_prime
    si64Matrix perm_a_prime(len,1);
    si64Matrix next_perm(len,1),prev_perm(len,1);
    prev_perm=perm_a;
    for(int id=0;id<3;id++){
        if(pIdx==id){
            permutate(pIdx, prev_perm, next_perm, prev_permutation);
        }
        else if(pIdx==((id+2)%3)){
            permutate(pIdx, prev_perm, next_perm, next_permutation);
        }
        //DEBUG
        else if(pIdx==(id+1)%3){
            std::fill_n(next_perm.mShares[0].data(),len,0);
            std::fill_n(next_perm.mShares[1].data(),len,0);

        }
        prev_perm=reshare_matrix(pIdx, next_perm, (id+1)%3, enc, runtime);
        
    }
    perm_a_prime=prev_perm;

    //3 - revealall perm_a_prime , get inverse permutation
    i64Matrix perm_a_prime_plain(len,1),perm_a_prime_inv_plain(len,1);
    enc.revealAll(runtime, perm_a_prime, perm_a_prime_plain).get();
    permutation_inverse(perm_a_prime_plain, perm_a_prime_inv_plain);


    //4 - plain_permutate:perm_a_prime_plain(perm_b)=perm_b_prime
    si64Matrix perm_b_prime(len,1);
    permutate(pIdx, perm_b, perm_b_prime, perm_a_prime_inv_plain);


    //5 - unshuffle:pi^{-1}(perm_b_prime)=perm_b
    prev_perm=perm_b_prime;
    for(int id=2;id>=0;id--){
        if(pIdx==id){
            permutate(pIdx, prev_perm, next_perm, prev_inverse_permutation);
            
        }
        else if(pIdx==((id+2)%3)){
            permutate(pIdx, prev_perm, next_perm, next_inverse_permutation);
            
        }
        //DEBUG
        else if(pIdx==(id+1)%3){
            std::fill_n(next_perm.mShares[0].data(),len,0);
            std::fill_n(next_perm.mShares[1].data(),len,0);

        }
        prev_perm=reshare_matrix(pIdx, next_perm, (id+1)%3, enc, runtime);
        
    }
    perm=prev_perm;

    return ;
}

void genPerm(int pIdx, si64Matrix &k, si64Matrix &perm, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    int n = k.rows();
    int cols = k.cols();
    int bitcount = 64*cols;

    sbMatrix k_bool(n,bitcount);
    arith2bool(pIdx, k, k_bool, enc, eval, runtime);
    si64Matrix k_j_arith_vector(n*bitcount, 1);
    getBitKey(pIdx, k_bool, k_j_arith_vector, enc, eval, runtime);

    si64Matrix k_0(n,1);
    k_0.mShares[0] = k_j_arith_vector.mShares[0].block(0, 0, n, 1);
    k_0.mShares[1] = k_j_arith_vector.mShares[1].block(0, 0, n, 1);
    //getBitKey(pIdx, k_bool, 0, k_0, enc, eval, runtime);
    
    si64Matrix pre_perm(n,1);
    genBitPerm(pIdx, k_0, pre_perm, enc, eval, runtime);

    
    for(int d=1;d<bitcount;d++){
        si64Matrix k_j(n,1),k_j_prime(n,1);

        k_j.mShares[0] = k_j_arith_vector.mShares[0].block(d*n, 0, n, 1);
        k_j.mShares[1] = k_j_arith_vector.mShares[1].block(d*n, 0, n, 1);

        //DEBUG
        // i64Matrix k_j_plain(n,1);
        // enc.revealAll(runtime, k_j, k_j_plain).get();
        // std::cout << "k_j: " << std::endl;
        // for(size_t i=0;i<n;i++){
        //     std::cout << k_j_plain(i,0) << " ";
        // }
        // std::cout << std::endl;
        //----k_j correct
        applyPerm(pIdx, pre_perm, k_j,k_j_prime, enc, eval, runtime);
        //DEBUG
        // i64Matrix k_j_prime_plain(n,1);
        // enc.revealAll(runtime, k_j_prime, k_j_prime_plain).get();
        // std::cout <<"d: " << d << " k_j_prime: " << std::endl;
        // for(size_t i=0;i<n;i++){
        //     std::cout << k_j_prime_plain(i,0) << " ";
        // }
        // std::cout << std::endl;
        //----k_j_prime correct

        si64Matrix next_perm(n,1);
        genBitPerm(pIdx, k_j_prime, next_perm, enc, eval, runtime);
        //DEBUG
        // i64Matrix next_perm_plain(n,1);
        // enc.revealAll(runtime, next_perm, next_perm_plain).get();
        // std::cout << "next_perm: " << std::endl;
        // for(size_t i=0;i<n;i++){
        //     std::cout << next_perm_plain(i,0) << " ";
        // }
        // std::cout << std::endl;
        //----next_perm correct
        composePerm(pIdx, pre_perm, next_perm, perm, enc, eval, runtime);
        //DEBUG
        // i64Matrix perm_plain(n,1);
        // enc.revealAll(runtime, perm, perm_plain).get();
        // std::cout << "perm: " << std::endl;
        // for(size_t i=0;i<n;i++){
        //     std::cout << perm_plain(i,0) << " ";
        // }
        // std::cout << std::endl;
        pre_perm=perm;
    }
    
    return;
}

void genPerm_bool(int pIdx, sbMatrix &k, si64Matrix &perm, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    int n = k.rows();
    int bitcount = k.bitCount();

    si64Matrix k_0(n,1);
    //getBitKey(pIdx, k, 0, k_0, enc, eval, runtime);
    si64Matrix k_j_arith_vector(n*bitcount, 1);
    getBitKey(pIdx, k, k_j_arith_vector, enc, eval, runtime);
    k_0.mShares[0] = k_j_arith_vector.mShares[0].block(0, 0, n, 1);
    k_0.mShares[1] = k_j_arith_vector.mShares[1].block(0, 0, n, 1);
    
    si64Matrix pre_perm(n,1);
    genBitPerm(pIdx, k_0, pre_perm, enc, eval, runtime);

    
    for(int d=1;d<bitcount;d++){
        si64Matrix k_j(n,1),k_j_prime(n,1);
        //getBitKey(pIdx, k, d, k_j, enc, eval, runtime);
        k_j.mShares[0] = k_j_arith_vector.mShares[0].block(d*n, 0, n, 1);
        k_j.mShares[1] = k_j_arith_vector.mShares[1].block(d*n, 0, n, 1);
        applyPerm(pIdx, pre_perm, k_j,k_j_prime, enc, eval, runtime);


        si64Matrix next_perm(n,1);
        genBitPerm(pIdx, k_j_prime, next_perm, enc, eval, runtime);
        composePerm(pIdx, pre_perm, next_perm, perm, enc, eval, runtime);

        pre_perm=perm;
    }
    
    return;
}

void concat_k_v(int pIdx, sbMatrix& k, sbMatrix& v, sbMatrix& res){

    //assume k、v less than 32bit
    for(size_t i=0;i<k.rows();i++){
        res.mShares[0](i, 0) = (k.mShares[0](i, 0) << 32) | (v.mShares[0](i, 0)&0xFFFFFFFF);
        res.mShares[1](i, 0) = (k.mShares[1](i, 0) << 32) | (v.mShares[1](i, 0)&0xFFFFFFFF);
    }
    return;
}

// void genPerm_kv(int pIdx, si64Matrix &k,si64Matrix &v, si64Matrix &perm, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
//     int n = k.rows();

//     sbMatrix k_bool(n,64);
//     arith2bool(pIdx, k, k_bool, enc, eval, runtime);
//     sbMatrix v_bool(n,64);
//     arith2bool(pIdx, v, v_bool, enc, eval, runtime);

//     //k和v拼接
//     sbMatrix key_concat(n,64);
//     concat_k_v(pIdx, k_bool, v_bool, key_concat);

//     si64Matrix k_0(n,1);
//     getBitKey(pIdx, key_concat, 0, k_0, enc, eval, runtime);
    
//     si64Matrix pre_perm(n,1);
//     genBitPerm(pIdx, k_0, pre_perm, enc, eval, runtime);

//     for(int d=1;d<64;d++){
//         si64Matrix k_j(n,1),k_j_prime(n,1);
//         getBitKey(pIdx, key_concat, d, k_j, enc, eval, runtime);
//         //DEBUG
//         // i64Matrix k_j_plain(n,1);
//         // enc.revealAll(runtime, k_j, k_j_plain).get();
//         // std::cout << "k_j: " << std::endl;
//         // for(size_t i=0;i<n;i++){
//         //     std::cout << k_j_plain(i,0) << " ";
//         // }
//         // std::cout << std::endl;
//         //----k_j correct
//         applyPerm(pIdx, pre_perm, k_j,k_j_prime, enc, eval, runtime);
//         //DEBUG
//         // i64Matrix k_j_prime_plain(n,1);
//         // enc.revealAll(runtime, k_j_prime, k_j_prime_plain).get();
//         // std::cout <<"d: " << d << " k_j_prime: " << std::endl;
//         // for(size_t i=0;i<n;i++){
//         //     std::cout << k_j_prime_plain(i,0) << " ";
//         // }
//         // std::cout << std::endl;
//         //----k_j_prime correct

//         si64Matrix next_perm(n,1);
//         genBitPerm(pIdx, k_j_prime, next_perm, enc, eval, runtime);
//         //DEBUG
//         // i64Matrix next_perm_plain(n,1);
//         // enc.revealAll(runtime, next_perm, next_perm_plain).get();
//         // std::cout << "next_perm: " << std::endl;
//         // for(size_t i=0;i<n;i++){
//         //     std::cout << next_perm_plain(i,0) << " ";
//         // }
//         // std::cout << std::endl;
//         //----next_perm correct
//         composePerm(pIdx, pre_perm, next_perm, perm, enc, eval, runtime);
//         //DEBUG
//         // i64Matrix perm_plain(n,1);
//         // enc.revealAll(runtime, perm, perm_plain).get();
//         // std::cout << "perm: " << std::endl;
//         // for(size_t i=0;i<n;i++){
//         //     std::cout << perm_plain(i,0) << " ";
//         // }
//         // std::cout << std::endl;
//         pre_perm=perm;
//     }
    
//     return;
// }





