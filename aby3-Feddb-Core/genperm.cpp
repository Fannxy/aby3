#include "genperm.h"
#include "utils.h"

using namespace oc;
using namespace aby3;

// This function implements the reshare procedure where party P_{i-1} and P_{i+1} 
// reshare [[x]] to party P_i by adding random secret sharing of 0
void reshare(int pIdx, si64& x, int targetPartyIdx, Sh3Encryptor& enc, Sh3Runtime& runtime) {
    
    // Calculate the indices of the three parties involved
    int prevPartyIdx = (targetPartyIdx + 2) % 3;  // P_{i-1}
    int nextPartyIdx = (targetPartyIdx + 1) % 3;  // P_{i+1}
    
    auto& comm = runtime.mComm;
    i64Matrix r(1,1);
    si64Matrix r_share(1,1);
    r(0,0)=0;
    
    // Step 1: Generate random shares r_1, r_2, r_3 such that r_1 + r_2 + r_3 = 0
    if (pIdx == 0) {
        enc.localIntMatrix(runtime, r, r_share).get();
    } else {
        enc.remoteIntMatrix(runtime, r_share).get();
    }
    
    if (pIdx == prevPartyIdx) {    
        //DEBUG
        x.mData[0] = x.mData[0] + r_share.mShares[0](0,0);
        x.mData[1] = x.mData[1] + r_share.mShares[1](0,0);
        std::cout << "pidx: " << pIdx << " x: " << x.mData[0] << " " << x.mData[1] << std::endl;

        comm.mNext.asyncSendCopy(&x.mData[0], sizeof(i64));
        //std::cout << "sended " << x.mData[0] << " to p1"  << std::endl;

    } else if (pIdx == nextPartyIdx) {//p2
        //DEBUG
        //std::cout << "pIdx: " << pIdx << std::endl;

        x.mData[0] = x.mData[0] + r_share.mShares[0](0,0);
        x.mData[1] = x.mData[1] + r_share.mShares[1](0,0);
        std::cout << "pidx: " << pIdx << " x: " << x.mData[0] << " " << x.mData[1] << std::endl;

        comm.mPrev.asyncSendCopy(&x.mData[1], sizeof(i64));
        // std::cout << "DEBUG1: sizeof(x.mData[1]) = " << sizeof(x.mData[1]) << std::endl;
        // std::cout << "sended " << x.mData[1] << " to p1"  << std::endl;

    }
    else if(pIdx == targetPartyIdx){
        // std::cout << "pIdx: " << pIdx << std::endl;
        // std::cout << "recving x.mData[1] from p0"  << std::endl;
        si64 x_prime;
        comm.mPrev.recv(&x_prime.mData[0], sizeof(i64));  // 从 prevPartyIdx 接收
        // std::cout << "recved x.mData[0] from prevPartyIdx: " << x_prime.mData[0] << std::endl;
        // std::cout << "recving x.mData[1] from nextPartyIdx"  << std::endl;
        comm.mNext.recv(&x_prime.mData[1], sizeof(i64));  // 从 nextPartyIdx 接收
        // std::cout << "recved x.mData[1] from nextPartyIdx: " << x_prime.mData[1] << std::endl;
        x.mData[1] = x_prime.mData[0];  // 来自 prevPartyIdx 的数据
        x.mData[0] = x_prime.mData[1];  // 来自 nextPartyIdx 的数据
        std::cout << "pidx: " << pIdx << " x: " << x.mData[0] << " " << x.mData[1] << std::endl;
    }

    i64 x_plain;
    enc.revealAll(runtime, x, x_plain).get();
    std::cout << "pidx: " << pIdx << " x: " << x_plain << std::endl;

    return;
}

void reshare_matrix(int pIdx, si64Matrix& x, int targetPartyIdx, Sh3Encryptor& enc, Sh3Runtime& runtime){
    int n = x.rows();
    //todo:1->n
    for(size_t i=0;i<n;i++){
        si64 x_i;
        x_i.mData[0] = x.mShares[0](i,0);
        x_i.mData[1] = x.mShares[1](i,0);

        i64 prev_x_plain;
        enc.revealAll(runtime, x_i, prev_x_plain).get();
        std::cout << "pidx: " << pIdx << " prev_x(i,0): " << prev_x_plain << std::endl;

        reshare(pIdx, x_i, targetPartyIdx, enc, runtime);
        x.mShares[0](i,0) = x_i.mData[0];
        x.mShares[1](i,0) = x_i.mData[1];
        i64 x_plain;
        enc.revealAll(runtime, x(i,0), x_plain).get();
        std::cout << "pidx: " << pIdx << " x(i,0): " << x_plain << std::endl;
    }
    
    return;
}


void getBitKey(int pIdx, sbMatrix &k_bool, int d, si64Matrix &k_j_arith, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    sbMatrix k_j(k_bool.rows(),64);
    
    for (size_t i = 0; i < k_bool.rows(); i++)
    {
        auto word = d / 64;
        auto offset = d % 64;

        k_j.mShares[0](i, word) = (k_bool.mShares[0](i, word) >> offset) & 1;
        k_j.mShares[1](i, word) = (k_bool.mShares[1](i, word) >> offset) & 1;
    }

    bool2arith(pIdx, k_j, k_j_arith, enc, eval, runtime);
    return;
}

void genBitPerm(int pIdx, si64Matrix &k_j, si64Matrix &perm, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    int n = k_j.rows();

    si64Matrix f0(n,1),f1(n,1),one(n,1);
    for(size_t i=0;i<n;i++){
        //one=1
       if(pIdx==0){
        one.mShares[0](i,0) = 1;
        one.mShares[1](i,0) = 0;
       }
       else if(pIdx==1){
        one.mShares[0](i,0) = 0;
        one.mShares[1](i,0) = 1;
       }
       else{
        one.mShares[0](i,0) = 0;
        one.mShares[1](i,0) = 0;
       }

       //f1=k_j
       f1.mShares[0](i,0) = k_j.mShares[0](i,0);
       f1.mShares[1](i,0) = k_j.mShares[1](i,0);

       //f0=1-k_j
       f0.mShares[0](i,0) = one.mShares[0](i,0) - k_j.mShares[0](i,0);
       f0.mShares[1](i,0) = one.mShares[1](i,0) - k_j.mShares[1](i,0);
    }
    si64Matrix s0(n,1),s1(n,1),s_sub(n,1);
    prefixsum(pIdx, f0, s0);
    prefixsum_with_initial_elements(pIdx, f1, s1, s0);
    //DEBUG
    i64Matrix s0_plain(n,1),s1_plain(n,1);
    enc.revealAll(runtime, s0, s0_plain).get();
    enc.revealAll(runtime, s1, s1_plain).get();
    std::cout << "s0: " << std::endl;
    for(size_t i=0;i<n;i++){
        std::cout << s0_plain(i,0) << " ";
    }
    std::cout << std::endl;
    std::cout << "s1: " << std::endl;
    for(size_t i=0;i<n;i++){
        std::cout << s1_plain(i,0) << " ";
    }
    std::cout << std::endl;
    //----s0,s1 correct

    for(size_t i=0;i<n;i++){
        s_sub.mShares[0](i,0) = s1.mShares[0](i,0) - s0.mShares[0](i,0);
        s_sub.mShares[1](i,0) = s1.mShares[1](i,0) - s0.mShares[1](i,0);
    }

    si64Matrix t(n,1);
    cipher_mul_seq(pIdx, k_j, s_sub, t, eval, enc, runtime);

    for(size_t i=0;i<n;i++){
        perm.mShares[0](i,0) = s0.mShares[0](i,0) + t.mShares[0](i,0)-one.mShares[0](i,0);
        perm.mShares[1](i,0) = s0.mShares[1](i,0) + t.mShares[1](i,0)-one.mShares[1](i,0);
    }
    
    return ;
}

void applyPerm(int pIdx, si64Matrix &perm, si64Matrix &k_j, si64Matrix &k_j_prime, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = k_j.rows();

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
    //todo:1->3
    for(size_t i=0;i<3;i++){
        if(pIdx==i){
            permutate(pIdx, prev_perm, next_perm, prev_permutation);
            //DEBUG
            std::cout << "pidx: " << pIdx << " prev_permutation: " << std::endl;
            for(size_t i=0;i<len;i++){
                std::cout << prev_permutation[i] << " ";
            }
            std::cout << std::endl;
        }
        else if(pIdx==((i+2)%3)){
            permutate(pIdx, prev_perm, next_perm, next_permutation);
            //DEBUG
            std::cout << "i: " << i << std::endl;
            std::cout << "pidx: " << pIdx << " next_permutation: " << std::endl;
            for(size_t i=0;i<len;i++){
                std::cout << next_permutation[i] << " ";
            }
            std::cout << std::endl;
        }
        //DEBUG
        else if(pIdx==(i+1)%3){
            for(size_t i=0;i<len;i++){
                std::cout << "pidx: " << pIdx << " next_perm(i,0): " << next_perm.mShares[0](i,0) << " " << next_perm.mShares[1](i,0) << std::endl;
            }
        }
        reshare_matrix(pIdx, next_perm, (i+1)%3, enc, runtime);
        prev_perm=next_perm;
        //DEBUG
        i64Matrix tmp(len,1);
        enc.revealAll(runtime, prev_perm, tmp).get();
        std::cout << " prev_perm: " << std::endl;
        for(size_t i=0;i<len;i++){
            std::cout << "pidx: " << pIdx << " " << tmp(i,0) << " ";
        }
        std::cout << std::endl;
    }
    perm_prime=prev_perm;
    
    std::cout << "finish shuffle: " <<pIdx << std::endl;

    //3-shuffle:pi(k_j)=k_j_perm
    si64Matrix k_j_perm(len,1);
    prev_perm=k_j;
    for(size_t i=0;i<3;i++){
        if(pIdx==i){
            permutate(pIdx, prev_perm, next_perm, prev_permutation);
        }
        else if(pIdx==(i+2)%3){
            permutate(pIdx, prev_perm, next_perm, next_permutation);
        }
        //reshare_matrix(pIdx, next_perm, (i+1)%3, enc, runtime);
        prev_perm=next_perm;
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
    for(size_t i=0;i<3;i++){
        if(pIdx==i){
            permutate(pIdx, prev_perm, next_perm, prev_permutation);
        }
        else if(pIdx==(i+2)%3){
            permutate(pIdx, prev_perm, next_perm, next_permutation);
        }
        reshare_matrix(pIdx, next_perm, (i+1)%3, enc, runtime);
        prev_perm=next_perm;
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
    for(size_t i=2;i>=0;i--){
        if(pIdx==i){
            permutate(pIdx, prev_perm, next_perm, prev_inverse_permutation);
        }
        else if(pIdx==(i+2)%3){
            permutate(pIdx, prev_perm, next_perm, next_inverse_permutation);
        }
        reshare_matrix(pIdx, next_perm, (i+1)%3, enc, runtime);
        prev_perm=next_perm;
    }
    perm=prev_perm;

    return ;
}

void genPerm(int pIdx, si64Matrix &k, si64Matrix &perm, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    int n = k.rows();

    sbMatrix k_bool(n,64);
    arith2bool(pIdx, k, k_bool, enc, eval, runtime);

    si64Matrix k_0(n,1);
    getBitKey(pIdx, k_bool, 0, k_0, enc, eval, runtime);
    
    //DEBUG
    // i64Matrix k_0_plain(n,1);
    // enc.revealAll(runtime, k_0, k_0_plain).get();
    // std::cout << "k_0: " << std::endl;
    // for(size_t i=0;i<n;i++){
    //     std::cout << k_0_plain(i,0) << " ";
    // }
    // std::cout << std::endl;
    //----k_0 correct

    si64Matrix pre_perm(n,1);
    genBitPerm(pIdx, k_0, pre_perm, enc, eval, runtime);
    //DEBUG
    // i64Matrix pre_perm_plain(n,1);
    // enc.revealAll(runtime, pre_perm, pre_perm_plain).get();
    // std::cout << "pre_perm: " << std::endl;
    // for(size_t i=0;i<n;i++){
    //     std::cout << pre_perm_plain(i,0) << " ";
    // }
    // std::cout << std::endl;
    //----pre_perm correct
//d
    for(int d=1;d<8;d++){
        si64Matrix k_j(n,1),k_j_prime(n,1);
        getBitKey(pIdx, k_bool, d, k_j, enc, eval, runtime);
        std::cout << "d: " << d << std::endl;
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
        i64Matrix k_j_prime_plain(n,1);
        enc.revealAll(runtime, k_j_prime, k_j_prime_plain).get();
        std::cout << "k_j_prime: " << std::endl;
        for(size_t i=0;i<n;i++){
            std::cout << k_j_prime_plain(i,0) << " ";
        }
        std::cout << std::endl;
        std::cout.flush();
        //----k_j_prime 

        si64Matrix next_perm(n,1);
        genBitPerm(pIdx, k_j_prime, next_perm, enc, eval, runtime);
        composePerm(pIdx, pre_perm, next_perm, perm, enc, eval, runtime);
        pre_perm=perm;
    }
    
    return;
}




