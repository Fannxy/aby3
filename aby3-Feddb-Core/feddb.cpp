#include "feddb.h"
#include "utils.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3-GORAM-Core/Sort.h"
#include "../aby3-GORAM-Core/Basics.h"

using namespace oc;
using namespace aby3;


//si64matrix: shuffle from genperm
void shuffle(int pIdx, si64Matrix& T, si64Matrix& Tres,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime ){
    
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = T.rows();

    // generate the prev, next - correlated randomness.
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    //shuffle
    si64Matrix next_perm(len,1),prev_perm(len,1);
    prev_perm=T;

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
    Tres=prev_perm;

    return ;

    }


//sbmatrix: effcient-shuffle from GORAM 
void shuffle(int pIdx, sbMatrix& T, sbMatrix &Tres, 
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

        efficient_shuffle(T, pIdx, Tres, enc, eval, runtime);
        return ;
    }

void shuffle(int pIdx, std::vector<si64Matrix>& T, std::vector<si64Matrix>& Tres,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    size_t len = T.size();
    size_t unit_len = T[0].rows();

    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    
    // generate the prev, next - correlated randomness.
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(unit_len, prev_permutation, prevSeed);
    get_permutation(unit_len, next_permutation, nextSeed);

    //shuffle
    std::vector<si64Matrix> next_perm(len),prev_perm(len);
    for(size_t i=0; i<len; i++){
        prev_perm[i].resize(unit_len, 1);
        next_perm[i].resize(unit_len, 1);
        prev_perm[i]=T[i];
    }

    for(size_t i=0;i<3;i++){
        if(pIdx==i){
            permutate(pIdx, prev_perm, next_perm, prev_permutation);
            
        }
        else if(pIdx==((i+2)%3)){
            permutate(pIdx, prev_perm, next_perm, next_permutation);
            
        }
        else if(pIdx==(i+1)%3){
            for(size_t j=0; j<len; j++){
                std::fill_n(next_perm[j].mShares[0].data(),unit_len,0);
                std::fill_n(next_perm[j].mShares[1].data(),unit_len,0);
            }
        }
        
        for(size_t j=0; j<len; j++){
            prev_perm[j]=reshare_matrix(pIdx, next_perm[j], (i+1)%3, enc, runtime);
        }

    }

    for(size_t i=0; i<len; i++){
        Tres[i].resize(unit_len, 1);
        Tres[i]=prev_perm[i];
    }

    return ;

    }


void persist_plain(int pIdx, const std::string &table_name, std::vector<i64Matrix> &T){
    std::string filename = "./aby3-Feddb-tmpfile/plain/" + table_name + "_" + std::to_string(pIdx) + ".txt";
    
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing" << std::endl;
        return;
    }

    int cols = T.size();
    outFile << "Matrix cols: " << cols << std::endl;
    int rows = T[0].rows();
    outFile << "Matrix rows: " << rows << std::endl;


    for(size_t i = 0; i < T.size(); i++){
        outFile << "Matrix_" << i << " values:" << std::endl;
        for (int j = 0; j < rows; j++) {
            outFile << T[i](j, 0) << std::endl;
        }
    }

    outFile.close();

    return;
}

    
void read_plain(int pIdx, const std::string &table_name, std::vector<i64Matrix> &T){
    std::string filename = "./aby3-Feddb-tmpfile/plain/" + table_name + "_" + std::to_string(pIdx) + ".txt";
    
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for reading" << std::endl;
        return;
    }

    std::string dummy;
    int cols;
    inFile >> dummy >> dummy >> cols;
    T.resize(cols);
    
    int rows;
    inFile >> dummy >> dummy >> rows;
    for(size_t i=0; i<T.size(); i++){
        T[i].resize(rows, 1);
        inFile >> dummy >> dummy;
        for (int j = 0; j < rows; j++) {
            inFile >> T[i](j, 0);
        }
    }

    inFile.close();
    return;
}

void oblivious_idx_select(int pIdx, si64Matrix &v, si64Matrix &idx, si64Matrix &result,              
        Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    size_t v_len = v.rows();
    size_t idx_len = idx.rows();

    //step zero: 将v和idx拼接成t
    size_t t_len=v_len+idx_len;
    si64Matrix t(t_len, 1);

    for(size_t i=0; i<v_len; i++){
        if(pIdx == 0){
            t.mShares[0](i, 0) = 0;
            t.mShares[1](i, 0) = i;
        }
        if(pIdx == 1){
            t.mShares[0](i, 0) = 0;
            t.mShares[1](i, 0) = 0;
        }
        if(pIdx == 2){
            t.mShares[0](i, 0) = i;
            t.mShares[1](i, 0) = 0;
        }
    }
    for(size_t i=v_len; i<t_len; i++){
        t.mShares[0](i, 0) = idx.mShares[0](i-v_len, 0);
        t.mShares[1](i, 0) = idx.mShares[1](i-v_len, 0);
    }


    //step one: 对t进行genperm得到sigma
    si64Matrix sigma(t_len, 1);
    genPerm(pIdx, t, sigma, enc, eval, runtime);    

    //step two: prefixsum_{-1}(v) && u
    si64Matrix prefix_inv(v_len, v.cols());
    prefixsum_inv(pIdx, v, prefix_inv);
    si64Matrix u(t_len, v.cols());
    for(size_t i=0; i<v_len; i++){
        u.mShares[0](i, 0) = prefix_inv.mShares[0](i, 0);
        u.mShares[1](i, 0) = prefix_inv.mShares[1](i, 0);
    }
    for(size_t i=v_len; i<t_len; i++){
        u.mShares[0](i, 0) = 0;
        u.mShares[1](i, 0) = 0;
    }

    //DEBUG
    // i64Matrix u_plain(t_len, v.cols());
    // enc.revealAll(runtime, u, u_plain).get();
    // std::cout << "u_plain: " << std::endl;
    // for(size_t i=0; i<t_len; i++){  
    //     std::cout << u_plain(i, 0) << " ";
    // }
    // std::cout << std::endl;  
    // std::cout.flush();      
    //----u_plain correct

    //---------------------------------apply permutation start
    //step three: 对u进行sigma的permutation
        //step 1 [[rsigma]] = pi([[sigma]])
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = sigma.rows();

    //  generate the permutations. pi
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);
    std::vector<size_t> prev_inverse_permutation;
    std::vector<size_t> next_inverse_permutation;
    get_inverse_permutation(prev_permutation, prev_inverse_permutation);
    get_inverse_permutation(next_permutation, next_inverse_permutation);

    //shuffle:pi(sigma)
    si64Matrix rsigma(len,1);
    si64Matrix next_perm(len,1),prev_perm(len,1);
    prev_perm=sigma;
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
    rsigma=prev_perm;

        //step 2 [[u1]] = pi([[u]])
    //shuffle:pi(u)
    si64Matrix u1(len,1);
    prev_perm=u;
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
    u1=prev_perm;

        // Step 3: 恢复rsigma至明文，并得到rsigma^{-1}
    i64Matrix rsigma_plain(len, 1);
    enc.revealAll(runtime, rsigma, rsigma_plain).get();
    
    i64Matrix rsigma_inverse_plain(len, 1);
    permutation_inverse(rsigma_plain, rsigma_inverse_plain);
    
        // Step 5: 根据rsigma对u1进行置换得到u_prime = rsigma (u1)= sigma * pi^{-1} * pi([[u1]]) = sigma(u1)
    si64Matrix u_prime(len, 1);
    permutate(pIdx, u1, u_prime, rsigma_plain);

    //---------------------------------apply permutation finished
    
    //step four: prefixsum_(u_prime)
    si64Matrix u_prime_prefix(len, v.cols());
    prefixsum(pIdx, u_prime, u_prime_prefix);

    //step five: unapply permutation
    //---------------------------------unapply permutation start
        //step 1:[[u1_]]=rsigma^{-1}(u_prime)
    si64Matrix u_1_(len, 1);
    permutate(pIdx, u_prime_prefix, u_1_, rsigma_inverse_plain);

       //step 2: [[u']]=pi^{-1}([[u1_]])
    si64Matrix u_prime_(len, 1);
    //unshuffle
    prev_perm=u_1_;
    for(int id=2;id>=0;id--){
        if(pIdx==id){
            permutate(pIdx, prev_perm, next_perm, prev_inverse_permutation);
            
        }
        else if(pIdx==((id+2)%3)){
            permutate(pIdx, prev_perm, next_perm, next_inverse_permutation);
            
        }
        else if(pIdx==(id+1)%3){
            std::fill_n(next_perm.mShares[0].data(),len,0);
            std::fill_n(next_perm.mShares[1].data(),len,0);

        }
        prev_perm=reshare_matrix(pIdx, next_perm, (id+1)%3, enc, runtime);
        
    }
    u_prime_=prev_perm;
    //---------------------------------unapply permutation finished

    //step six:get result
    result.resize(idx_len, v.cols());
    for (size_t i = 0; i < idx_len; i++) {
        result.mShares[0](i, 0) = u_prime_.mShares[0](i+v_len, 0);
        result.mShares[1](i, 0) = u_prime_.mShares[1](i+v_len, 0);
    }

    return; 
}

//agg: sum(valcol) group by keycol : only val:non-zero entries
//data[0]: keycol & data[1]: valcol
void index_agg(int pIdx, si64Matrix &equalFlag, i64Matrix &idx,std::vector<si64Matrix> &data,std::vector<si64Matrix> &finalRes,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime ){
    
        int rows=data[0].rows();

        si64Matrix keycol,valcol;
        keycol.resize(rows, 1);
        valcol.resize(rows, 1);
        //sorted : acoording to idx
        for(int i=0; i<rows ;i++){
            i64 idx_i=idx(i, 0);
            keycol.mShares[0](i, 0) = data[0].mShares[0](idx_i, 0);
            keycol.mShares[1](i, 0) = data[0].mShares[1](idx_i, 0);
            //keycol(i,0) = data[0](idx_i, 0);
            valcol.mShares[0](i, 0) = data[1].mShares[0](idx_i, 0);
            valcol.mShares[1](i, 0) = data[1].mShares[1](idx_i, 0);
            //valcol(i,0) = data[1](idx_i, 0);
        }

        
        for(int i=0; i<(rows-1) ;i++){
            // valcol(i,0) = leftVal * (1 - eqFlag);
            // valcol(i+1,0) = rightVal + leftVal * (eqFlag);
            si64Matrix leftval(1,1),rightval(1,1);
            leftval.mShares[0](0,0)=valcol.mShares[0](i,0);
            leftval.mShares[1](0,0)=valcol.mShares[1](i,0);
            rightval.mShares[0](0,0)=valcol.mShares[0](i+1,0);
            rightval.mShares[1](0,0)=valcol.mShares[1](i+1,0);

            si64Matrix new_leftval(1,1),new_rightval(1,1);

            i64Matrix one(1,1);
            one(0,0)=1;
            si64Matrix not_eqFlag(1,1),eqFlag(1,1),oneShared(1,1);
            if (pIdx == 0) {
                enc.localIntMatrix(runtime, one, oneShared).get();
            } else {
                enc.remoteIntMatrix(runtime, oneShared).get();
            }

            eqFlag.mShares[0](0,0)=equalFlag.mShares[0](i+1,0);
            eqFlag.mShares[1](0,0)=equalFlag.mShares[1](i+1,0);
            not_eqFlag=oneShared-eqFlag;

            cipher_mul(pIdx, leftval, not_eqFlag, new_leftval, eval, enc, runtime);

            si64Matrix leftval_times_eqFlag(1,1);
            cipher_mul(pIdx, leftval, eqFlag, leftval_times_eqFlag, eval, enc, runtime);
            new_rightval = rightval + leftval_times_eqFlag;

            valcol.mShares[0](i,0) = new_leftval.mShares[0](0,0);
            valcol.mShares[1](i,0) = new_leftval.mShares[1](0,0);
            valcol.mShares[0](i+1,0) = new_rightval.mShares[0](0,0);
            valcol.mShares[1](i+1,0) = new_rightval.mShares[1](0,0);

        }


        sbMatrix valcol_sb(rows,64), zeroFlag_sb(rows,1);
        i64Matrix zero(rows,1);
        for(int i=0;i<rows;i++){
            zero(i,0)=0;
        }
        arith2bool(pIdx, valcol, valcol_sb, enc, eval, runtime);
        bool_cipher_eq(pIdx, valcol_sb, zero, zeroFlag_sb, enc, eval, runtime);


        si64Matrix zeroFlag(rows,1);
        bool2arith(pIdx, zeroFlag_sb, zeroFlag, enc, eval, runtime);


        std::vector<si64Matrix> res(3), shuffledRes(3);
        res[0]=keycol;
        res[1]=valcol;
        res[2]=zeroFlag;

        shuffle(pIdx, res, shuffledRes, enc, eval, runtime);
        zeroFlag=shuffledRes[2];
        i64Matrix zeroFlag_plain(rows,1);
        enc.revealAll(runtime, zeroFlag, zeroFlag_plain).get();

        int resRows=0;
        for(int i=0;i<rows;i++){
            if(zeroFlag_plain(i,0)==0){
                resRows++;
            }
        }

        finalRes[0].resize(resRows, 1);
        finalRes[1].resize(resRows, 1);

        int resIdx=0;
        for(int i=0;i<rows;i++){
            if(zeroFlag_plain(i,0)==0){
                finalRes[0].mShares[0](resIdx,0)=shuffledRes[0].mShares[0](i,0);
                finalRes[0].mShares[1](resIdx,0)=shuffledRes[0].mShares[1](i,0);
                finalRes[1].mShares[0](resIdx,0)=shuffledRes[1].mShares[0](i,0);
                finalRes[1].mShares[1](resIdx,0)=shuffledRes[1].mShares[1](i,0);
                resIdx++;
            }
        }
        
        return;

    }


void group_by_common(int pIdx, si64Matrix &key, si64Matrix &val, 
    si64Matrix &key_g, si64Matrix &val_g, si64Matrix &e, si64Matrix &perm_GN, si64Matrix &key_GN, si64Matrix &key_out,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int rows=key.rows();
    si64Matrix null_share(rows, 1),one_share(rows, 1);
    set_const_share(pIdx, std::numeric_limits<i64>::max(), null_share, enc, eval, runtime);
    set_const_share(pIdx, 1, one_share, enc, eval, runtime);
    
    //step 1
    si64Matrix perm(rows, 1);
    genPerm_kv(pIdx, key, val, perm, enc, eval, runtime);

    //step 2
    key_g.resize(rows, 1);
    val_g.resize(rows, 1);
    applyPerm(pIdx, perm, key, key_g, enc, eval, runtime);
    applyPerm(pIdx, perm, val, val_g, enc, eval, runtime);

    //step 3
    sbMatrix f(rows-1, 64);
    compare_consecutive_rows_arith(pIdx, key_g, f, enc, eval, runtime);
    si64Matrix f_si(rows-1, 1);
    bool2arith(pIdx, f, f_si, enc, eval, runtime);

    //step 4
    si64Matrix e_partial(rows-1, 1);
    e_partial = one_share - f_si;
    e.resize(rows, 1);
    for(size_t i=0;i<rows-1;i++){
        e.mShares[0](i, 0) = e_partial.mShares[0](i, 0);
        e.mShares[1](i, 0) = e_partial.mShares[1](i, 0);
    }
    e.mShares[0](rows-1, 0) = one_share.mShares[0](rows-1, 0) ;
    e.mShares[1](rows-1, 0) = one_share.mShares[1](rows-1, 0) ;
    

    //step 5
    si64Matrix tmp(rows, 1);
    cipher_mul(pIdx, e, key_g-null_share, tmp, eval, enc, runtime);
    key_GN.resize(rows, 1);
    key_GN = tmp + null_share;

    //step 6
    tmp = one_share-e;
    genPerm(pIdx, tmp, perm_GN, enc, eval, runtime);

    //step 7
    key_out.resize(rows, 1);
    applyPerm(pIdx, perm_GN, key_GN, key_out, enc, eval, runtime);

    return;
}

void group_by_common_without_val(int pIdx, si64Matrix &key, 
    si64Matrix &key_g, si64Matrix &e, si64Matrix &perm_GN, si64Matrix &key_GN, si64Matrix &key_out,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int rows=key.rows();
    si64Matrix null_share(rows, 1),one_share(rows, 1);
    set_const_share(pIdx, std::numeric_limits<i64>::max(), null_share, enc, eval, runtime);
    set_const_share(pIdx, 1, one_share, enc, eval, runtime);
    
    //step 1
    si64Matrix perm(rows, 1);
    genPerm(pIdx, key, perm, enc, eval, runtime);

    //step 2
    key_g.resize(rows, 1);
    applyPerm(pIdx, perm, key, key_g, enc, eval, runtime);

    //step 3
    sbMatrix f(rows-1, 64);
    compare_consecutive_rows_arith(pIdx, key_g, f, enc, eval, runtime);
    si64Matrix f_si(rows-1, 1);
    bool2arith(pIdx, f, f_si, enc, eval, runtime);

    //step 4
    si64Matrix e_partial(rows-1, 1);
    e_partial = one_share - f_si;
    e.resize(rows, 1);
    for(size_t i=0;i<rows-1;i++){
        e.mShares[0](i, 0) = e_partial.mShares[0](i, 0);
        e.mShares[1](i, 0) = e_partial.mShares[1](i, 0);
    }
    e.mShares[0](rows-1, 0) = one_share.mShares[0](rows-1, 0) ;
    e.mShares[1](rows-1, 0) = one_share.mShares[1](rows-1, 0) ;
    

    //step 5
    si64Matrix tmp(rows, 1);
    cipher_mul(pIdx, e, key_g-null_share, tmp, eval, enc, runtime);
    key_GN.resize(rows, 1);
    key_GN = tmp + null_share;

    //step 6
    tmp = one_share-e;
    genPerm(pIdx, tmp, perm_GN, enc, eval, runtime);

    //step 7
    key_out.resize(rows, 1);
    applyPerm(pIdx, perm_GN, key_GN, key_out, enc, eval, runtime);

    return;
}

void group_count(int pIdx, si64Matrix &key, 
    si64Matrix &key_out, si64Matrix &c, 
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

    int rows=key.rows();
    key_out.resize(rows, 1);
    c.resize(rows, 1);

    //step 1
    si64Matrix key_g(rows, 1);
    si64Matrix e(rows, 1);
    si64Matrix perm_GN(rows, 1);
    si64Matrix key_GN(rows, 1);
    group_by_common_without_val(pIdx, key, key_g, e, perm_GN, key_GN, key_out, enc, eval, runtime);

    //step 2
    si64Matrix x(rows, 1);
    i64Matrix idx(rows, 1);
    si64Matrix idxShared(rows, 1);
    for(size_t i=0;i<rows;i++){
        idx(i, 0) = i+1;
    }
    if(pIdx==0){
        enc.localIntMatrix(runtime, idx, idxShared).get();
    }
    else{
        enc.remoteIntMatrix(runtime, idxShared).get();
    }
    
    si64Matrix mShared(rows, 1);
    set_const_share(pIdx, rows, mShared, enc, eval, runtime);
    si64Matrix tmp(rows, 1);
    cipher_mul(pIdx, e, idxShared-mShared, tmp, eval, enc, runtime);
    x = tmp + mShared;

    //step 3
    si64Matrix y(rows, 1);
    applyPerm(pIdx, perm_GN, x, y, enc, eval, runtime);

    //step 4
    prefixsum_inv(pIdx, y, c);

    return ;
}


void group_sum(int pIdx, si64Matrix &key, si64Matrix &val, 
    si64Matrix &key_out, si64Matrix &sum,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){ 
    
    int rows=key.rows();
    key_out.resize(rows, 1);
    sum.resize(rows, 1);

    si64Matrix key_g(rows, 1);
    si64Matrix val_g(rows, 1);
    si64Matrix e(rows, 1);
    si64Matrix perm_GN(rows, 1);
    si64Matrix key_GN(rows, 1);
    //step 1
    group_by_common(pIdx, key, val, key_g, val_g, e, perm_GN, key_GN, key_out, enc, eval, runtime);

    //step 2
    si64Matrix w(rows, 1);
    prefixsum(pIdx, val_g, w);

    //step 3
    si64Matrix w_m(rows, 1);
    si64Matrix x(rows,1);
    for(size_t i=0;i<rows;i++){
        w_m.mShares[0](i, 0) = w.mShares[0](rows-1, 0);
        w_m.mShares[1](i, 0) = w.mShares[1](rows-1, 0);
    }
    si64Matrix tmp(rows, 1);
    cipher_mul(pIdx, e, w-w_m, tmp, eval, enc, runtime);
    x = tmp + w_m;

    //step 4
    si64Matrix y(rows, 1);
    applyPerm(pIdx, perm_GN, x, y, enc, eval, runtime);

    //step 5
    prefixsum_inv(pIdx, y, sum);
    
    return ;
}

void group_max(int pIdx, si64Matrix &key, si64Matrix &val, 
    si64Matrix &key_out, si64Matrix &max,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int rows=key.rows();
    key_out.resize(rows, 1);
    max.resize(rows, 1);
    
    si64Matrix key_g(rows, 1);
    si64Matrix val_g(rows, 1);
    si64Matrix e(rows, 1);
    si64Matrix perm_GN(rows, 1);
    si64Matrix key_GN(rows, 1);
    //step 1
    group_by_common(pIdx, key, val, key_g, val_g, e, perm_GN, key_GN, key_out, enc, eval, runtime);

    //step 2
    si64Matrix x(rows, 1);
    si64Matrix zero(rows, 1);
    si64Matrix tmp(rows, 1);
    set_const_share(pIdx, 0, zero, enc, eval, runtime);
    cipher_mul(pIdx, e, val_g-zero, tmp, eval, enc, runtime);
    x = tmp + zero;

    //step 3
    applyPerm(pIdx, perm_GN, x, max, enc, eval, runtime);
    
    return ;
}

void group_min(int pIdx, si64Matrix &key, si64Matrix &val, 
    si64Matrix &key_out, si64Matrix &min,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int rows=key.rows();
    key_out.resize(rows, 1);
    min.resize(rows, 1);

    si64Matrix key_g(rows, 1);
    si64Matrix val_g(rows, 1);
    si64Matrix e(rows, 1);
    si64Matrix perm_GN(rows, 1);
    si64Matrix key_GN(rows, 1);
    //step 1
    group_by_common(pIdx, key, val, key_g, val_g, e, perm_GN, key_GN, key_out, enc, eval, runtime);

    //step 2.1
    si64Matrix e_prime(rows, 1);
    for(size_t i=1;i<rows;i++){
        e_prime.mShares[0](i, 0) = e.mShares[0](i-1, 0);
        e_prime.mShares[1](i, 0) = e.mShares[1](i-1, 0);
    }
    switch(pIdx){
        case 0:
            e_prime.mShares[0](0, 0) = 0;
            e_prime.mShares[1](0, 0) = 0;
            break;
        case 1:
            e_prime.mShares[0](0, 0) = 1;
            e_prime.mShares[1](0, 0) = 0;
            break;
        case 2:
            e_prime.mShares[0](0, 0) = 0;
            e_prime.mShares[1](0, 0) = 1;
            break;
        default:
            THROW_RUNTIME_ERROR("group_min: pIdx out of range.");
    }

    
    si64Matrix g(rows, 1);
    si64Matrix one(rows, 1);
    set_const_share(pIdx, 1, one, enc, eval, runtime);
    g = one - e_prime;

    //step 2.2
    si64Matrix x(rows, 1);
    si64Matrix tmp(rows, 1);
    si64Matrix zero(rows, 1);
    set_const_share(pIdx, 0, zero, enc, eval, runtime);
    cipher_mul(pIdx, e_prime, val_g-zero, tmp, eval, enc, runtime);
    x = tmp + zero;

    //step 3
    si64Matrix perm_out(rows, 1);
    genPerm(pIdx, g, perm_out, enc, eval, runtime);

    //step 4
    applyPerm(pIdx, perm_out, x, min, enc, eval, runtime);
    return ;
}

//subfunc for join
void augment_table(int pIdx, std::vector<si64Matrix> &T_1, std::vector<si64Matrix> &T_2,
    std::vector<sbMatrix> &T_1_auged, std::vector<sbMatrix> &T_2_auged,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    size_t len1 = T_1[0].rows();
    size_t len2 = T_2[0].rows();

    //step 1:T_c: concatenate T_1 and T_2
    std::vector<si64Matrix> T_c(2);
    for(size_t i=0;i<2;i++){
        T_c[i].resize(len1+len2, 64);
        std::memcpy(T_c[i].mShares[0].data(), T_1[i].mShares[0].data(), len1 * sizeof(T_1[i].mShares[0](0, 0)));
        std::memcpy(T_c[i].mShares[1].data(), T_1[i].mShares[1].data(), len1 * sizeof(T_1[i].mShares[1](0, 0)));
        std::memcpy(T_c[i].mShares[0].data() + len1, T_2[i].mShares[0].data(), len2 * sizeof(T_2[i].mShares[0](0, 0)));
        std::memcpy(T_c[i].mShares[1].data() + len1, T_2[i].mShares[1].data(), len2 * sizeof(T_2[i].mShares[1](0, 0)));

    }

    //T_c_bool: convert T_c to bool
    std::vector<sbMatrix> T_c_bool(3);
    for(size_t i=0;i<2;i++){
        T_c_bool[i].resize(len1+len2, 64);
        arith2bool(pIdx, T_c[i], T_c_bool[i], enc, eval, runtime);
    }

    //T_c_bool[2]: tid
    T_c_bool[2].resize(len1+len2, 1);
    sbMatrix zero(len1, 1);
    sbMatrix one(len2, 1);
    bool_init_false(pIdx, zero);
    bool_init_true(pIdx, one);
    
    std::memcpy(T_c_bool[2].mShares[0].data(), zero.mShares[0].data(), len1 * sizeof(zero.mShares[0](0, 0)));
    std::memcpy(T_c_bool[2].mShares[1].data(), zero.mShares[1].data(), len1 * sizeof(zero.mShares[1](0, 0)));
    std::memcpy(T_c_bool[2].mShares[0].data() + len1, one.mShares[0].data(), len2 * sizeof(one.mShares[0](0, 0)));
    std::memcpy(T_c_bool[2].mShares[1].data() + len1, one.mShares[1].data(), len2 * sizeof(one.mShares[1](0, 0)));
    
    //step 2: sort T_c: j, tid(转成了bool进行拼接排序)
    sbMatrix entry_key(len1+len2, 64);
    for(size_t i=0; i<len1+len2; i++){
        entry_key.mShares[0](i, 0) = (T_c_bool[0].mShares[0](i, 0) << 33) | ((T_c_bool[2].mShares[0](i, 0) & 0x3) << 31) | (T_c_bool[1].mShares[0](i, 0)&0x7FFFFFFF);
        entry_key.mShares[1](i, 0) = (T_c_bool[0].mShares[1](i, 0) << 33) | ((T_c_bool[2].mShares[1](i, 0) & 0x3) << 31) | (T_c_bool[1].mShares[1](i, 0)&0x7FFFFFFF);
    }

    sbMatrix entry_key_sorted(len1+len2, 64);
    odd_even_merge_sort(entry_key, entry_key_sorted, pIdx, enc, eval, runtime);

    sbMatrix j_sb(len1+len2, 64);
    sbMatrix tid_sb(len1+len2, 1);
    sbMatrix d_sb(len1+len2, 64);

    for(size_t i=0; i< len1+len2; i++){
        j_sb.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 33) & 0x7FFFFFFF;
        j_sb.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 33) & 0x7FFFFFFF;
        
        tid_sb.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 31) & 0x3;
        tid_sb.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 31) & 0x3;

        d_sb.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, 0)& 0x7FFFFFFF;
        d_sb.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, 0)& 0x7FFFFFFF;

    }

    //step 3: full-dimension
    sbMatrix alpha_1(len1+len2, 64), alpha_2(len1+len2, 64);
    //downward scan: 
    for(size_t i=0; i<len1+len2; i++){
        sbMatrix j_i(1, 64), tid_i(1, 64);
        j_i.mShares[0](0, 0) = j_sb.mShares[0](i, 0);
        j_i.mShares[1](0, 0) = j_sb.mShares[1](i, 0);
        tid_i.mShares[0](0, 0) = tid_sb.mShares[0](i, 0);
        tid_i.mShares[1](0, 0) = tid_sb.mShares[1](i, 0);

        //join_key_same
        sbMatrix same_attr(1, 1), not_same_attr(1,1);
        if(i == 0){
            bool_init_false(pIdx, same_attr);
        } else {
            sbMatrix j_i_prev(1, 64);
            j_i_prev.mShares[0](0, 0) = j_sb.mShares[0](i-1, 0);
            j_i_prev.mShares[1](0, 0) = j_sb.mShares[1](i-1, 0);
            bool_cipher_eq(pIdx, j_i, j_i_prev, same_attr, enc, eval, runtime);
        }
        sbMatrix one_bit(1, 1);
        bool_init_true(pIdx, one_bit);
        bool_cipher_sub(pIdx, one_bit, same_attr, not_same_attr, enc, eval, runtime);

        //tid
        sbMatrix zero64(1, 64);
        sbMatrix one64(1, 64);
        bool_init_false(pIdx, zero64);
        bool_init_true(pIdx, one64);
        sbMatrix is_table_1(1, 1);
        bool_cipher_eq(pIdx, tid_i, zero64, is_table_1, enc, eval, runtime);
        sbMatrix is_table_2(1, 1);
        bool_cipher_eq(pIdx, tid_i, one64, is_table_2, enc, eval, runtime);

        //condition
        sbMatrix cond1, cond2, cond3, cond4;
        bool_cipher_and(pIdx, is_table_1, not_same_attr, cond1, enc, eval, runtime);
        bool_cipher_and(pIdx, is_table_1, same_attr, cond2, enc, eval, runtime);
        bool_cipher_and(pIdx, is_table_2, not_same_attr, cond3, enc, eval, runtime);
        bool_cipher_and(pIdx, is_table_2, same_attr, cond4, enc, eval, runtime);

        // alpha1
        sbMatrix current_alpha_1(1, 64);
        if(i == 0) {
            current_alpha_1 = zero64;
        } else {
            current_alpha_1.mShares[0](0, 0) = alpha_1.mShares[0](i-1, 0);
            current_alpha_1.mShares[1](0, 0) = alpha_1.mShares[1](i-1, 0);
        }

        sbMatrix alpha_1_plus(1, 64);
        bool_cipher_add(pIdx, current_alpha_1, one64, alpha_1_plus, enc, eval, runtime);
        
        sbMatrix term1(1, 64), term2(1, 64), term3(1, 64), term4(1, 64);
        cond1.resize(1, 64);
        cond2.resize(1, 64);
        cond3.resize(1, 64);
        cond4.resize(1, 64);
        cond1.mShares[0](0, 0) = (cond1.mShares[0](0, 0) == 1) ? -1 : -0;
        cond1.mShares[1](0, 0) = (cond1.mShares[1](0, 0) == 1) ? -1 : -0;
        cond2.mShares[0](0, 0) = (cond2.mShares[0](0, 0) == 1) ? -1 : -0;
        cond2.mShares[1](0, 0) = (cond2.mShares[1](0, 0) == 1) ? -1 : -0;
        cond3.mShares[0](0, 0) = (cond3.mShares[0](0, 0) == 1) ? -1 : -0;
        cond3.mShares[1](0, 0) = (cond3.mShares[1](0, 0) == 1) ? -1 : -0;
        cond4.mShares[0](0, 0) = (cond4.mShares[0](0, 0) == 1) ? -1 : -0;
        cond4.mShares[1](0, 0) = (cond4.mShares[1](0, 0) == 1) ? -1 : -0;

        bool_cipher_and(pIdx, cond1, one64, term1, enc, eval, runtime);
        bool_cipher_and(pIdx, cond2, alpha_1_plus, term2, enc, eval, runtime);
        bool_cipher_and(pIdx, zero64, cond3, term3, enc, eval, runtime);
        bool_cipher_and(pIdx, current_alpha_1, cond4, term4, enc, eval, runtime);
       
        alpha_1.mShares[0](i, 0) = term1.mShares[0](0, 0) ^ term2.mShares[0](0, 0)^ term3.mShares[0](0, 0)^ term4.mShares[0](0, 0);
        alpha_1.mShares[1](i, 0) = term1.mShares[1](0, 0) ^ term2.mShares[1](0, 0)^ term3.mShares[1](0, 0)^ term4.mShares[1](0, 0);
        

        // alpha2
        sbMatrix current_alpha_2(1, 64);
        if(i == 0) {
            current_alpha_2 = zero64;
        } else {
            current_alpha_2.mShares[0](0, 0) = alpha_2.mShares[0](i-1, 0);
            current_alpha_2.mShares[1](0, 0) = alpha_2.mShares[1](i-1, 0);
        }

        sbMatrix alpha_2_plus(1, 64);
        bool_cipher_add(pIdx, current_alpha_2, one64, alpha_2_plus, enc, eval, runtime);

        sbMatrix term5(1, 64), term6(1, 64), term7(1, 64), term8(1, 64);
        bool_cipher_and(pIdx, zero64, cond1, term5, enc, eval, runtime);
        bool_cipher_and(pIdx, zero64, cond2, term6, enc, eval, runtime);
        bool_cipher_and(pIdx, one64, cond3, term7, enc, eval, runtime);
        bool_cipher_and(pIdx, alpha_2_plus, cond4, term8, enc, eval, runtime);
        

        alpha_2.mShares[0](i, 0) = term5.mShares[0](0, 0) ^ term6.mShares[0](0, 0)^ term7.mShares[0](0, 0)^ term8.mShares[0](0, 0);
        alpha_2.mShares[1](i, 0) = term5.mShares[1](0, 0) ^ term6.mShares[1](0, 0)^ term7.mShares[1](0, 0)^ term8.mShares[1](0, 0);
        
    }
    

    //upward scan: 
    for(int i=len1+len2-1; i>=0; i--){
        sbMatrix j_i(1, 64), tid_i(1, 64);
        j_i.mShares[0](0, 0) = j_sb.mShares[0](i, 0);
        j_i.mShares[1](0, 0) = j_sb.mShares[1](i, 0);
        tid_i.mShares[0](0, 0) = tid_sb.mShares[0](i, 0);
        tid_i.mShares[1](0, 0) = tid_sb.mShares[1](i, 0);

        //join_key_same
        sbMatrix same_attr(1, 1), not_same_attr(1,1);
        if(i == len1+len2-1){
            bool_init_false(pIdx, same_attr);
        } else {
            sbMatrix j_i_prev(1, 64);
            j_i_prev.mShares[0](0, 0) = j_sb.mShares[0](i+1, 0);
            j_i_prev.mShares[1](0, 0) = j_sb.mShares[1](i+1, 0);
            bool_cipher_eq(pIdx, j_i, j_i_prev, same_attr, enc, eval, runtime);
        }
        sbMatrix one_bit(1, 1);
        bool_init_true(pIdx, one_bit);
        bool_cipher_sub(pIdx, one_bit, same_attr, not_same_attr, enc, eval, runtime);
        
        //tid
        sbMatrix zero64(1, 64);
        sbMatrix one64(1, 64);
        bool_init_false(pIdx, zero64);
        bool_init_true(pIdx, one64);
        sbMatrix is_table_1(1, 1);
        bool_cipher_eq(pIdx, tid_i, zero64, is_table_1, enc, eval, runtime);
        sbMatrix is_table_2(1, 1);
        bool_cipher_eq(pIdx, tid_i, one64, is_table_2, enc, eval, runtime);

        //condition
        sbMatrix cond1, cond2, cond3, cond4;
        bool_cipher_and(pIdx, is_table_1, not_same_attr, cond1, enc, eval, runtime);
        bool_cipher_and(pIdx, is_table_1, same_attr, cond2, enc, eval, runtime);
        bool_cipher_and(pIdx, is_table_2, not_same_attr, cond3, enc, eval, runtime);
        bool_cipher_and(pIdx, is_table_2, same_attr, cond4, enc, eval, runtime);   
        cond1.resize(1, 64);
        cond2.resize(1, 64);
        cond3.resize(1, 64);
        cond4.resize(1, 64);
        cond1.mShares[0](0, 0) = (cond1.mShares[0](0, 0) == 1) ? -1 : -0;
        cond1.mShares[1](0, 0) = (cond1.mShares[1](0, 0) == 1) ? -1 : -0;
        cond2.mShares[0](0, 0) = (cond2.mShares[0](0, 0) == 1) ? -1 : -0;
        cond2.mShares[1](0, 0) = (cond2.mShares[1](0, 0) == 1) ? -1 : -0;
        cond3.mShares[0](0, 0) = (cond3.mShares[0](0, 0) == 1) ? -1 : -0;
        cond3.mShares[1](0, 0) = (cond3.mShares[1](0, 0) == 1) ? -1 : -0;
        cond4.mShares[0](0, 0) = (cond4.mShares[0](0, 0) == 1) ? -1 : -0;
        cond4.mShares[1](0, 0) = (cond4.mShares[1](0, 0) == 1) ? -1 : -0;

        // alpha1
        sbMatrix pre_alpha_1(1, 64);
        if(i == len1+len2-1) {
            pre_alpha_1.mShares[0](0, 0) = alpha_1.mShares[0](len1+len2-1, 0);
            pre_alpha_1.mShares[1](0, 0) = alpha_1.mShares[1](len1+len2-1, 0);
        } else {
            pre_alpha_1.mShares[0](0, 0) = alpha_1.mShares[0](i+1, 0);
            pre_alpha_1.mShares[1](0, 0) = alpha_1.mShares[1](i+1, 0);
        }

        sbMatrix current_alpha_1(1, 64);
        current_alpha_1.mShares[0](0, 0) = alpha_1.mShares[0](i, 0);
        current_alpha_1.mShares[1](0, 0) = alpha_1.mShares[1](i, 0);

        sbMatrix term1(1, 64), term2(1, 64), term3(1, 64), term4(1, 64);
        bool_cipher_and(pIdx, current_alpha_1, cond1, term1, enc, eval, runtime);
        bool_cipher_and(pIdx, pre_alpha_1, cond2, term2, enc, eval, runtime);
        bool_cipher_and(pIdx, current_alpha_1, cond3, term3, enc, eval, runtime);
        bool_cipher_and(pIdx, current_alpha_1, cond4, term4, enc, eval, runtime);

        alpha_1.mShares[0](i, 0) = term1.mShares[0](0, 0) ^ term2.mShares[0](0, 0)^ term3.mShares[0](0, 0)^ term4.mShares[0](0, 0);
        alpha_1.mShares[1](i, 0) = term1.mShares[1](0, 0) ^ term2.mShares[1](0, 0)^ term3.mShares[1](0, 0)^ term4.mShares[1](0, 0);
        
        // alpha2
        sbMatrix pre_alpha_2(1, 64);
        if(i == len1+len2-1) {
            pre_alpha_2.mShares[0](0, 0) = alpha_2.mShares[0](len1+len2-1, 0);
            pre_alpha_2.mShares[1](0, 0) = alpha_2.mShares[1](len1+len2-1, 0);
        } else {
            pre_alpha_2.mShares[0](0, 0) = alpha_2.mShares[0](i+1, 0);
            pre_alpha_2.mShares[1](0, 0) = alpha_2.mShares[1](i+1, 0);
        }

        sbMatrix current_alpha_2(1, 64);
        current_alpha_2.mShares[0](0, 0) = alpha_2.mShares[0](i, 0);
        current_alpha_2.mShares[1](0, 0) = alpha_2.mShares[1](i, 0);

        sbMatrix term5(1, 64), term6(1, 64), term7(1, 64), term8(1, 64);
        bool_cipher_and(pIdx, pre_alpha_2, cond1, term5, enc, eval, runtime);
        bool_cipher_and(pIdx, pre_alpha_2, cond2, term6, enc, eval, runtime);
        bool_cipher_and(pIdx, current_alpha_2, cond3, term7, enc, eval, runtime);
        bool_cipher_and(pIdx, pre_alpha_2, cond4, term8, enc, eval, runtime);

        alpha_2.mShares[0](i, 0) = term5.mShares[0](0, 0) ^ term6.mShares[0](0, 0)^ term7.mShares[0](0, 0)^ term8.mShares[0](0, 0);
        alpha_2.mShares[1](i, 0) = term5.mShares[1](0, 0) ^ term6.mShares[1](0, 0)^ term7.mShares[1](0, 0)^ term8.mShares[1](0, 0);
        
    }

    // // //DEBUG
    // i64Matrix alpha_1_plain(len1+len2, 1), alpha_2_plain(len1+len2, 1);
    // enc.revealAll(runtime, alpha_1, alpha_1_plain).get();
    // enc.revealAll(runtime, alpha_2, alpha_2_plain).get();
    // if(pIdx == 0){
    //     std::cout << "alpha_1_plain: " << std::endl;
    //     for(size_t i=0; i<len1+len2; i++){
    //         std::cout << alpha_1_plain(i, 0) << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "alpha_2_plain: " << std::endl;
    //     for(size_t i=0; i<len1+len2; i++){
    //         std::cout << alpha_2_plain(i, 0) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // //-------alpha_1、alpha_2 correct

    //step 4: sort by: tid
    for(size_t i=0; i<len1+len2; i++){
        //tid如果只取一位的话，无符号数的最高位特殊，排序是反的
        entry_key.mShares[0](i, 0) = (tid_sb.mShares[0](i, 0)& 0x3) << 62 | (j_sb.mShares[0](i, 0) & 0x1FFFFFF) << 37 | (d_sb.mShares[0](i, 0) & 0x1FFFFFF) << 12 | (alpha_1.mShares[0](i, 0) & 0x3F) << 6 | (alpha_2.mShares[0](i, 0) & 0x3F) ;
        entry_key.mShares[1](i, 0) = (tid_sb.mShares[1](i, 0)& 0x3) << 62 | (j_sb.mShares[1](i, 0) & 0x1FFFFFF) << 37 | (d_sb.mShares[1](i, 0) & 0x1FFFFFF) << 12 | (alpha_1.mShares[1](i, 0) & 0x3F) << 6 | (alpha_2.mShares[1](i, 0) & 0x3F) ;
    }

    odd_even_merge_sort(entry_key, entry_key_sorted, pIdx, enc, eval, runtime);

    sbMatrix j_sorted(len1+len2, 64), d_sorted(len1+len2, 64),tid_sorted(len1+len2, 64);
    sbMatrix alpha_1_sorted(len1+len2, 64), alpha_2_sorted(len1+len2, 64);
 
    for(size_t i=0; i< len1+len2; i++){
        j_sorted.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 37) & 0x1FFFFFF;
        j_sorted.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 37) & 0x1FFFFFF;
        
        d_sorted.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 12) & 0x1FFFFFF;
        d_sorted.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 12) & 0x1FFFFFF;

        tid_sorted.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 62) & 0x3;
        tid_sorted.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 62) & 0x3;

        alpha_1_sorted.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 6) & 0x3F;
        alpha_1_sorted.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 6) & 0x3F;
        alpha_2_sorted.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) & 0x3F);
        alpha_2_sorted.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) & 0x3F);

    }
   
    //step 5: abstract from T_c
    T_1_auged.resize(4);
    T_2_auged.resize(4);
    for(size_t i=0; i<4; i++){
        T_1_auged[i].resize(len1, 64);
        T_2_auged[i].resize(len2, 64);
    }
    
    std::memcpy(T_1_auged[0].mShares[0].data(), j_sorted.mShares[0].data(), len1 * sizeof(j_sorted.mShares[0](0, 0)));
    std::memcpy(T_1_auged[0].mShares[1].data(), j_sorted.mShares[1].data(), len1 * sizeof(j_sorted.mShares[1](0, 0)));
    std::memcpy(T_1_auged[1].mShares[0].data(), d_sorted.mShares[0].data(), len1 * sizeof(d_sorted.mShares[0](0, 0)));
    std::memcpy(T_1_auged[1].mShares[1].data(), d_sorted.mShares[1].data(), len1 * sizeof(d_sorted.mShares[1](0, 0)));
    std::memcpy(T_1_auged[2].mShares[0].data(), alpha_1_sorted.mShares[0].data(), len1 * sizeof(alpha_1_sorted.mShares[0](0, 0)));
    std::memcpy(T_1_auged[2].mShares[1].data(), alpha_1_sorted.mShares[1].data(), len1 * sizeof(alpha_1_sorted.mShares[1](0, 0)));
    std::memcpy(T_1_auged[3].mShares[0].data(), alpha_2_sorted.mShares[0].data(), len1 * sizeof(alpha_2_sorted.mShares[0](0, 0)));
    std::memcpy(T_1_auged[3].mShares[1].data(), alpha_2_sorted.mShares[1].data(), len1 * sizeof(alpha_2_sorted.mShares[1](0, 0)));

    std::memcpy(T_2_auged[0].mShares[0].data(), j_sorted.mShares[0].data()+len1, len2 * sizeof(j_sorted.mShares[0](0, 0)));
    std::memcpy(T_2_auged[0].mShares[1].data(), j_sorted.mShares[1].data()+len1, len2 * sizeof(j_sorted.mShares[1](0, 0)));
    std::memcpy(T_2_auged[1].mShares[0].data(), d_sorted.mShares[0].data()+len1, len2 * sizeof(d_sorted.mShares[0](0, 0)));
    std::memcpy(T_2_auged[1].mShares[1].data(), d_sorted.mShares[1].data()+len1, len2 * sizeof(d_sorted.mShares[1](0, 0)));
    std::memcpy(T_2_auged[2].mShares[0].data(), alpha_1_sorted.mShares[0].data()+len1, len2 * sizeof(alpha_1_sorted.mShares[0](0, 0)));
    std::memcpy(T_2_auged[2].mShares[1].data(), alpha_1_sorted.mShares[1].data()+len1, len2 * sizeof(alpha_1_sorted.mShares[1](0, 0)));
    std::memcpy(T_2_auged[3].mShares[0].data(), alpha_2_sorted.mShares[0].data()+len1, len2 * sizeof(alpha_2_sorted.mShares[0](0, 0)));
    std::memcpy(T_2_auged[3].mShares[1].data(), alpha_2_sorted.mShares[1].data()+len1, len2 * sizeof(alpha_2_sorted.mShares[1](0, 0)));
    
    return ;
}

void oblivious_expand(int pIdx, std::vector<sbMatrix> &T, std::vector<sbMatrix> &A, i64 tid,
    i64Matrix &s_plain,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    
    int n = T[0].rows();
    //step 1: 计算f(x)和s
    sbMatrix g(n, 64);
    if(tid == 0){
        g = T[3];
    } else {
        g = T[2];
    }

    sbMatrix zero(n, 64);
    sbMatrix flag(n, 1);
    bool_init_false(pIdx, zero);
    bool_cipher_eq(pIdx, g, zero, flag, enc, eval, runtime);

    //flag标记了是否g(x)=0，不用把j,d设为空值，这样后续排序不对

    //f(x)
    sbMatrix fx(n, 64);
    sbMatrix s(1, 64);
    bool_init_false(pIdx, s);
    for(size_t i=0; i<n; i++){
        sbMatrix f_i(1, 64), g_i(1, 64), zero_i(1, 64);
        sbMatrix flag_i(1, 1);
        bool_init_false(pIdx, zero_i);
        flag_i.mShares[0](0, 0) = flag.mShares[0](i, 0);
        flag_i.mShares[1](0, 0) = flag.mShares[1](i, 0);
       
        bool_cipher_selector(pIdx, flag_i, zero_i, s, f_i, enc, eval, runtime);

        fx.mShares[0](i, 0) = f_i.mShares[0](0, 0);
        fx.mShares[1](i, 0) = f_i.mShares[1](0, 0);
        g_i.mShares[0](0, 0) = g.mShares[0](i, 0);
        g_i.mShares[1](0, 0) = g.mShares[1](i, 0);
        
        bool_cipher_add(pIdx, s, g_i, s, enc, eval, runtime);
    }

    s_plain.resize(1, 1);
    enc.revealAll(runtime, s, s_plain).get();

    // //DEBUG
    // i64Matrix fx_plain(n, 1);
    // enc.revealAll(runtime, fx, fx_plain).get();
    // if(pIdx == 0){
    //     std::cout << "fx_plain: " << std::endl;
    //     for(size_t i=0; i<n; i++){
    //         std::cout << fx_plain(i, 0) << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // //-----correct

    //step 2: obli_distribute
    A.resize(4);
    sbMatrix flag_sorted_auged(s_plain(0,0), 1);
    oblivious_distribute(pIdx, T, flag, fx, s_plain, A, flag_sorted_auged, enc, eval, runtime);

    //step 3:fill down
    std::vector<sbMatrix> px(4);
    for(size_t k=0; k<4; k++){
        px[k].resize(1, 64);
        bool_init_i64(pIdx, std::numeric_limits<i64>::max(), px[k], enc, eval, runtime);
    }

    for(size_t i=0; i<s_plain(0,0); i++){
        sbMatrix flag_i(1, 1);
        flag_i.mShares[0](0, 0) = flag_sorted_auged.mShares[0](i, 0);
        flag_i.mShares[1](0, 0) = flag_sorted_auged.mShares[1](i, 0);

        sbMatrix one_share(1, 1);
        bool_init_i64(pIdx, 1, one_share, enc, eval, runtime);

        sbMatrix cond(1, 1);
        bool_cipher_eq(pIdx, flag_i, one_share, cond, enc, eval, runtime);

        for(size_t k=0; k<4; k++){
            sbMatrix A_k_i(1, 64);
            sbMatrix A_k_i_new(1, 64);
            A_k_i.mShares[0](0, 0) = A[k].mShares[0](i, 0);
            A_k_i.mShares[1](0, 0) = A[k].mShares[1](i, 0);

            bool_cipher_selector(pIdx, cond, px[k], A_k_i, A_k_i_new, enc, eval, runtime);

            A[k].mShares[0](i, 0) = A_k_i_new.mShares[0](0, 0);
            A[k].mShares[1](i, 0) = A_k_i_new.mShares[1](0, 0);
            px[k] = A_k_i_new;
        }

    }

    //     //DEBUG
    // for(size_t i=0;i<4;i++){
    //     i64Matrix A_plain(s_plain(0,0), 1);
    //     enc.revealAll(runtime, A[i], A_plain).get();
    //     if(pIdx == 0){
    //         std::cout << "A[" << i << "]_plain: " << std::endl;
    //         for(size_t j=0; j<s_plain(0,0); j++){
    //             std::cout << A_plain(j, 0) << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }
    // //-----correct



    return ;
}

void oblivious_distribute(int pIdx, std::vector<sbMatrix> &T_prime, sbMatrix &flag, sbMatrix &fx, i64Matrix &s_plain, 
    std::vector<sbMatrix> &A, sbMatrix &flag_sorted_auged,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int n = T_prime[0].rows();
    i64 m = s_plain(0, 0);

    for(size_t i=0;i<4;i++){
        T_prime[i].resize(n, 64);
    }

    //step 1: sort T_prime(转成了bool进行拼接排序)
    sbMatrix entry_key(n, 64);
    for(size_t i=0; i<n; i++){
        entry_key.mShares[0](i, 0) = 0 | (flag.mShares[0](i, 0) & 0x1) << 62| (T_prime[0].mShares[0](i, 0) & 0x7FFFF)<< 43| (T_prime[1].mShares[0](i, 0) & 0x7FFFF) << 24 | (fx.mShares[0](i, 0) & 0xFFF) << 12 | (T_prime[2].mShares[0](i, 0) & 0x3F) << 6 | (T_prime[3].mShares[0](i, 0) & 0x3F);
        entry_key.mShares[1](i, 0) = 0 | (flag.mShares[1](i, 0) & 0x1) << 62| (T_prime[0].mShares[1](i, 0) & 0x7FFFF)<< 43| (T_prime[1].mShares[1](i, 0) & 0x7FFFF) << 24 | (fx.mShares[1](i, 0) & 0xFFF) << 12 | (T_prime[2].mShares[1](i, 0) & 0x3F) << 6 | (T_prime[3].mShares[1](i, 0) & 0x3F);
    }

    sbMatrix entry_key_sorted(n, 64);
    odd_even_merge_sort(entry_key, entry_key_sorted, pIdx, enc, eval, runtime);

    sbMatrix j(n, 64), d(n, 64), fx_prime(n, 64), flag_prime(n, 1);
    sbMatrix alpha_1(n, 64), alpha_2(n, 64);
    for(size_t i=0; i< n; i++){
        flag_prime.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 62) & 0x1;
        flag_prime.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 62) & 0x1;

        j.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 43) & 0x7FFFF;
        j.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 43) & 0x7FFFF;
        
        d.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 24) & 0x7FFFF;
        d.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 24) & 0x7FFFF;

        fx_prime.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 12) & 0xFFF;
        fx_prime.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 12) & 0xFFF;

        alpha_1.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 6) & 0x3F;
        alpha_1.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 6) & 0x3F;

        alpha_2.mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) ) & 0x3F;
        alpha_2.mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) ) & 0x3F;
    }

    std::vector<sbMatrix> T_sorted(4);
    for(size_t i=0; i<4; i++){
        T_sorted[i].resize(n, 64);
    }
    T_sorted[0] = j ;
    T_sorted[1] = d;
    T_sorted[2] = alpha_1;
    T_sorted[3] = alpha_2;

    sbMatrix fx_sorted_auged(m, 64);
    flag_sorted_auged.resize(m, 1);
    for(size_t i=0; i<4; i++){
        A[i].resize(m, 64);
    }
    
    //step 2 填充空值得到A, fx_sorted_auged
    if(m >= n){
        for(size_t i=0; i<4; i++){
            std::memcpy(A[i].mShares[0].data(), T_sorted[i].mShares[0].data(), n * sizeof(T_sorted[i].mShares[0](0, 0)));
            std::memcpy(A[i].mShares[1].data(), T_sorted[i].mShares[1].data(), n * sizeof(T_sorted[i].mShares[1](0, 0)));
        }
        std::memcpy(fx_sorted_auged.mShares[0].data(), fx_prime.mShares[0].data(), n * sizeof(fx_prime.mShares[0](0, 0)));
        std::memcpy(fx_sorted_auged.mShares[1].data(), fx_prime.mShares[1].data(), n * sizeof(fx_prime.mShares[1](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[0].data(), flag_prime.mShares[0].data(), n * sizeof(flag_prime.mShares[0](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[1].data(), flag_prime.mShares[1].data(), n * sizeof(flag_prime.mShares[1](0, 0)));

        sbMatrix null_share_m(m-n, 64);
        sbMatrix zero_m(m-n, 64);
        bool_init_false(pIdx, zero_m);
        sbMatrix one_share_m(m-n, 1);
        bool_init_true(pIdx, one_share_m);
        bool_init_i64(pIdx, std::numeric_limits<i64>::max(), null_share_m, enc, eval, runtime);
        
        for(size_t i=0; i<4; i++){
            std::memcpy(A[i].mShares[0].data() + n, null_share_m.mShares[0].data(), (m-n) * sizeof(null_share_m.mShares[0](0, 0)));
            std::memcpy(A[i].mShares[1].data() + n, null_share_m.mShares[1].data(), (m-n) * sizeof(null_share_m.mShares[1](0, 0)));
        }
        std::memcpy(fx_sorted_auged.mShares[0].data() + n, zero_m.mShares[0].data(), (m-n) * sizeof(zero_m.mShares[0](0, 0)));
        std::memcpy(fx_sorted_auged.mShares[1].data() + n, zero_m.mShares[1].data(), (m-n) * sizeof(zero_m.mShares[1](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[0].data() + n, one_share_m.mShares[0].data(), (m-n) * sizeof(one_share_m.mShares[0](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[1].data() + n, one_share_m.mShares[1].data(), (m-n) * sizeof(one_share_m.mShares[1](0, 0)));
    }
    else{
        for(size_t i=0; i<4; i++){
            std::memcpy(A[i].mShares[0].data(), T_sorted[i].mShares[0].data(), m * sizeof(T_sorted[i].mShares[0](0, 0)));
            std::memcpy(A[i].mShares[1].data(), T_sorted[i].mShares[1].data(), m * sizeof(T_sorted[i].mShares[1](0, 0)));
        }
        std::memcpy(fx_sorted_auged.mShares[0].data(), fx_prime.mShares[0].data(), m * sizeof(fx_prime.mShares[0](0, 0)));
        std::memcpy(fx_sorted_auged.mShares[1].data(), fx_prime.mShares[1].data(), m * sizeof(fx_prime.mShares[1](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[0].data(), flag_prime.mShares[0].data(), m * sizeof(flag_prime.mShares[0](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[1].data(), flag_prime.mShares[1].data(), m * sizeof(flag_prime.mShares[1](0, 0)));
    }


    //step 3: distribute loop
    // 计算 j = 2^⌈log₂(m)⌉ - 1
    i64 step_size;
    i64 power = static_cast<i64>(std::ceil(std::log2(m)))-1;
    step_size = 1ULL << power;

    while(step_size >= 1){
        for(i64 i = m-1-step_size; i >= 0; i--){
            sbMatrix flag_i(1, 1), fx_i(1, 64);
            sbMatrix flag_i_plus_j(1, 1), fx_i_plus_j(1, 64);
            fx_i.mShares[0](0, 0) = fx_sorted_auged.mShares[0](i, 0);
            fx_i.mShares[1](0, 0) = fx_sorted_auged.mShares[1](i, 0);
            flag_i.mShares[0](0, 0) = flag_sorted_auged.mShares[0](i, 0);
            flag_i.mShares[1](0, 0) = flag_sorted_auged.mShares[1](i, 0);

            fx_i_plus_j.mShares[0](0, 0) = fx_sorted_auged.mShares[0](i+step_size, 0);
            fx_i_plus_j.mShares[1](0, 0) = fx_sorted_auged.mShares[1](i+step_size, 0);
            flag_i_plus_j.mShares[0](0, 0) = flag_sorted_auged.mShares[0](i+step_size, 0);
            flag_i_plus_j.mShares[1](0, 0) = flag_sorted_auged.mShares[1](i+step_size, 0);

            sbMatrix i_plus_j(1, 64);
            bool_init_i64(pIdx, i+step_size-1, i_plus_j, enc, eval, runtime);

            //目标idx超过i+j则交换，否则值不变
            sbMatrix cond1(1, 1);
            bool_cipher_lt(pIdx, i_plus_j, fx_i, cond1, enc, eval, runtime);

            sbMatrix A_k_i(1, 64), A_k_i_plus_j(1, 64);
            sbMatrix A_k_i_new(1, 64), A_k_i_plus_new(1, 64);
            for(size_t k=0; k<4; k++){
                A_k_i.mShares[0](0, 0) = A[k].mShares[0](i, 0);
                A_k_i.mShares[1](0, 0) = A[k].mShares[1](i, 0);
                A_k_i_plus_j.mShares[0](0, 0) = A[k].mShares[0](i+step_size, 0);
                A_k_i_plus_j.mShares[1](0, 0) = A[k].mShares[1](i+step_size, 0);

                bool_cipher_selector(pIdx, cond1, A_k_i_plus_j, A_k_i, A_k_i_new, enc, eval, runtime);
            
                A[k].mShares[0](i, 0) = A_k_i_new.mShares[0](0, 0);
                A[k].mShares[1](i, 0) = A_k_i_new.mShares[1](0, 0);

                bool_cipher_selector(pIdx, cond1, A_k_i, A_k_i_plus_j, A_k_i_plus_new, enc, eval, runtime);

                A[k].mShares[0](i+step_size, 0) = A_k_i_plus_new.mShares[0](0, 0);
                A[k].mShares[1](i+step_size, 0) = A_k_i_plus_new.mShares[1](0, 0);

            }

            sbMatrix flag_i_new(1, 1), fx_i_new(1, 64);
            bool_cipher_selector(pIdx, cond1, fx_i_plus_j, fx_i, fx_i_new, enc, eval, runtime);
            bool_cipher_selector(pIdx, cond1, flag_i_plus_j, flag_i, flag_i_new, enc, eval, runtime);
            
            fx_sorted_auged.mShares[0](i, 0) = fx_i_new.mShares[0](0, 0);
            fx_sorted_auged.mShares[1](i, 0) = fx_i_new.mShares[1](0, 0);
            flag_sorted_auged.mShares[0](i, 0) = flag_i_new.mShares[0](0, 0);
            flag_sorted_auged.mShares[1](i, 0) = flag_i_new.mShares[1](0, 0);

            sbMatrix flag_i_plus_new(1, 1), fx_i_plus_new(1, 64);
            bool_cipher_selector(pIdx, cond1, fx_i, fx_i_plus_j, fx_i_plus_new, enc, eval, runtime);
            bool_cipher_selector(pIdx, cond1, flag_i, flag_i_plus_j, flag_i_plus_new, enc, eval, runtime);
         
            fx_sorted_auged.mShares[0](i+step_size, 0) = fx_i_plus_new.mShares[0](0, 0);
            fx_sorted_auged.mShares[1](i+step_size, 0) = fx_i_plus_new.mShares[1](0, 0);
            flag_sorted_auged.mShares[0](i+step_size, 0) = flag_i_plus_new.mShares[0](0, 0);
            flag_sorted_auged.mShares[1](i+step_size, 0) = flag_i_plus_new.mShares[1](0, 0);
        }
        step_size = step_size / 2;
    }

    
    // //DEBUG
    // for(size_t i=0;i<4;i++){
    //     i64Matrix A_plain(m, 1);
    //     enc.revealAll(runtime, A[i], A_plain).get();
    //     if(pIdx == 0){
    //         std::cout << "A[" << i << "]_plain: " << std::endl;
    //         for(size_t j=0; j<m; j++){
    //             std::cout << A_plain(j, 0) << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }
    // //-----correct

    return ;
}

void align_table(int pIdx, std::vector<sbMatrix> &T, std::vector<sbMatrix> &T_aligned,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int m=T[0].rows();
    sbMatrix q(m, 64), k(m, 64);
    bool_init_false(pIdx, q);
    bool_init_false(pIdx, k);

    //step 1: set e.ii
    for(size_t i=0; i<m; i++){
        //j是否相同: same_attr
        sbMatrix j_i(1, 64);
        j_i.mShares[0](0, 0) = T[0].mShares[0](i, 0);
        j_i.mShares[1](0, 0) = T[0].mShares[1](i, 0);

        sbMatrix same_attr(1, 1), not_same_attr(1,1);
        if(i == 0){
            bool_init_false(pIdx, same_attr);
        } else {
            sbMatrix j_i_prev(1, 64);
            j_i_prev.mShares[0](0, 0) = T[0].mShares[0](i-1, 0);
            j_i_prev.mShares[1](0, 0) = T[0].mShares[1](i-1, 0);
            bool_cipher_eq(pIdx, j_i, j_i_prev, same_attr, enc, eval, runtime);
        }

        sbMatrix zero(1, 64);
        bool_init_false(pIdx, zero);
        sbMatrix one(1, 64);
        bool_init_true(pIdx, one);

        sbMatrix current_q(1, 64),current_k(1, 64);
        if(i == 0) {
            current_q = zero;
            current_k = zero;
        } else {
            current_q.mShares[0](0, 0) = q.mShares[0](i-1, 0);
            current_q.mShares[1](0, 0) = q.mShares[1](i-1, 0);
            current_k.mShares[0](0, 0) = k.mShares[0](i-1, 0);
            current_k.mShares[1](0, 0) = k.mShares[1](i-1, 0);
        }

        sbMatrix q_plus(1, 64);
        bool_cipher_add(pIdx, current_q, one, q_plus, enc, eval, runtime);
        sbMatrix k_plus(1, 64);
        bool_cipher_add(pIdx, current_k, one, k_plus, enc, eval, runtime);

        //if same :q++, k不变 else :q=0,k=0
        sbMatrix q_mid(1, 64), k_mid(1, 64);
        bool_cipher_selector(pIdx, same_attr, q_plus, zero, q_mid, enc, eval, runtime);
        bool_cipher_selector(pIdx, same_attr, current_k, zero, k_mid, enc, eval, runtime);
   
        //if  q_new>alpha_1-1 : q=0 ,k++ ; else q=0,k不变
        sbMatrix cond(1, 1);
        sbMatrix alpha_1_i(1, 64),alpha_1_minus_one(1, 64);
        alpha_1_i.mShares[0](0, 0) = T[2].mShares[0](i, 0);
        alpha_1_i.mShares[1](0, 0) = T[2].mShares[1](i, 0);
        bool_cipher_sub(pIdx, alpha_1_i, one, alpha_1_minus_one, enc, eval, runtime);
        bool_cipher_lt(pIdx, alpha_1_minus_one, q_mid, cond, enc, eval, runtime);

        sbMatrix q_new(1, 64), k_new(1, 64);
        bool_cipher_selector(pIdx, cond, zero, q_mid, q_new, enc, eval, runtime);
        bool_cipher_selector(pIdx, cond, k_plus, k_mid, k_new, enc, eval, runtime);


        q.mShares[0](i, 0) = q_new.mShares[0](0, 0);
        q.mShares[1](i, 0) = q_new.mShares[1](0, 0);
        k.mShares[0](i, 0) = k_new.mShares[0](0, 0);
        k.mShares[1](i, 0) = k_new.mShares[1](0, 0);

    }

    //这个好像必须得转成si64才能算算术乘法
    si64Matrix q_si(m, 1), k_si(m, 1);
    si64Matrix ii_si(m, 1);
    si64Matrix alpha_2_si(m, 1);
    bool2arith(pIdx, q, q_si, enc, eval, runtime);
    bool2arith(pIdx, k, k_si, enc, eval, runtime);
    bool2arith(pIdx, T[3], alpha_2_si, enc, eval, runtime);
    cipher_mul(pIdx, q_si, alpha_2_si, ii_si, eval, enc, runtime);
    ii_si = ii_si + k_si;

    // //DEBUG
    // i64Matrix q_plain(m, 1), k_plain(m, 1), ii_plain(m, 1);
    // enc.revealAll(runtime, q, q_plain).get();
    // enc.revealAll(runtime, k, k_plain).get();
    // enc.revealAll(runtime, ii_si, ii_plain).get();
    // if(pIdx == 0){
    //     std::cout<<"q:"<<std::endl;
    //     for(size_t i=0; i<m; i++){
    //         std::cout<<q_plain(i, 0)<<" ";
    //     }
    //     std::cout<<std::endl;
    //     std::cout<<"k:"<<std::endl;
    //     for(size_t i=0; i<m; i++){
    //         std::cout<<k_plain(i, 0)<<" ";
    //     }
    //     std::cout<<std::endl;
    //     std::cout<<"ii:"<<std::endl;
    //     for(size_t i=0; i<m; i++){
    //         std::cout<<ii_plain(i, 0)<<" ";
    //     }
    //     std::cout<<std::endl;    
    // }
    // //-------correct

    //step 2: sort: j,ii
    sbMatrix  ii_sb(m, 64);
    arith2bool(pIdx, ii_si, ii_sb, enc, eval, runtime);

    sbMatrix entry_key(m, 64);
    for(size_t i=0; i< m; i++){
        entry_key.mShares[0](i, 0) = (T[0].mShares[0](i, 0) & 0x1FFFFFF) << 39| (ii_sb.mShares[0](i, 0) & 0x3FFF) << 25 | (T[1].mShares[0](i, 0)&0x1FFFFFF);
        entry_key.mShares[1](i, 0) = (T[0].mShares[1](i, 0) & 0x1FFFFFF) << 39| (ii_sb.mShares[1](i, 0) & 0x3FFF) << 25 | (T[1].mShares[1](i, 0)&0x1FFFFFF);
    }

    sbMatrix entry_key_sorted(m, 64);
    odd_even_merge_sort(entry_key, entry_key_sorted, pIdx, enc, eval, runtime);

    T_aligned.resize(2);
    T_aligned[0].resize(m, 64);
    T_aligned[1].resize(m, 64);

    for(size_t i=0; i< m; i++){
        T_aligned[0].mShares[0](i, 0) = (entry_key_sorted.mShares[0](i, 0) >> 39) & 0x1FFFFFF;
        T_aligned[0].mShares[1](i, 0) = (entry_key_sorted.mShares[1](i, 0) >> 39) & 0x1FFFFFF;
        
        T_aligned[1].mShares[0](i, 0) = entry_key_sorted.mShares[0](i, 0)& 0x1FFFFFF;
        T_aligned[1].mShares[1](i, 0) = entry_key_sorted.mShares[1](i, 0)& 0x1FFFFFF;
    }

    return ;
    
}

void join(int pIdx, std::vector<si64Matrix> &T_1, std::vector<si64Matrix> &T_2, std::vector<si64Matrix> &T_joined,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    std::vector<sbMatrix> T_1_auged(4),T_2_auged(4);
    augment_table(pIdx, T_1, T_2, T_1_auged, T_2_auged, enc, eval, runtime);

    std::vector<sbMatrix> T_1_expanded(4);
    i64Matrix output_size(1, 1);
    oblivious_expand(pIdx, T_1_auged,  T_1_expanded, 0, output_size, enc, eval, runtime);
    std::vector<sbMatrix> T_2_expanded(4);
    oblivious_expand(pIdx, T_2_auged,  T_2_expanded, 1, output_size, enc, eval, runtime);
    i64 m = output_size(0, 0);

    std::vector<sbMatrix> T_2_aligned(2);
    align_table(pIdx,T_2_expanded, T_2_aligned, enc, eval, runtime);

    T_joined.resize(3);
    for(size_t i=0; i<3; i++){
        T_joined[i].resize(m, 1);
    }
    bool2arith(pIdx, T_1_expanded[0], T_joined[0], enc, eval, runtime);
    bool2arith(pIdx, T_1_expanded[1], T_joined[1], enc, eval, runtime);
    bool2arith(pIdx, T_2_aligned[1], T_joined[2], enc, eval, runtime);

    return ;
}

