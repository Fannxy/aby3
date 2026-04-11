#include "feddb.h"
#include "utils.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3-GORAM-Core/Sort.h"
#include "../aby3-GORAM-Core/Basics.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>


using namespace oc;
using namespace aby3;



void both2cipher(int pIdx, std::vector<int> &plain_cols_idx, std::vector<int> &cipher_cols_idx, 
    std::vector<i64Matrix> &input_plain_cols, std::vector<si64Matrix> &input_cipher_cols,
    std::vector<si64Matrix> &input_cols,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    size_t input_plain_cols_size = plain_cols_idx.size();
    size_t input_cipher_cols_size = cipher_cols_idx.size();
    
    size_t input_size = input_plain_cols_size + input_cipher_cols_size;
    input_cols.resize(input_size);
        
    
    for(size_t i=0 ; i < input_cipher_cols_size ; i++){
        input_cols[cipher_cols_idx[i]] = input_cipher_cols[i];
    }
    
    size_t input_unit_size ;
    if(input_plain_cols_size > 0){
        input_unit_size= input_plain_cols[0].rows();
    }else if (input_cipher_cols_size > 0){
        input_unit_size= input_cipher_cols[0].rows();
    }else{
        throw std::runtime_error("INPUT is empty");
    }
    
  
    i64Matrix input_plain_concat(input_unit_size*input_plain_cols_size, 1);
    si64Matrix input_concat(input_unit_size*input_plain_cols_size, 1);
    for(size_t i=0; i<input_plain_cols_size; i++){
        std::memcpy(input_plain_concat.data() + i*input_unit_size, input_plain_cols[i].data(), input_unit_size * sizeof(i64));
    }

    if (pIdx == 0){
        enc.localIntMatrix(runtime, input_plain_concat, input_concat).get();
    } else {
        enc.remoteIntMatrix(runtime, input_concat).get();
    }

    for (size_t s=0; s<2; s++){
        for (size_t i=0 ; i < input_plain_cols_size ; i++){
            input_cols[plain_cols_idx[i]].mShares[s] = input_concat.mShares[s].block(i*input_unit_size, 0, input_unit_size, 1);
        }
    }

    // for (size_t i=0 ; i < input_plain_cols_size ; i++){
    //     si64Matrix input_i(input_unit_size, 1);
    //     if (pIdx == 0){
    //         enc.localIntMatrix(runtime, input_plain_cols[i], input_i).get();
    //     } else {
    //         enc.remoteIntMatrix(runtime, input_i).get();
    //     }
    //     input_cols[plain_cols_idx[i]] = input_i;
    // }

    
    return ;
}

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
//TODO：memcpy
    // for(size_t i=0; i<len; i++){
    //     Tres[i].resize(unit_len, 1);
    //     Tres[i]=prev_perm[i];
    // }
    Tres = prev_perm;

    return ;

    }


void persist_plain(int pIdx, const std::string &table_name, std::vector<i64Matrix> &T){
    std::string filename = "./aby3-Feddb-tmpfile/plain/" + table_name + "_" + std::to_string(pIdx) + ".bin";

    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error: Unable to open file " << filename << " for writing" << std::endl;
        return;
    }

    int cols = T.size();
    int rows = T[0].rows();
    fwrite(&cols, sizeof(int), 1, fp);
    fwrite(&rows, sizeof(int), 1, fp);

    for(size_t i = 0; i < T.size(); i++){
        fwrite(T[i].data(), sizeof(T[i](0, 0)), rows, fp);
    }

    fclose(fp);
    return;
}

void read_plain(int pIdx, const std::string &table_name, std::vector<i64Matrix> &T){
    std::string filename = "./aby3-Feddb-tmpfile/plain/" + table_name + "_" + std::to_string(pIdx) + ".bin";

    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        std::cerr << "Error: Unable to open file " << filename << " for reading" << std::endl;
        return;
    }

    int cols, rows;
    fread(&cols, sizeof(int), 1, fp);
    fread(&rows, sizeof(int), 1, fp);
    T.resize(cols);

    for(size_t i = 0; i < T.size(); i++){
        T[i].resize(rows, 1);
        fread(T[i].data(), sizeof(T[i](0, 0)), rows, fp);
    }

    fclose(fp);
    return;
}

void oblivious_idx_select(int pIdx, si64Matrix &v, si64Matrix &idx, si64Matrix &result,              
        Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    size_t v_len = v.rows();
    size_t idx_len = idx.rows();

    //step zero: 将v和idx拼接成t
    size_t t_len=v_len+idx_len;
    si64Matrix t(t_len, 1);

    switch(pIdx){
        case 0:
            std::fill_n(t.mShares[0].data(), v_len, (i64)0);
            for(size_t i=0; i<v_len; i++) t.mShares[1](i, 0) = i;
            break;
        case 1:
            std::fill_n(t.mShares[0].data(), v_len, (i64)0);
            std::fill_n(t.mShares[1].data(), v_len, (i64)0);
            break;
        case 2:
            for(size_t i=0; i<v_len; i++) t.mShares[0](i, 0) = i;
            std::fill_n(t.mShares[1].data(), v_len, (i64)0);
            break;
    }
   
    std::memcpy(t.mShares[0].data() + v_len, idx.mShares[0].data(), idx_len * sizeof(idx.mShares[0](0, 0)));
    std::memcpy(t.mShares[1].data() + v_len, idx.mShares[1].data(), idx_len * sizeof(idx.mShares[1](0, 0)));



    //step one: 对t进行genperm得到sigma
    si64Matrix sigma(t_len, 1);
    genPerm(pIdx, t, sigma, enc, eval, runtime);    

    //step two: prefixsum_{-1}(v) && u
    si64Matrix prefix_inv(v_len, v.cols());
    prefixsum_inv(pIdx, v, prefix_inv);
    si64Matrix u(t_len, v.cols());
    std::memcpy(u.mShares[0].data(), prefix_inv.mShares[0].data(), v_len * sizeof(prefix_inv.mShares[0](0, 0)));
    std::memcpy(u.mShares[1].data(), prefix_inv.mShares[1].data(), v_len * sizeof(prefix_inv.mShares[1](0, 0)));

    std::fill_n(u.mShares[0].data() + v_len, idx_len, 0);
    std::fill_n(u.mShares[1].data() + v_len, idx_len, 0);

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
    
        // Step 4: 根据rsigma对u1进行置换得到u_prime = rsigma (u1)= sigma * pi^{-1} * pi([[u1]]) = sigma(u1)
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
    std::memcpy(result.mShares[0].data(), u_prime_.mShares[0].data() + v_len, idx_len * sizeof(u_prime_.mShares[0](0, 0)));
    std::memcpy(result.mShares[1].data(), u_prime_.mShares[1].data() + v_len, idx_len * sizeof(u_prime_.mShares[1](0, 0)));


    return; 
}

void oblivious_idx_select(int pIdx, std::vector<si64Matrix> &v_vec, si64Matrix &idx, std::vector<si64Matrix> &results,
        Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    size_t num_cols = v_vec.size();
    if(num_cols == 0) return;

    size_t v_len = v_vec[0].rows();
    size_t idx_len = idx.rows();
    size_t t_len = v_len + idx_len;

    //step zero: 构造t（和单列版本一样，t只取决于v_len和idx，与v的内容无关）
    si64Matrix t(t_len, 1);
    switch(pIdx){
        case 0:
            std::fill_n(t.mShares[0].data(), v_len, (i64)0);
            for(size_t i=0; i<v_len; i++) t.mShares[1](i, 0) = i;
            break;
        case 1:
            std::fill_n(t.mShares[0].data(), v_len, (i64)0);
            std::fill_n(t.mShares[1].data(), v_len, (i64)0);
            break;
        case 2:
            for(size_t i=0; i<v_len; i++) t.mShares[0](i, 0) = i;
            std::fill_n(t.mShares[1].data(), v_len, (i64)0);
            break;
    }
    std::memcpy(t.mShares[0].data() + v_len, idx.mShares[0].data(), idx_len * sizeof(idx.mShares[0](0, 0)));
    std::memcpy(t.mShares[1].data() + v_len, idx.mShares[1].data(), idx_len * sizeof(idx.mShares[1](0, 0)));

    //step one: 对t进行genperm得到sigma（只需一次）
    si64Matrix sigma(t_len, 1);
    genPerm(pIdx, t, sigma, enc, eval, runtime);

    //生成pi（公共随机置换，只需一次）
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = sigma.rows();

    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);
    std::vector<size_t> prev_inverse_permutation;
    std::vector<size_t> next_inverse_permutation;
    get_inverse_permutation(prev_permutation, prev_inverse_permutation);
    get_inverse_permutation(next_permutation, next_inverse_permutation);

    //计算rsigma = pi(sigma)（只需一次）
    si64Matrix rsigma(len, 1);
    si64Matrix next_perm_s(len, 1), prev_perm_s(len, 1);
    prev_perm_s = sigma;
    for(size_t i = 0; i < 3; i++){
        if(pIdx == (int)i){
            permutate(pIdx, prev_perm_s, next_perm_s, prev_permutation);
        }
        else if(pIdx == ((i+2)%3)){
            permutate(pIdx, prev_perm_s, next_perm_s, next_permutation);
        }
        else if(pIdx == (int)((i+1)%3)){
            std::fill_n(next_perm_s.mShares[0].data(), len, 0);
            std::fill_n(next_perm_s.mShares[1].data(), len, 0);
        }
        prev_perm_s = reshare_matrix(pIdx, next_perm_s, (i+1)%3, enc, runtime);
    }
    rsigma = prev_perm_s;

    // 恢复rsigma至明文，并得到rsigma^{-1}（只需一次）
    i64Matrix rsigma_plain(len, 1);
    enc.revealAll(runtime, rsigma, rsigma_plain).get();
    i64Matrix rsigma_inverse_plain(len, 1);
    permutation_inverse(rsigma_plain, rsigma_inverse_plain);

    //step two: 对每一列v_vec[col]进行 prefixsum_inv -> 构造u_vec
    std::vector<si64Matrix> u_vec(num_cols);
    for(size_t col = 0; col < num_cols; col++){
        si64Matrix prefix_inv(v_len, 1);
        prefixsum_inv(pIdx, v_vec[col], prefix_inv);
        u_vec[col].resize(t_len, 1);
        std::memcpy(u_vec[col].mShares[0].data(), prefix_inv.mShares[0].data(), v_len * sizeof(i64));
        std::memcpy(u_vec[col].mShares[1].data(), prefix_inv.mShares[1].data(), v_len * sizeof(i64));
        std::fill_n(u_vec[col].mShares[0].data() + v_len, idx_len, 0);
        std::fill_n(u_vec[col].mShares[1].data() + v_len, idx_len, 0);
    }

    //step three: [[u1_vec]] = pi([[u_vec]])  对所有列一起做
    std::vector<si64Matrix> u1_vec(num_cols);
    std::vector<si64Matrix> prev_perm_vec(num_cols), next_perm_vec(num_cols);
    for(size_t col = 0; col < num_cols; col++){
        prev_perm_vec[col].resize(len, 1);
        next_perm_vec[col].resize(len, 1);
        prev_perm_vec[col] = u_vec[col];
    }
    for(size_t i = 0; i < 3; i++){
        if(pIdx == (int)i){
            permutate(pIdx, prev_perm_vec, next_perm_vec, prev_permutation);
        }
        else if(pIdx == ((i+2)%3)){
            permutate(pIdx, prev_perm_vec, next_perm_vec, next_permutation);
        }
        else if(pIdx == (int)((i+1)%3)){
            for(size_t col = 0; col < num_cols; col++){
                std::fill_n(next_perm_vec[col].mShares[0].data(), len, 0);
                std::fill_n(next_perm_vec[col].mShares[1].data(), len, 0);
            }
        }
        for(size_t col = 0; col < num_cols; col++){
            prev_perm_vec[col] = reshare_matrix(pIdx, next_perm_vec[col], (i+1)%3, enc, runtime);
        }
    }
    u1_vec = prev_perm_vec;

    // step four: u_prime = rsigma(u1), 对所有列一起做; 再分别计算 prefixsum
    std::vector<si64Matrix> u_prime_vec(num_cols);
    for(size_t col = 0; col < num_cols; col++){
        u_prime_vec[col].resize(len, 1);
        permutate(pIdx, u1_vec[col], u_prime_vec[col], rsigma_plain);
    }

    std::vector<si64Matrix> u_prime_prefix_vec(num_cols);
    for(size_t col = 0; col < num_cols; col++){
        u_prime_prefix_vec[col].resize(len, 1);
        prefixsum(pIdx, u_prime_vec[col], u_prime_prefix_vec[col]);
    }

    //step five: unapply permutation
    // [[u1_]] = rsigma^{-1}(u_prime_prefix), 对所有列一起做
    std::vector<si64Matrix> u_1_vec(num_cols);
    for(size_t col = 0; col < num_cols; col++){
        u_1_vec[col].resize(len, 1);
        permutate(pIdx, u_prime_prefix_vec[col], u_1_vec[col], rsigma_inverse_plain);
    }

    // [[u']] = pi^{-1}([[u1_]]), 对所有列一起做
    for(size_t col = 0; col < num_cols; col++){
        prev_perm_vec[col] = u_1_vec[col];
    }
    for(int id = 2; id >= 0; id--){
        if(pIdx == id){
            permutate(pIdx, prev_perm_vec, next_perm_vec, prev_inverse_permutation);
        }
        else if(pIdx == ((id+2)%3)){
            permutate(pIdx, prev_perm_vec, next_perm_vec, next_inverse_permutation);
        }
        else if(pIdx == (id+1)%3){
            for(size_t col = 0; col < num_cols; col++){
                std::fill_n(next_perm_vec[col].mShares[0].data(), len, 0);
                std::fill_n(next_perm_vec[col].mShares[1].data(), len, 0);
            }
        }
        for(size_t col = 0; col < num_cols; col++){
            prev_perm_vec[col] = reshare_matrix(pIdx, next_perm_vec[col], (id+1)%3, enc, runtime);
        }
    }

    //step six: get results
    results.resize(num_cols);
    for(size_t col = 0; col < num_cols; col++){
        results[col].resize(idx_len, 1);
        std::memcpy(results[col].mShares[0].data(), prev_perm_vec[col].mShares[0].data() + v_len, idx_len * sizeof(i64));
        std::memcpy(results[col].mShares[1].data(), prev_perm_vec[col].mShares[1].data() + v_len, idx_len * sizeof(i64));
    }

    return;
}

//agg: sum(valcol) group by keycol : only val:non-zero entries
//data_key为多列group_key_cols， 一列data_val
void index_agg(int pIdx, si64Matrix &equalFlag, i64Matrix &idx,std::vector<si64Matrix> &data_key,si64Matrix &data_val,std::vector<si64Matrix> &finalRes,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime ){
    
        int rows=data_key[0].rows();

        std::vector<si64Matrix> keycols(data_key.size());
        si64Matrix valcol(rows, 1);

        for(size_t i=0; i<data_key.size(); i++){
            keycols[i].resize(rows, 1);
        }
        valcol.resize(rows, 1);

        //step 1: sort according to idx
        for(int i=0; i<rows ;i++){
            i64 idx_i=idx(i, 0);
            for(size_t j=0; j<data_key.size(); j++){
                keycols[j].mShares[0](i, 0) = data_key[j].mShares[0](idx_i, 0);
                keycols[j].mShares[1](i, 0) = data_key[j].mShares[1](idx_i, 0);
            }
            //keycols[j](i,0) = data_key[j](idx_i, 0);
            valcol.mShares[0](i, 0) = data_val.mShares[0](idx_i, 0);
            valcol.mShares[1](i, 0) = data_val.mShares[1](idx_i, 0);
            //valcol(i,0) = data_val(idx_i, 0);
        }

        //step 2: valcol sum 
        //noteq
        si64Matrix not_equalFlag(rows, 1), oneShared(rows, 1);
        set_const_share(pIdx, 1, oneShared, enc, eval, runtime);
        not_equalFlag = oneShared - equalFlag;

        // [Original serial scan — O(n) cipher_mul rounds, commented out]
        // for(int i=0; i<(rows-1) ;i++){
        //     si64Matrix leftval(1,1),rightval(1,1);
        //     leftval.mShares[0](0,0)=valcol.mShares[0](i,0);
        //     leftval.mShares[1](0,0)=valcol.mShares[1](i,0);
        //     rightval.mShares[0](0,0)=valcol.mShares[0](i+1,0);
        //     rightval.mShares[1](0,0)=valcol.mShares[1](i+1,0);
        //     si64Matrix new_leftval(1,1),new_rightval(1,1);
        //     si64Matrix not_eqFlag(1,1),eqFlag(1,1);
        //     eqFlag.mShares[0](0,0)=equalFlag.mShares[0](i+1,0);
        //     eqFlag.mShares[1](0,0)=equalFlag.mShares[1](i+1,0);
        //     not_eqFlag.mShares[0](0,0) = not_equalFlag.mShares[0](i+1,0);
        //     not_eqFlag.mShares[1](0,0) = not_equalFlag.mShares[1](i+1,0);
        //     si64Matrix mul_left(2, 1), mul_right(2, 1), mul_result(2, 1);
        //     mul_left.mShares[0](0,0) = leftval.mShares[0](0,0);
        //     mul_left.mShares[1](0,0) = leftval.mShares[1](0,0);
        //     mul_left.mShares[0](1,0) = leftval.mShares[0](0,0);
        //     mul_left.mShares[1](1,0) = leftval.mShares[1](0,0);
        //     mul_right.mShares[0](0,0) = not_eqFlag.mShares[0](0,0);
        //     mul_right.mShares[1](0,0) = not_eqFlag.mShares[1](0,0);
        //     mul_right.mShares[0](1,0) = eqFlag.mShares[0](0,0);
        //     mul_right.mShares[1](1,0) = eqFlag.mShares[1](0,0);
        //     cipher_mul(pIdx, mul_left, mul_right, mul_result, eval, enc, runtime);
        //     new_leftval.mShares[0](0,0) = mul_result.mShares[0](0,0);
        //     new_leftval.mShares[1](0,0) = mul_result.mShares[1](0,0);
        //     si64Matrix leftval_times_eqFlag(1,1);
        //     leftval_times_eqFlag.mShares[0](0,0) = mul_result.mShares[0](1,0);
        //     leftval_times_eqFlag.mShares[1](0,0) = mul_result.mShares[1](1,0);
        //     new_rightval = rightval + leftval_times_eqFlag;
        //     valcol.mShares[0](i,0) = new_leftval.mShares[0](0,0);
        //     valcol.mShares[1](i,0) = new_leftval.mShares[1](0,0);
        //     valcol.mShares[0](i+1,0) = new_rightval.mShares[0](0,0);
        //     valcol.mShares[1](i+1,0) = new_rightval.mShares[1](0,0);
        // }

        // [Optimized: segmented parallel prefix sum — O(log n) cipher_mul rounds]
        // Monoid: (val, flag) ⊗ (val', flag') = (val' + val * flag', flag * flag')
        // scan_val accumulates the prefix sum; scan_eq propagates segment flags.
        si64Matrix scan_val(rows, 1), scan_eq(rows, 1);
        // Initialize scan_val = valcol
        std::memcpy(scan_val.mShares[0].data(), valcol.mShares[0].data(), rows * sizeof(i64));
        std::memcpy(scan_val.mShares[1].data(), valcol.mShares[1].data(), rows * sizeof(i64));
        scan_eq = equalFlag;

        // Up-sweep: parallel prefix scan with stride doubling
        for (int d = 1; d < rows; d <<= 1) {
            int count = rows - d;  // number of active elements this level

            // Layout: rows [0..count-1]       = scan_val[i-d] * scan_eq[i]   (value accumulation)
            //         rows [count..2*count-1]  = scan_eq[i-d]  * scan_eq[i]   (flag propagation)
            si64Matrix mul_left(2 * count, 1), mul_right(2 * count, 1), mul_result(2 * count, 1);

            for (int s = 0; s < 2; s++) {
                std::memcpy(mul_left.mShares[s].data(),
                            scan_val.mShares[s].data(), count * sizeof(i64));
                std::memcpy(mul_right.mShares[s].data(),
                            scan_eq.mShares[s].data() + d, count * sizeof(i64));
                std::memcpy(mul_left.mShares[s].data() + count,
                            scan_eq.mShares[s].data(), count * sizeof(i64));
                std::memcpy(mul_right.mShares[s].data() + count,
                            scan_eq.mShares[s].data() + d, count * sizeof(i64));
            }

            cipher_mul(pIdx, mul_left, mul_right, mul_result, eval, enc, runtime);

            // Write back results
            for (int s = 0; s < 2; s++) {
                // scan_val[d..rows-1] += mul_result[0..count-1]
                for (int k = 0; k < count; k++)
                    scan_val.mShares[s](d + k, 0) += mul_result.mShares[s](k, 0);
                // scan_eq[d..rows-1] = mul_result[count..2*count-1]
                std::memcpy(scan_eq.mShares[s].data() + d,
                            mul_result.mShares[s].data() + count, count * sizeof(i64));
            }
        }

 
        // Now zero out non-segment-tail positions: valcol[i] = scan_val[i] * not_eq[i+1]
        // The last row is always a segment tail: valcol[rows-1] = scan_val[rows-1]
        if (rows > 1) {
            int n = rows - 1;
            si64Matrix clear_left(n, 1), clear_right(n, 1), clear_result(n, 1);
            for (int s = 0; s < 2; s++) {
                std::memcpy(clear_left.mShares[s].data(),
                            scan_val.mShares[s].data(), n * sizeof(i64));
                std::memcpy(clear_right.mShares[s].data(),
                            not_equalFlag.mShares[s].data() + 1, n * sizeof(i64));
            }
            cipher_mul(pIdx, clear_left, clear_right, clear_result, eval, enc, runtime);
            for (int s = 0; s < 2; s++) {
                std::memcpy(valcol.mShares[s].data(),
                            clear_result.mShares[s].data(), n * sizeof(i64));
            }
        }
        valcol.mShares[0](rows - 1, 0) = scan_val.mShares[0](rows - 1, 0);
        valcol.mShares[1](rows - 1, 0) = scan_val.mShares[1](rows - 1, 0);

        //只留下非零行
        sbMatrix zeroFlag_sb(rows,1);
        si64Matrix zeroShared(rows,1);
        set_const_share(pIdx, 0, zeroShared, enc, eval, runtime);

        cipher_eq(pIdx, valcol, zeroShared, zeroFlag_sb, eval, runtime);
        
        si64Matrix zeroFlag(rows,1);
        bool2arith(pIdx, zeroFlag_sb, zeroFlag, enc, eval, runtime);


        std::vector<si64Matrix> res(data_key.size()+2), shuffledRes(data_key.size()+2);
        for(size_t i=0; i<data_key.size(); i++){
            res[i]=keycols[i];
        }
        res[data_key.size()]=valcol;
        res[data_key.size()+1]=zeroFlag;

        shuffle(pIdx, res, shuffledRes, enc, eval, runtime);
        zeroFlag=shuffledRes[data_key.size()+1];
        i64Matrix zeroFlag_plain(rows,1);
        enc.revealAll(runtime, zeroFlag, zeroFlag_plain).get();

        // Build index of non-zero rows once
        std::vector<int> keepIdx;
        keepIdx.reserve(rows);
        for (int i = 0; i < rows; i++) {
            if (zeroFlag_plain(i, 0) == 0)
                keepIdx.push_back(i);
        }
        int resRows = (int)keepIdx.size();

        // Compact all columns (keys + val) uniformly
        size_t totalCols = data_key.size() + 1;  // key cols + val col
        finalRes.resize(totalCols);
        for (size_t j = 0; j < totalCols; j++) {
            finalRes[j].resize(resRows, 1);
            for (int s = 0; s < 2; s++) {
                for (int r = 0; r < resRows; r++)
                    finalRes[j].mShares[s](r, 0) = shuffledRes[j].mShares[s](keepIdx[r], 0);
            }
        }


        return;

    }



void group_by_common(int pIdx, si64Matrix &key, si64Matrix &val, 
    si64Matrix &key_g, si64Matrix &val_g, si64Matrix &e, si64Matrix &perm_GN, si64Matrix &key_GN, si64Matrix &key_out,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int rows = key.rows();
    int key_cols = key.cols();
    //std::cout<<"rows: "<<rows<<", cols: "<<key_cols<<std::endl;
    si64Matrix null_share(rows, key_cols),one_share(rows, 1);
    set_const_share(pIdx, std::numeric_limits<i64>::max(), null_share, enc, eval, runtime);
    set_const_share(pIdx, 1, one_share, enc, eval, runtime);
    
    //step 1
    si64Matrix perm(rows, 1);
    si64Matrix k_v_concat(rows, key_cols + 1);
    //倒序
    k_v_concat.mShares[0].col(0) = val.mShares[0].col(0);
    k_v_concat.mShares[1].col(0) = val.mShares[1].col(0);
    for(size_t j=1; j<=key_cols; j++){
        k_v_concat.mShares[0].col(j) = key.mShares[0].col(key_cols-j);
        k_v_concat.mShares[1].col(j) = key.mShares[1].col(key_cols-j);
    }
    genPerm(pIdx, k_v_concat, perm, enc, eval, runtime);


    //step 2
    key_g.resize(rows, key_cols);
    val_g.resize(rows, 1);
    applyPerm(pIdx, perm, key, key_g, enc, eval, runtime);
    applyPerm(pIdx, perm, val, val_g, enc, eval, runtime);

    //step 3
    sbMatrix f(rows-1, 64);
    sbMatrix f_vector((rows-1)*key_cols, 1);
    //std::cout<<"1"<<std::endl;
    compare_consecutive_rows_arith(pIdx, key_g, f_vector, enc, eval, runtime);
    
    // //DEBUG
    // i64Matrix f_vector_plain((rows-1)*key_cols, 1);
    // enc.revealAll(runtime, f_vector, f_vector_plain).get();
    // if(pIdx == 0){
    //     std::cout<<"f_vector_plain: "<<std::endl;
    //     for(int i=0;i<rows-1;i++){
    //         for(int j=0;j<key_cols;j++){
    //             std::cout<<f_vector_plain(i*key_cols+j, 0)<<" ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }
    // //--correct

 
    // 初始化 f 为全1
    f.resize(rows-1, 1);
    bool_init_true(pIdx, f);
    
    // 对每一列，将 f 与 f_vector 中对应的值进行 AND
    for(int j=0;j<key_cols;j++){
        sbMatrix f_vector_col(rows-1, 1);
        // 从 f_vector 中提取第 j 列的值（即每行的第 j 个元素）
        for(int i=0;i<rows-1;i++){
            f_vector_col.mShares[0](i, 0) = f_vector.mShares[0](i*key_cols+j, 0);
            f_vector_col.mShares[1](i, 0) = f_vector.mShares[1](i*key_cols+j, 0);
        }
        // 使用安全的 AND 电路
        sbMatrix f_new(rows-1, 1);
        bool_cipher_and(pIdx, f, f_vector_col, f_new, enc, eval, runtime);
        f = f_new;
    }

    si64Matrix f_si(rows-1, 1);
    bool2arith(pIdx, f, f_si, enc, eval, runtime);

    // //DEBUG
    // i64Matrix f_plain(rows-1, 1);
    // enc.revealAll(runtime, f, f_plain).get();
    // if(pIdx == 0){
    //     std::cout<<"f_plain: "<<std::endl;
    //     for(int i=0;i<rows-1;i++){
    //         std::cout<<f_plain(i, 0)<<" ";
    //     }
    // }
    // //
    // i64Matrix f_si_plain(rows-1, 1);
    // enc.revealAll(runtime, f_si, f_si_plain).get();
    // if(pIdx == 0){
    //     std::cout<<"f_si_plain: "<<std::endl;
    //     for(int i=0;i<rows-1;i++){
    //         std::cout<<f_si_plain(i, 0)<<" ";
    //     }
    // }

    //step 4
    si64Matrix e_partial(rows-1, 1);
    e_partial = one_share - f_si;
    e.resize(rows, 1);
    std::memcpy(e.mShares[0].data(), e_partial.mShares[0].data(), (rows-1) * sizeof(e_partial.mShares[0](0, 0)));
    std::memcpy(e.mShares[1].data(), e_partial.mShares[1].data(), (rows-1) * sizeof(e_partial.mShares[1](0, 0)));

    e.mShares[0](rows-1, 0) = one_share.mShares[0](rows-1, 0) ;
    e.mShares[1](rows-1, 0) = one_share.mShares[1](rows-1, 0) ;

    // //DEBUG
    // i64Matrix e_plain(rows, 1);
    // enc.revealAll(runtime, e, e_plain).get();
    // if(pIdx == 0){
    //     std::cout<<"e_plain: "<<std::endl;
    //     for(int i=0;i<rows;i++){
    //         std::cout<<e_plain(i, 0)<<" ";
    //     }
    // }
    // //
    

    //step 5
    si64Matrix tmp(rows, key_cols);
    si64Matrix e_ex(rows, key_cols);
    arith_cols_expand(e, e_ex);
    cipher_mul(pIdx, e_ex, key_g-null_share, tmp, eval, enc, runtime);

    // //DEBUG
    // i64Matrix tmp_plain(rows, key_cols);
    // enc.revealAll(runtime, tmp, tmp_plain).get();
    // if(pIdx == 0){
    //     std::cout<<"tmp_plain: "<<std::endl;
    //     for(int i=0;i<rows;i++){
    //         for(int j=0;j<key_cols;j++){
    //             std::cout<<tmp_plain(i, j)<<" ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }
    // //
    key_GN.resize(rows, key_cols);
    key_GN = tmp + null_share;

    //step 6
    tmp = one_share-e;
    genPerm(pIdx, tmp, perm_GN, enc, eval, runtime);

    //step 7
    key_out.resize(rows, key_cols);
    applyPerm(pIdx, perm_GN, key_GN, key_out, enc, eval, runtime);

    return;
}


void group_by_common_without_val(int pIdx, si64Matrix &key, 
    si64Matrix &key_g, si64Matrix &e, si64Matrix &perm_GN, si64Matrix &key_GN, si64Matrix &key_out,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int rows = key.rows();
    int cols = key.cols();

    
    si64Matrix null_share(rows, cols),one_share(rows, 1);
    set_const_share(pIdx, std::numeric_limits<i64>::max(), null_share, enc, eval, runtime);
    set_const_share(pIdx, 1, one_share, enc, eval, runtime);

    key_g.resize(rows, cols);
    key_GN.resize(rows, cols);
    key_out.resize(rows, cols);
    
    //step 1
    si64Matrix perm(rows, 1);
    genPerm(pIdx, key, perm, enc, eval, runtime);

    //step 2
    applyPerm(pIdx, perm, key, key_g, enc, eval, runtime);

    //step 3
    sbMatrix f(rows-1, 64);
    sbMatrix f_vector((rows-1)*cols, 1);
    compare_consecutive_rows_arith(pIdx, key_g, f_vector, enc, eval, runtime);
     //eq_result(rows*cols, 1))

    f.resize(rows-1, 1);
    bool_init_true(pIdx, f);
    
    // 对每一列，将 f 与 f_vector 中对应的值进行 AND
    for(int j=0;j<cols;j++){
        sbMatrix f_vector_col(rows-1, 1);
        // 从 f_vector 中提取第 j 列的值（即每行的第 j 个元素）
        for(int i=0;i<rows-1;i++){
            f_vector_col.mShares[0](i, 0) = f_vector.mShares[0](i*cols+j, 0);
            f_vector_col.mShares[1](i, 0) = f_vector.mShares[1](i*cols+j, 0);
        }
        // 使用安全的 AND 电路
        sbMatrix f_new(rows-1, 1);
        bool_cipher_and(pIdx, f, f_vector_col, f_new, enc, eval, runtime);
        f = f_new;
    }
 
    si64Matrix f_si(rows-1, 1);
    bool2arith(pIdx, f, f_si, enc, eval, runtime);

    //step 4
    si64Matrix e_partial(rows-1, 1);
    e_partial = one_share - f_si;
    e.resize(rows, 1);
    std::memcpy(e.mShares[0].data(), e_partial.mShares[0].data(), (rows-1) * sizeof(e_partial.mShares[0](0, 0)));
    std::memcpy(e.mShares[1].data(), e_partial.mShares[1].data(), (rows-1) * sizeof(e_partial.mShares[1](0, 0)));

    e.mShares[0](rows-1, 0) = one_share.mShares[0](rows-1, 0) ;
    e.mShares[1](rows-1, 0) = one_share.mShares[1](rows-1, 0) ;
    

    //step 5
    si64Matrix tmp(rows, cols);
    si64Matrix e_ex(rows, cols);
    arith_cols_expand(e, e_ex);
    cipher_mul(pIdx, e_ex, key_g-null_share, tmp, eval, enc, runtime);
    key_GN = tmp + null_share;

    //step 6
    tmp = one_share-e;
    genPerm(pIdx, tmp, perm_GN, enc, eval, runtime);

    //step 7
    applyPerm(pIdx, perm_GN, key_GN, key_out, enc, eval, runtime);

    return;
}


void group_count(int pIdx, std::vector<si64Matrix> &key, 
    std::vector<si64Matrix> &key_out, si64Matrix &c, 
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

    int rows = key[0].rows();
    int cols = key.size();
    key_out.resize(cols);
    // 初始化 key_out 中每个矩阵的大小
    for(int j=0; j<cols; j++){
        key_out[j].resize(rows, 1);
    }
    c.resize(rows, 1);

    si64Matrix key_concat(rows, cols);
    //step 0:concat key（逆序）
    for(int j=0;j<cols;j++){
        key_concat.mShares[0].col(j) = key[cols-1-j].mShares[0].col(0);
        key_concat.mShares[1].col(j) = key[cols-1-j].mShares[1].col(0);
    }

    //step 1
    si64Matrix key_g(rows, cols);
    si64Matrix key_GN(rows, cols);
    si64Matrix key_out_concat(rows, cols);

    si64Matrix e(rows, 1);
    si64Matrix perm_GN(rows, 1);
    group_by_common_without_val(pIdx, key_concat, key_g, e, perm_GN, key_GN, key_out_concat, enc, eval, runtime);

    for(int j=0;j<cols;j++){
        key_out[cols-1-j].mShares[0].col(0) = key_out_concat.mShares[0].col(j);
        key_out[cols-1-j].mShares[1].col(0) = key_out_concat.mShares[1].col(j);
    }

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


void group_sum(int pIdx, std::vector<si64Matrix> &key, si64Matrix &val, 
    std::vector<si64Matrix> &key_out, si64Matrix &sum,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){ 
    
    int rows = key[0].rows();
    int cols = key.size();
    key_out.resize(cols);
    // 初始化 key_out 中每个矩阵的大小
    for(int j=0; j<cols; j++){
        key_out[j].resize(rows, 1);
    }
    sum.resize(rows, 1);

    si64Matrix key_concat(rows, cols);
    //step 0:concat key（逆序）
    for(int j=0;j<cols;j++){
        key_concat.mShares[0].col(j) = key[cols-1-j].mShares[0].col(0);
        key_concat.mShares[1].col(j) = key[cols-1-j].mShares[1].col(0);
    }

    si64Matrix key_g(rows, cols);
    si64Matrix key_GN(rows, cols);
    si64Matrix key_out_concat(rows, cols);

    si64Matrix val_g(rows, 1);
    si64Matrix e(rows, 1);
    si64Matrix perm_GN(rows, 1);
    //step 1
    group_by_common(pIdx, key_concat, val, key_g, val_g, e, perm_GN, key_GN, key_out_concat, enc, eval, runtime);

    for(int j=0;j<cols;j++){
        key_out[cols-1-j].mShares[0].col(0) = key_out_concat.mShares[0].col(j);
        key_out[cols-1-j].mShares[1].col(0) = key_out_concat.mShares[1].col(j);
    }


    //step 2
    si64Matrix w(rows, 1);
    prefixsum(pIdx, val_g, w);

    //step 3
    si64Matrix w_m(rows, 1);
    si64Matrix x(rows,1);
    std::fill_n(w_m.mShares[0].data(), rows, w.mShares[0](rows-1, 0));
    std::fill_n(w_m.mShares[1].data(), rows, w.mShares[1](rows-1, 0));
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

void group_max(int pIdx, std::vector<si64Matrix> &key, si64Matrix &val, 
    std::vector<si64Matrix> &key_out, si64Matrix &max,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int rows = key[0].rows();
    int cols = key.size();
    key_out.resize(cols);
    // 初始化 key_out 中每个矩阵的大小
    for(int j=0; j<cols; j++){
        key_out[j].resize(rows, 1);
    }
    max.resize(rows, 1);

    si64Matrix key_concat(rows, cols);
    //step 0:concat key（逆序）
    for(int j=0;j<cols;j++){
        key_concat.mShares[0].col(j) = key[cols-1-j].mShares[0].col(0);
        key_concat.mShares[1].col(j) = key[cols-1-j].mShares[1].col(0);
    }


    //step 1
    si64Matrix key_g(rows, cols);
    si64Matrix key_GN(rows, cols);
    si64Matrix key_out_concat(rows, cols);

    si64Matrix val_g(rows, 1);
    si64Matrix e(rows, 1);
    si64Matrix perm_GN(rows, 1);
    group_by_common(pIdx, key_concat, val, key_g, val_g, e, perm_GN, key_GN, key_out_concat, enc, eval, runtime);

    for(int j=0;j<cols;j++){
        key_out[cols-1-j].mShares[0].col(0) = key_out_concat.mShares[0].col(j);
        key_out[cols-1-j].mShares[1].col(0) = key_out_concat.mShares[1].col(j);
    }

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

void group_min(int pIdx, std::vector<si64Matrix> &key, si64Matrix &val, 
    std::vector<si64Matrix> &key_out, si64Matrix &min,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int rows = key[0].rows();
    int cols = key.size();
    key_out.resize(cols);
    // 初始化 key_out 中每个矩阵的大小
    for(int j=0; j<cols; j++){
        key_out[j].resize(rows, 1);
    }
    min.resize(rows, 1);

    si64Matrix key_concat(rows, cols);
    //step 0:concat key（逆序）
    for(int j=0;j<cols;j++){
        key_concat.mShares[0].col(j) = key[cols-1-j].mShares[0].col(0);
        key_concat.mShares[1].col(j) = key[cols-1-j].mShares[1].col(0);
    }

    si64Matrix key_g(rows, cols);
    si64Matrix key_GN(rows, cols);
    si64Matrix key_out_concat(rows, cols);

    si64Matrix val_g(rows, 1);
    si64Matrix e(rows, 1);
    si64Matrix perm_GN(rows, 1);
    //step 1
    group_by_common(pIdx, key_concat, val, key_g, val_g, e, perm_GN, key_GN, key_out_concat, enc, eval, runtime);

    for(int j=0;j<cols;j++){
        key_out[cols-1-j].mShares[0].col(0) = key_out_concat.mShares[0].col(j);
        key_out[cols-1-j].mShares[1].col(0) = key_out_concat.mShares[1].col(j);
    }

    //step 2.1
    si64Matrix e_prime(rows, 1);
    std::memcpy(e_prime.mShares[0].data() + 1, e.mShares[0].data(), (rows-1) * sizeof(i64));
    std::memcpy(e_prime.mShares[1].data() + 1, e.mShares[1].data(), (rows-1) * sizeof(i64));
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
//optimized:sbMatrix + permutation sort + parallel scan
void augment_table(int pIdx, std::vector<si64Matrix> &T_1_key, std::vector<si64Matrix> &T_1_other, 
    std::vector<si64Matrix> &T_2_key, std::vector<si64Matrix> &T_2_other,
    std::vector<si64Matrix> &T_1_auged, std::vector<si64Matrix> &T_2_auged,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    // // 时间测量变量（按需保留，当前禁用）
    // auto total_start = std::chrono::high_resolution_clock::now();
    // std::vector<std::pair<std::string, double>> timing_results;

    int key_num = T_1_key.size();
    int other_1_num = T_1_other.size();
    int other_2_num = T_2_other.size();
    int max_other_num = std::max(other_1_num, other_2_num);
    
    size_t len1 = T_1_key[0].rows();
    size_t len2 = T_2_key[0].rows();

    //step 1:concatenate T_1 and T_2 & sort
    // auto t1 = std::chrono::high_resolution_clock::now();
        //step 1.1: tid tag
    si64Matrix tid_0(len1, 1),tid_1(len2, 1);
    set_const_share(pIdx, 0, tid_0, enc, eval, runtime);
    set_const_share(pIdx, 1, tid_1, enc, eval, runtime);

        //step 1.2: entry_key concat
    int entry_key_cols = key_num+max_other_num+1;
    si64Matrix entry_key_si(len1+len2, entry_key_cols);

    //other列（max_other_num列）
    int other_area_cols = entry_key_cols - key_num - 1;

    for(int s = 0; s < 2; ++s) {
        // 1. 批量清零
        entry_key_si.mShares[s].leftCols(other_area_cols).setZero();

        // 2. 填充 T_1_other (占前 len1 行)
        for(int i = 0; i < other_1_num; ++i) {
            int target_col = (entry_key_cols - key_num - 2) - i;
            if (target_col < 0) break;
            entry_key_si.mShares[s].block(0, target_col, len1, 1) = T_1_other[i].mShares[s];
        }

        // 3. 填充 T_2_other (占接下来的 len2 行)
        for(int i = 0; i < other_2_num; ++i) {
            int target_col = (entry_key_cols - key_num - 2) - i;
            if (target_col < 0) break;
            entry_key_si.mShares[s].block(len1, target_col, len2, 1) = T_2_other[i].mShares[s];
        }
    }

    //tid列
    for(int s = 0; s < 2; ++s) {
        entry_key_si.mShares[s].block(0, entry_key_cols-key_num-1, len1, 1) = tid_0.mShares[s];
        entry_key_si.mShares[s].block(len1, entry_key_cols-key_num-1, len2, 1) = tid_1.mShares[s];
    }

    //key列
    for(int i = 0; i < key_num; ++i) {
        int target_col = entry_key_cols - 1 - i;
        for(int s = 0; s < 2; ++s) {
            //T_1_key[i]变成简单的vector
            entry_key_si.mShares[s].block(0, target_col, len1, 1) = T_1_key[i].mShares[s];
            entry_key_si.mShares[s].block(len1, target_col, len2, 1) = T_2_key[i].mShares[s];
        }
    }
    
    si64Matrix entry_key_perm(len1+len2, key_num+1);
    entry_key_perm.mShares[0] = entry_key_si.mShares[0].block(0, entry_key_cols-key_num-1, len1+len2, key_num + 1);
    entry_key_perm.mShares[1] = entry_key_si.mShares[1].block(0, entry_key_cols-key_num-1, len1+len2, key_num + 1);
    // auto t_step1_assign_end = std::chrono::high_resolution_clock::now();
    // double step1_assign_time = std::chrono::duration<double, std::milli>(t_step1_assign_end - t1).count();
    // timing_results.push_back({"augment_step1.1_assign_before_sort", step1_assign_time});

    si64Matrix perm(len1+len2, 1);
    auto t_genperm_begin = std::chrono::high_resolution_clock::now();
    //genPerm(pIdx, entry_key_si, perm, enc, eval, runtime);
    genPerm(pIdx, entry_key_perm, perm, enc, eval, runtime);
    // auto t_genperm_end = std::chrono::high_resolution_clock::now();
    // double step1_genperm_time = std::chrono::duration<double, std::milli>(t_genperm_end - t1).count();
    // timing_results.push_back({"augment_step1.2_genperm", step1_genperm_time});

    si64Matrix entry_key_sorted_si(len1+len2, entry_key_cols);
    //auto t_applyperm_begin = std::chrono::high_resolution_clock::now();
    applyPerm(pIdx, perm, entry_key_si, entry_key_sorted_si, enc, eval, runtime);
    // auto t_applyperm_end = std::chrono::high_resolution_clock::now();
    // double step1_applyperm_time = std::chrono::duration<double, std::milli>(t_applyperm_end - t1).count();
    // timing_results.push_back({"augment_step1.3_applyperm", step1_applyperm_time});

    int entry_key_bitcount = 64*key_num+64*max_other_num+64;
    sbMatrix entry_key_sorted(len1+len2, entry_key_bitcount);
    arith2bool(pIdx, entry_key_sorted_si, entry_key_sorted, enc, eval, runtime);

    sbMatrix j_sb(len1+len2, 64*key_num);
    sbMatrix tid_sb(len1+len2, 64);
    sbMatrix d_sb(len1+len2, 64*max_other_num);

    for(size_t s = 0; s < 2; ++s){
        size_t rows = len1+len2;
        size_t src_stride = entry_key_sorted.mShares[s].cols();
        for(size_t i=0; i< rows; i++){
            i64* src_row = entry_key_sorted.mShares[s].data() + i * src_stride;
            //d: columns [0, max_other_num)
            std::memcpy(d_sb.mShares[s].data() + i * max_other_num, src_row, max_other_num * sizeof(i64));
            //tid: column [max_other_num]
            tid_sb.mShares[s].data()[i] = src_row[max_other_num];
            //j: columns [max_other_num+1, max_other_num+1+key_num)
            std::memcpy(j_sb.mShares[s].data() + i * key_num, src_row + max_other_num + 1, key_num * sizeof(i64));
        }
    }
    // auto t2 = std::chrono::high_resolution_clock::now();
    // double step1_sort_result_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"augment_step1_perm_sort_result", step1_sort_result_time});

    //step 2: full-dimension
        //step 2.1: downward scan
    // t1 = std::chrono::high_resolution_clock::now();
    // auto t_same_attr_begin = t1;
    sbMatrix alpha_1(len1+len2, 64), alpha_2(len1+len2, 64);
        //same_attr & not_same_attr:
    sbMatrix false_matrix(1,1);
    bool_init_false(pIdx, false_matrix);
    sbMatrix same_attr_partial(len1+len2-1, 1);
    compare_consecutive_rows_bool(pIdx, j_sb, same_attr_partial, enc, eval, runtime);
 
    
    sbMatrix same_attr(len1+len2, 1);
    same_attr.mShares[0](0, 0) = false_matrix.mShares[0](0, 0);
    same_attr.mShares[1](0, 0) = false_matrix.mShares[1](0, 0);
    std::memcpy(same_attr.mShares[0].data() + 1, same_attr_partial.mShares[0].data(), (len1+len2-1) * sizeof(same_attr_partial.mShares[0](0, 0)));
    std::memcpy(same_attr.mShares[1].data() + 1, same_attr_partial.mShares[1].data(), (len1+len2-1) * sizeof(same_attr_partial.mShares[1](0, 0)));
        //1bit-> 64bits
    same_attr.resize(len1+len2, 64);

    for(size_t i=0; i<len1+len2; i++){
        same_attr.mShares[0](i, 0) = same_attr.mShares[0](i, 0) == 1 ? -1 : -0;
        same_attr.mShares[1](i, 0) = same_attr.mShares[1](i, 0) == 1 ? -1 : -0;
    }
        //not_same_attr:
    sbMatrix not_same_attr(len1+len2, 64);
    bool_cipher_not(pIdx, same_attr, not_same_attr);

    // 记录 same_attr 相关时间（比较 + 按列 AND 聚合）
    // auto t_same_attr_end = std::chrono::high_resolution_clock::now();
    // double step2_same_attr_time = std::chrono::duration<double, std::milli>(t_same_attr_end - t_same_attr_begin).count();
    // timing_results.push_back({"augment_step2.1_same_attr", step2_same_attr_time});
   

        //is_table_1 & is_table_2:
    //auto t_is_table_begin = std::chrono::high_resolution_clock::now();
    sbMatrix zero64_vector(len1+len2, 64);
    bool_init_false(pIdx, zero64_vector);
    sbMatrix is_table_1(len1+len2, 1);
    bool_cipher_eq(pIdx, tid_sb, zero64_vector, is_table_1, enc, eval, runtime);
    is_table_1.resize(len1+len2, 64);

    for(size_t i=0; i<len1+len2; i++){
        is_table_1.mShares[0](i, 0) = is_table_1.mShares[0](i, 0) == 1 ? -1 : -0;
        is_table_1.mShares[1](i, 0) = is_table_1.mShares[1](i, 0) == 1 ? -1 : -0;
    }
    sbMatrix is_table_2(len1+len2, 64);
    bool_cipher_not(pIdx, is_table_1, is_table_2);
    // auto t_is_table_end = std::chrono::high_resolution_clock::now();
    // double step2_is_table_time = std::chrono::duration<double, std::milli>(t_is_table_end - t_is_table_begin).count();
    // timing_results.push_back({"augment_step2.1_is_table", step2_is_table_time});

        //condition
    sbMatrix is_table_vector(4*(len1+len2), 64), same_attr_vector(4*(len1+len2), 64);
    std::memcpy(is_table_vector.mShares[0].data(), is_table_1.mShares[0].data(), (len1+len2) * sizeof(is_table_1.mShares[0](0, 0)));
    std::memcpy(is_table_vector.mShares[1].data(), is_table_1.mShares[1].data(), (len1+len2) * sizeof(is_table_1.mShares[1](0, 0)));
    std::memcpy(same_attr_vector.mShares[0].data(), not_same_attr.mShares[0].data(), (len1+len2) * sizeof(not_same_attr.mShares[0](0, 0)));
    std::memcpy(same_attr_vector.mShares[1].data(), not_same_attr.mShares[1].data(), (len1+len2) * sizeof(not_same_attr.mShares[1](0, 0)));
    std::memcpy(is_table_vector.mShares[0].data() + (len1+len2), is_table_1.mShares[0].data(), (len1+len2) * sizeof(is_table_1.mShares[0](0, 0)));
    std::memcpy(is_table_vector.mShares[1].data() + (len1+len2), is_table_1.mShares[1].data(), (len1+len2) * sizeof(is_table_1.mShares[1](0, 0)));
    std::memcpy(same_attr_vector.mShares[0].data() + (len1+len2), same_attr.mShares[0].data(), (len1+len2) * sizeof(same_attr.mShares[0](0, 0)));
    std::memcpy(same_attr_vector.mShares[1].data() + (len1+len2), same_attr.mShares[1].data(), (len1+len2) * sizeof(same_attr.mShares[1](0, 0)));
    std::memcpy(is_table_vector.mShares[0].data() + 2*(len1+len2), is_table_2.mShares[0].data(), (len1+len2) * sizeof(is_table_2.mShares[0](0, 0)));
    std::memcpy(is_table_vector.mShares[1].data() + 2*(len1+len2), is_table_2.mShares[1].data(), (len1+len2) * sizeof(is_table_2.mShares[1](0, 0)));
    std::memcpy(same_attr_vector.mShares[0].data() + 2*(len1+len2), not_same_attr.mShares[0].data(), (len1+len2) * sizeof(not_same_attr.mShares[0](0, 0)));
    std::memcpy(same_attr_vector.mShares[1].data() + 2*(len1+len2), not_same_attr.mShares[1].data(), (len1+len2) * sizeof(not_same_attr.mShares[1](0, 0)));
    std::memcpy(is_table_vector.mShares[0].data() + 3*(len1+len2), is_table_2.mShares[0].data(), (len1+len2) * sizeof(is_table_2.mShares[0](0, 0)));
    std::memcpy(is_table_vector.mShares[1].data() + 3*(len1+len2), is_table_2.mShares[1].data(), (len1+len2) * sizeof(is_table_2.mShares[1](0, 0)));
    std::memcpy(same_attr_vector.mShares[0].data() + 3*(len1+len2), same_attr.mShares[0].data(), (len1+len2) * sizeof(same_attr.mShares[0](0, 0)));
    std::memcpy(same_attr_vector.mShares[1].data() + 3*(len1+len2), same_attr.mShares[1].data(), (len1+len2) * sizeof(same_attr.mShares[1](0, 0)));
    sbMatrix condition(4*(len1+len2), 64);
    bool_cipher_and(pIdx, is_table_vector, same_attr_vector, condition, enc, eval, runtime);

        //for loop: update alpha1 & alpha2
    
    sbMatrix zero64(2, 64);
    sbMatrix one64(2, 64);
    bool_init_false(pIdx, zero64);
    bool_init_true(pIdx, one64);

    //auto t_downward_loop_begin = std::chrono::high_resolution_clock::now();

    //op3:parallel scan
    int N = 2*(len1 + len2);
    sbMatrix A(N, 64), B_pre(N, 64), B(N, 64);
    sbMatrix newA(N, 64), newB(N, 64), tmp(N, 64);

    for(int i=0; i<len1+len2; i++){
       //alpha1
       A.mShares[0](i, 0) = condition.mShares[0](i+(len1+len2), 0) ^ condition.mShares[0](i+3*(len1+len2), 0);
       A.mShares[1](i, 0) = condition.mShares[1](i+(len1+len2), 0) ^ condition.mShares[1](i+3*(len1+len2), 0);

       B_pre.mShares[0](i, 0) = condition.mShares[0](i, 0) ^ condition.mShares[0](i+1*(len1+len2), 0);
       B_pre.mShares[1](i, 0) = condition.mShares[1](i, 0) ^ condition.mShares[1](i+1*(len1+len2), 0);

       //alpha2
       A.mShares[0](i+len1+len2, 0) = condition.mShares[0](i+3*(len1+len2), 0);
       A.mShares[1](i+len1+len2, 0) = condition.mShares[1](i+3*(len1+len2), 0);

       B_pre.mShares[0](i+len1+len2, 0) = condition.mShares[0](i+3*(len1+len2), 0) ^ condition.mShares[0](i+2*(len1+len2), 0);
       B_pre.mShares[1](i+len1+len2, 0) = condition.mShares[1](i+3*(len1+len2), 0) ^ condition.mShares[1](i+2*(len1+len2), 0);
    }
    
    sbMatrix oneN(N, 64);
    bool_init_true(pIdx, oneN);
    bool_cipher_and(pIdx, B_pre, oneN, B, enc, eval, runtime);

    sbMatrix A_B_shift(2*N, 64),A_concat(2*N, 64);
    sbMatrix newA_tmp_concat(2*N, 64);

    int half = len1 + len2;
    for(int d = 1 ; d < half ; d<<=1){
       
        //初始单位元：A_shift=1(AND单位元)，B_shift=0(ADD/XOR单位元)
        switch (pIdx) {
            case 0:
                std::fill_n(A_B_shift.mShares[0].data(), N, 0);
                std::fill_n(A_B_shift.mShares[0].data() + N, N, 1);
                std::fill_n(A_B_shift.mShares[1].data(), N, 0);
                std::fill_n(A_B_shift.mShares[1].data() + N, N, 0);
                break;
            case 1:
                std::fill_n(A_B_shift.mShares[0].data(), N, -1);
                std::fill_n(A_B_shift.mShares[0].data() + N, N, 1);
                std::fill_n(A_B_shift.mShares[1].data(), N, 0);
                std::fill_n(A_B_shift.mShares[1].data() + N, N, 1);
                break;
            case 2:
                std::fill_n(A_B_shift.mShares[0].data(), N, 0);
                std::fill_n(A_B_shift.mShares[0].data() + N, N, 0);
                std::fill_n(A_B_shift.mShares[1].data(), N, -1);
                std::fill_n(A_B_shift.mShares[1].data() + N, N, 1);
                break;
            default:
                break;
        }
        
        
        int cnt = half - d; // elements to shift-copy per half
        if(cnt > 0){
            for(size_t s = 0; s < 2; ++s){
                // first half: A[0..cnt) -> A_B_shift[d..half), B[0..cnt) -> A_B_shift[N+d..N+half)
                std::memcpy(A_B_shift.mShares[s].data() + d,          A.mShares[s].data(),        cnt * sizeof(i64));
                std::memcpy(A_B_shift.mShares[s].data() + N + d,      B.mShares[s].data(),        cnt * sizeof(i64));
                // second half: A[half..half+cnt) -> A_B_shift[half+d..N), B[half..half+cnt) -> A_B_shift[N+half+d..2N)
                std::memcpy(A_B_shift.mShares[s].data() + half + d,   A.mShares[s].data() + half, cnt * sizeof(i64));
                std::memcpy(A_B_shift.mShares[s].data() + N+half + d, B.mShares[s].data() + half, cnt * sizeof(i64));
            }
        }
        

        std::memcpy(A_concat.mShares[0].data(), A.mShares[0].data(), N * sizeof(A.mShares[0](0, 0)));
        std::memcpy(A_concat.mShares[0].data() + N, A.mShares[0].data(), N * sizeof(A.mShares[0](0, 0)));
        std::memcpy(A_concat.mShares[1].data(), A.mShares[1].data(), N * sizeof(A.mShares[1](0, 0)));
        std::memcpy(A_concat.mShares[1].data() + N, A.mShares[1].data(), N * sizeof(A.mShares[1](0, 0)));
        

        bool_cipher_and(pIdx, A_concat, A_B_shift, newA_tmp_concat, enc, eval, runtime);
        std::memcpy(newA.mShares[0].data(), newA_tmp_concat.mShares[0].data(), N * sizeof(newA_tmp_concat.mShares[0](0, 0)));
        std::memcpy(newA.mShares[1].data(), newA_tmp_concat.mShares[1].data(), N * sizeof(newA_tmp_concat.mShares[1](0, 0)));
        std::memcpy(tmp.mShares[0].data(), newA_tmp_concat.mShares[0].data() + N, N * sizeof(newA_tmp_concat.mShares[0](0, 0)));
        std::memcpy(tmp.mShares[1].data(), newA_tmp_concat.mShares[1].data() + N, N * sizeof(newA_tmp_concat.mShares[1](0, 0)));
        bool_cipher_add(pIdx, tmp, B, newB, enc, eval, runtime);

        A = newA;
        B = newB;
    }   

    
    for(size_t s = 0; s < 2; ++s){
        std::memcpy(alpha_1.mShares[s].data(), B.mShares[s].data(), half * sizeof(i64));
        std::memcpy(alpha_2.mShares[s].data(), B.mShares[s].data() + half, half * sizeof(i64));
    }


    // auto t_downward_loop_end = std::chrono::high_resolution_clock::now();
    // double step2_downward_loop_time = std::chrono::duration<double, std::milli>(t_downward_loop_end - t_downward_loop_begin).count();
    // timing_results.push_back({"augment_step2.1_downward_loop", step2_downward_loop_time});

    // t2 = std::chrono::high_resolution_clock::now();
    // double step2_downward_scan_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"augment_step2.1_downward_scan", step2_downward_scan_time});

        //step 2.2: upward scan
    //t1 = std::chrono::high_resolution_clock::now(); 
        //same_inv_attr & not_same_inv_attr:
    sbMatrix same_inv_attr(len1+len2, 1);
    std::memcpy(same_inv_attr.mShares[0].data() , same_attr_partial.mShares[0].data(), (len1+len2-1) * sizeof(same_attr_partial.mShares[0](0, 0)));
    std::memcpy(same_inv_attr.mShares[1].data() , same_attr_partial.mShares[1].data(), (len1+len2-1) * sizeof(same_attr_partial.mShares[1](0, 0)));
    same_inv_attr.mShares[0](len1+len2-1, 0) = false_matrix.mShares[0](0, 0);
    same_inv_attr.mShares[1](len1+len2-1, 0) = false_matrix.mShares[1](0, 0);
    same_inv_attr.resize(len1+len2, 64);

    for(size_t i=0; i<len1+len2; i++){
        same_inv_attr.mShares[0](i, 0) = same_inv_attr.mShares[0](i, 0) == 1 ? -1 : -0;
        same_inv_attr.mShares[1](i, 0) = same_inv_attr.mShares[1](i, 0) == 1 ? -1 : -0;
    }
    sbMatrix not_same_inv_attr(len1+len2, 64);
    bool_cipher_not(pIdx, same_inv_attr, not_same_inv_attr);

        //condition
    sbMatrix same_inv_attr_vector(4*(len1+len2), 64);
    std::memcpy(same_inv_attr_vector.mShares[0].data(), not_same_inv_attr.mShares[0].data(), (len1+len2) * sizeof(not_same_inv_attr.mShares[0](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[1].data(), not_same_inv_attr.mShares[1].data(), (len1+len2) * sizeof(not_same_inv_attr.mShares[1](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[0].data() + (len1+len2), same_inv_attr.mShares[0].data(), (len1+len2) * sizeof(same_inv_attr.mShares[0](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[1].data() + (len1+len2), same_inv_attr.mShares[1].data(), (len1+len2) * sizeof(same_inv_attr.mShares[1](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[0].data() + 2*(len1+len2), not_same_inv_attr.mShares[0].data(), (len1+len2) * sizeof(not_same_inv_attr.mShares[0](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[1].data() + 2*(len1+len2), not_same_inv_attr.mShares[1].data(), (len1+len2) * sizeof(not_same_inv_attr.mShares[1](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[0].data() + 3*(len1+len2), same_inv_attr.mShares[0].data(), (len1+len2) * sizeof(same_inv_attr.mShares[0](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[1].data() + 3*(len1+len2), same_inv_attr.mShares[1].data(), (len1+len2) * sizeof(same_inv_attr.mShares[1](0, 0)));

    bool_cipher_and(pIdx, is_table_vector, same_inv_attr_vector, condition, enc, eval, runtime);

    //op3:upward parallel scan

    for(int i=0; i<len1+len2; i++){
        //alpha1
        A.mShares[0](i, 0) = condition.mShares[0](i+1*(len1+len2), 0) ;
        A.mShares[1](i, 0) = condition.mShares[1](i+1*(len1+len2), 0) ;
 
        B_pre.mShares[0](i, 0) = condition.mShares[0](i, 0) ^ condition.mShares[0](i+2*(len1+len2), 0) ^ condition.mShares[0](i+3*(len1+len2), 0);
        B_pre.mShares[1](i, 0) = condition.mShares[1](i, 0) ^ condition.mShares[1](i+2*(len1+len2), 0) ^ condition.mShares[1](i+3*(len1+len2), 0);
 
        //alpha2
        A.mShares[0](i+len1+len2, 0) = condition.mShares[0](i+1*(len1+len2), 0) ^ condition.mShares[0](i+3*(len1+len2), 0);
        A.mShares[1](i+len1+len2, 0) = condition.mShares[1](i+1*(len1+len2), 0) ^ condition.mShares[1](i+3*(len1+len2), 0);
 
        B_pre.mShares[0](i+len1+len2, 0) = condition.mShares[0](i, 0) ^ condition.mShares[0](i+2*(len1+len2), 0);
        B_pre.mShares[1](i+len1+len2, 0) = condition.mShares[1](i, 0) ^ condition.mShares[1](i+2*(len1+len2), 0);
     }
     
     sbMatrix alpha_1_2_down(N, 64);
     std::memcpy(alpha_1_2_down.mShares[0].data(), alpha_1.mShares[0].data(), (len1+len2) * sizeof(alpha_1.mShares[0](0, 0)));
     std::memcpy(alpha_1_2_down.mShares[1].data(), alpha_1.mShares[1].data(), (len1+len2) * sizeof(alpha_1.mShares[1](0, 0)));
     std::memcpy(alpha_1_2_down.mShares[0].data() + (len1+len2), alpha_2.mShares[0].data(), (len1+len2) * sizeof(alpha_2.mShares[0](0, 0)));
     std::memcpy(alpha_1_2_down.mShares[1].data() + (len1+len2), alpha_2.mShares[1].data(), (len1+len2) * sizeof(alpha_2.mShares[1](0, 0)));
     bool_cipher_and(pIdx, B_pre, alpha_1_2_down, B, enc, eval, runtime);

    sbMatrix A_rev(N, 64);
    sbMatrix B_rev(N, 64);
    {
        i64 half = len1 + len2;
        for(size_t s = 0; s < 2; ++s){
            // reverse first half [0, half) of A and B
            std::reverse_copy(A.mShares[s].data(), A.mShares[s].data() + half, A_rev.mShares[s].data());
            std::reverse_copy(B.mShares[s].data(), B.mShares[s].data() + half, B_rev.mShares[s].data());
            // reverse second half [half, N) of A and B
            std::reverse_copy(A.mShares[s].data() + half, A.mShares[s].data() + N, A_rev.mShares[s].data() + half);
            std::reverse_copy(B.mShares[s].data() + half, B.mShares[s].data() + N, B_rev.mShares[s].data() + half);
        }
    }

 
     for(int d = 1 ; d < half ; d<<=1){
        
         //初始单位元：A_shift=1(AND单位元)，B_shift=0(ADD/XOR单位元)
         switch (pIdx) {
            case 0:
                std::fill_n(A_B_shift.mShares[0].data(), N, 0);
                std::fill_n(A_B_shift.mShares[0].data() + N, N, 1);
                std::fill_n(A_B_shift.mShares[1].data(), N, 0);
                std::fill_n(A_B_shift.mShares[1].data() + N, N, 0);
                break;
            case 1:
                std::fill_n(A_B_shift.mShares[0].data(), N, -1);
                std::fill_n(A_B_shift.mShares[0].data() + N, N, 1);
                std::fill_n(A_B_shift.mShares[1].data(), N, 0);
                std::fill_n(A_B_shift.mShares[1].data() + N, N, 1);
                break;
            case 2:
                std::fill_n(A_B_shift.mShares[0].data(), N, 0);
                std::fill_n(A_B_shift.mShares[0].data() + N, N, 0);
                std::fill_n(A_B_shift.mShares[1].data(), N, -1);
                std::fill_n(A_B_shift.mShares[1].data() + N, N, 1);
                break;
            default:
                break;
        }
        
            
        int cnt = half - d;
        if(cnt > 0){
            for(size_t s = 0; s < 2; ++s){
                std::memcpy(A_B_shift.mShares[s].data() + d,          A_rev.mShares[s].data(),        cnt * sizeof(i64));
                std::memcpy(A_B_shift.mShares[s].data() + N + d,      B_rev.mShares[s].data(),        cnt * sizeof(i64));
                std::memcpy(A_B_shift.mShares[s].data() + half + d,   A_rev.mShares[s].data() + half, cnt * sizeof(i64));
                std::memcpy(A_B_shift.mShares[s].data() + N+half + d, B_rev.mShares[s].data() + half, cnt * sizeof(i64));
            }
        }
        

        std::memcpy(A_concat.mShares[0].data(), A_rev.mShares[0].data(), N * sizeof(A.mShares[0](0, 0)));
        std::memcpy(A_concat.mShares[0].data() + N, A_rev.mShares[0].data(), N * sizeof(A.mShares[0](0, 0)));
        std::memcpy(A_concat.mShares[1].data(), A_rev.mShares[1].data(), N * sizeof(A.mShares[1](0, 0)));
        std::memcpy(A_concat.mShares[1].data() + N, A_rev.mShares[1].data(), N * sizeof(A.mShares[1](0, 0)));
        

        bool_cipher_and(pIdx, A_concat, A_B_shift, newA_tmp_concat, enc, eval, runtime);
        std::memcpy(newA.mShares[0].data(), newA_tmp_concat.mShares[0].data(), N * sizeof(newA_tmp_concat.mShares[0](0, 0)));
        std::memcpy(newA.mShares[1].data(), newA_tmp_concat.mShares[1].data(), N * sizeof(newA_tmp_concat.mShares[1](0, 0)));
        std::memcpy(tmp.mShares[0].data(), newA_tmp_concat.mShares[0].data() + N, N * sizeof(newA_tmp_concat.mShares[0](0, 0)));
        std::memcpy(tmp.mShares[1].data(), newA_tmp_concat.mShares[1].data() + N, N * sizeof(newA_tmp_concat.mShares[1](0, 0)));

         //bool_cipher_and(pIdx, A_rev, A_shift, newA, enc, eval, runtime);
         //bool_cipher_and(pIdx, A_rev, B_shift, tmp, enc, eval, runtime);
         bool_cipher_add(pIdx, tmp, B_rev, newB, enc, eval, runtime);
 
         A_rev = newA;
         B_rev = newB;
     }   
 
     
    for(size_t s = 0; s < 2; ++s){
        std::reverse_copy(B_rev.mShares[s].data(), B_rev.mShares[s].data() + half, alpha_1.mShares[s].data());
        std::reverse_copy(B_rev.mShares[s].data() + half, B_rev.mShares[s].data() + 2*half, alpha_2.mShares[s].data());
    }
     

    // t2 = std::chrono::high_resolution_clock::now();
    // double step2_upward_scan_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"augment_step2.2_upward_scan", step2_upward_scan_time});


    //step 3: sort by: tid
    // t1 = std::chrono::high_resolution_clock::now();
    sbMatrix entry_key(len1+len2, entry_key_bitcount+64+64);
    entry_key_sorted.resize(len1+len2,entry_key_bitcount+64+64);

    for(int i=0; i<len1+len2; i++){
        //alpha2
        entry_key.mShares[0](i, 0) = alpha_2.mShares[0](i, 0);
        entry_key.mShares[1](i, 0) = alpha_2.mShares[1](i, 0);
        //alpha1
        entry_key.mShares[0](i, 1) = alpha_1.mShares[0](i, 0);
        entry_key.mShares[1](i, 1) = alpha_1.mShares[1](i, 0);
        //d
        for(size_t j=0; j<max_other_num; j++){
            entry_key.mShares[0](i, j+2) = d_sb.mShares[0](i, j);
            entry_key.mShares[1](i, j+2) = d_sb.mShares[1](i, j);
        }
        //j
        for(size_t j=0; j<key_num; j++){
            entry_key.mShares[0](i, j+max_other_num+2) = j_sb.mShares[0](i, j);
            entry_key.mShares[1](i, j+max_other_num+2) = j_sb.mShares[1](i, j);
        }
        //tid
        entry_key.mShares[0](i, max_other_num+2+key_num) = tid_sb.mShares[0](i, 0);
        entry_key.mShares[1](i, max_other_num+2+key_num) = tid_sb.mShares[1](i, 0);

    }

    //genperm according tid,j,d
    sbMatrix entry_key_perm_sb(len1+len2, (1+max_other_num+key_num)*64);
    for(int s = 0; s < 2; ++s) {
        for(size_t i=0; i<len1+len2; i++){
            for(size_t j=2; j<3+max_other_num+key_num; j++){
                entry_key_perm_sb.mShares[s](i, j-2) = entry_key.mShares[s](i, j);
            }
        }
    }

    genPerm_bool(pIdx, entry_key_perm_sb, perm, enc, eval, runtime);
    entry_key_si.resize(len1+len2, 3 + max_other_num + key_num);
    bool2arith(pIdx, entry_key, entry_key_si, enc, eval, runtime);

    entry_key_sorted_si.resize(len1+len2, 3 + max_other_num + key_num);
    applyPerm(pIdx, perm, entry_key_si, entry_key_sorted_si, enc, eval, runtime);
    // t2 = std::chrono::high_resolution_clock::now();
    // double step3_sort_result_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"augment_step3_sort_result", step3_sort_result_time});


    //step 4: split T_1 and T_2
    // t1 = std::chrono::high_resolution_clock::now();
    //拆分直接得到std::vector<si64Matrix> &T_1_auged, std::vector<si64Matrix> &T_2_auged
    T_1_auged.resize(4);
    T_2_auged.resize(4);
    T_1_auged[0].resize(len1, key_num);
    T_1_auged[1].resize(len1, max_other_num);
    T_2_auged[0].resize(len2, key_num);
    T_2_auged[1].resize(len2, max_other_num);
    for(size_t i=2; i<4; i++){
        T_1_auged[i].resize(len1, 1);
        T_2_auged[i].resize(len2, 1);
    }


    T_1_auged[3].mShares[0].block(0, 0, len1, 1) = entry_key_sorted_si.mShares[0].block(0, 0, len1, 1);
    T_2_auged[3].mShares[0].block(0, 0, len2, 1) = entry_key_sorted_si.mShares[0].block(len1, 0, len2, 1);
    T_1_auged[2].mShares[0].block(0, 0, len1, 1) = entry_key_sorted_si.mShares[0].block(0, 1, len1, 1);
    T_2_auged[2].mShares[0].block(0, 0, len2, 1) = entry_key_sorted_si.mShares[0].block(len1, 1, len2, 1);
    T_1_auged[1].mShares[0].block(0, 0, len1, max_other_num) = entry_key_sorted_si.mShares[0].block(0, 2, len1, max_other_num);
    T_2_auged[1].mShares[0].block(0, 0, len2, max_other_num) = entry_key_sorted_si.mShares[0].block(len1, 2, len2, max_other_num);
    T_1_auged[0].mShares[0].block(0, 0, len1, key_num) = entry_key_sorted_si.mShares[0].block(0, max_other_num+2, len1, key_num);
    T_2_auged[0].mShares[0].block(0, 0, len2, key_num) = entry_key_sorted_si.mShares[0].block(len1, max_other_num+2, len2, key_num);

    T_1_auged[3].mShares[1].block(0, 0, len1, 1) = entry_key_sorted_si.mShares[1].block(0, 0, len1, 1);
    T_2_auged[3].mShares[1].block(0, 0, len2, 1) = entry_key_sorted_si.mShares[1].block(len1, 0, len2, 1);
    T_1_auged[2].mShares[1].block(0, 0, len1, 1) = entry_key_sorted_si.mShares[1].block(0, 1, len1, 1);
    T_2_auged[2].mShares[1].block(0, 0, len2, 1) = entry_key_sorted_si.mShares[1].block(len1, 1, len2, 1);
    T_1_auged[1].mShares[1].block(0, 0, len1, max_other_num) = entry_key_sorted_si.mShares[1].block(0, 2, len1, max_other_num);
    T_2_auged[1].mShares[1].block(0, 0, len2, max_other_num) = entry_key_sorted_si.mShares[1].block(len1, 2, len2, max_other_num);
    T_1_auged[0].mShares[1].block(0, 0, len1, key_num) = entry_key_sorted_si.mShares[1].block(0, max_other_num+2, len1, key_num);
    T_2_auged[0].mShares[1].block(0, 0, len2, key_num) = entry_key_sorted_si.mShares[1].block(len1, max_other_num+2, len2, key_num);
    // t2 = std::chrono::high_resolution_clock::now();
    // double step4_split_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"augment_step4_split", step4_split_time});
    // auto total_end = std::chrono::high_resolution_clock::now();
    // double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    // timing_results.push_back({"augment_total", total_time});
    // if (pIdx == 0) {
    //     std::string filename = "./join_timing_results_role0.txt";
    //     std::ofstream outFile(filename, std::ios::app);
    //     if (outFile.is_open()) {
    //         outFile << "--- augment_table Internal Timing (Role " << pIdx << ") ---" << std::endl;
    //         outFile << std::fixed << std::setprecision(3);
    //         for(const auto& result : timing_results) {
    //             outFile << "  " << result.first << ": " << result.second << " ms" << std::endl;
    //         }
    //         outFile << "----------------------------------------" << std::endl;
    //         outFile << std::endl;
    //         outFile.close();
    //     }
    // }


    return;
}


/*
void augment_table(int pIdx, std::vector<si64Matrix> &T_1_key, std::vector<si64Matrix> &T_1_other, 
    std::vector<si64Matrix> &T_2_key, std::vector<si64Matrix> &T_2_other,
    std::vector<sbMatrix> &T_1_auged, std::vector<sbMatrix> &T_2_auged,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

    // 时间测量变量
    auto total_start = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<std::string, double>> timing_results;

    int key_num = T_1_key.size();
    int other_1_num = T_1_other.size();
    int other_2_num = T_2_other.size();
    int max_other_num = std::max(other_1_num, other_2_num);
    
    size_t len1 = T_1_key[0].rows();
    size_t len2 = T_2_key[0].rows();

    //------odd-even merge sort------

    // //step 1:T_c: concatenate T_1 and T_2
    // auto t1 = std::chrono::high_resolution_clock::now();
    // //T_c[0]: key(j) ; T_c[1]: other(d)
    // std::vector<si64Matrix> T_c(2);
    // T_c[0].resize(len1+len2, key_num);
    // T_c[1].resize(len1+len2, max_other_num);
 
    // // for(size_t i=0; i<len1; i++){
    // //     for(size_t j=0; j<key_num; j++){
    // //         T_c[0].mShares[0](i, j) = T_1_key[j].mShares[0](i, 0);
    // //         T_c[0].mShares[1](i, j) = T_1_key[j].mShares[1](i, 0);
    // //     }
    // //     for(size_t j=0; j<other_1_num; j++){
    // //         T_c[1].mShares[0](i, j) = T_1_other[j].mShares[0](i, 0);
    // //         T_c[1].mShares[1](i, j) = T_1_other[j].mShares[1](i, 0);
    // //     }
    // //     if(max_other_num > other_1_num){
    // //         for(size_t j=other_1_num; j<max_other_num; j++){
    // //             T_c[1].mShares[0](i, j) = 0;
    // //             T_c[1].mShares[1](i, j) = 0;
    // //         }
    // //     }
    // // }
    // //optimize: 按列写入的,使用col()访问
    // for(size_t j=0; j<key_num; j++){
    //     // 前 len1 行：使用 head(len1)
    //     T_c[0].mShares[0].col(j).head(len1) = T_1_key[j].mShares[0].col(0);
    //     T_c[0].mShares[1].col(j).head(len1) = T_1_key[j].mShares[1].col(0);
        
    //     // 后 len2 行：使用 segment(len1, len2)
    //     T_c[0].mShares[0].col(j).segment(len1, len2) = T_2_key[j].mShares[0].col(0);
    //     T_c[0].mShares[1].col(j).segment(len1, len2) = T_2_key[j].mShares[1].col(0);
    // }
    // for(size_t j=0; j<other_1_num; j++){
    //     T_c[1].mShares[0].col(j).head(len1) = T_1_other[j].mShares[0].col(0);
    //     T_c[1].mShares[1].col(j).head(len1) = T_1_other[j].mShares[1].col(0);
    // }
    // if(max_other_num > other_1_num){
    //     for(size_t j=other_1_num; j<max_other_num; j++){
    //         // T_c[1].mShares[0](i, j) = 0;
    //         // T_c[1].mShares[1](i, j) = 0;
    //         T_c[1].mShares[0].col(j).head(len1).setZero();
    //         T_c[1].mShares[1].col(j).head(len1).setZero();
    //     }
    // }
    // for(size_t j=0; j<other_2_num; j++){
    //     T_c[1].mShares[0].col(j).segment(len1, len2) = T_2_other[j].mShares[0].col(0);
    //     T_c[1].mShares[1].col(j).segment(len1, len2) = T_2_other[j].mShares[1].col(0);
    // }
    // if(max_other_num > other_2_num){
    //     for(size_t j=other_2_num; j<max_other_num; j++){
    //         // T_c[1].mShares[0](i, j) = 0;
    //         // T_c[1].mShares[1](i, j) = 0;
    //         T_c[1].mShares[0].col(j).segment(len1, len2).setZero();
    //         T_c[1].mShares[1].col(j).segment(len1, len2).setZero();
    //     }
    // }

    // // for(size_t i=len1; i<len1+len2; i++){
    // //     for(size_t j=0; j<key_num; j++){
    // //         T_c[0].mShares[0](i, j) = T_2_key[j].mShares[0](i-len1, 0);
    // //         T_c[0].mShares[1](i, j) = T_2_key[j].mShares[1](i-len1, 0);
    // //     }
    // //     for(size_t j=0; j<other_2_num; j++){
    // //         T_c[1].mShares[0](i, j) = T_2_other[j].mShares[0](i-len1, 0);
    // //         T_c[1].mShares[1](i, j) = T_2_other[j].mShares[1](i-len1, 0);
    // //     }
    // //     if(max_other_num > other_2_num){
    // //         for(size_t j=other_2_num; j<max_other_num; j++){
    // //             T_c[1].mShares[0](i, j) = 0;
    // //             T_c[1].mShares[1](i, j) = 0;
    // //         }
    // //     }
    // // }


    
    // //T_c_bool: convert T_c to bool
    // std::vector<sbMatrix> T_c_bool(3);
    // T_c_bool[0].resize(len1+len2, key_num*64);
    // T_c_bool[1].resize(len1+len2, max_other_num*64);
    // //TODO VECTOR：由于两个的列数不一样，所以先不合并吧
    // arith2bool(pIdx, T_c[0], T_c_bool[0], enc, eval, runtime);
    // arith2bool(pIdx, T_c[1], T_c_bool[1], enc, eval, runtime);

    // //T_c_bool[2]: tid
    // T_c_bool[2].resize(len1+len2, 1);
    // sbMatrix zero(len1, 1);
    // sbMatrix one(len2, 1);
    // bool_init_false(pIdx, zero);
    // bool_init_true(pIdx, one);
    
    // std::memcpy(T_c_bool[2].mShares[0].data(), zero.mShares[0].data(), len1 * sizeof(zero.mShares[0](0, 0)));
    // std::memcpy(T_c_bool[2].mShares[1].data(), zero.mShares[1].data(), len1 * sizeof(zero.mShares[1](0, 0)));
    // std::memcpy(T_c_bool[2].mShares[0].data() + len1, one.mShares[0].data(), len2 * sizeof(one.mShares[0](0, 0)));
    // std::memcpy(T_c_bool[2].mShares[1].data() + len1, one.mShares[1].data(), len2 * sizeof(one.mShares[1](0, 0)));
    // auto t2 = std::chrono::high_resolution_clock::now();
    // double step1_arith2bool_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"augment_step1_concat_T1_T2", step1_arith2bool_time});


     
    // //step 2: sort T_c: j, tid(转成了bool进行拼接排序)
    // t1 = std::chrono::high_resolution_clock::now();
    // int entry_key_bitcount = 64*key_num+64*max_other_num+64;
    // int entry_key_cols = key_num+max_other_num+1;
    // sbMatrix entry_key(len1+len2, entry_key_bitcount);
    // for(int j=entry_key_cols-1; j>=entry_key_cols-key_num; j--){
    //     entry_key.mShares[0].col(j) = T_c_bool[0].mShares[0].col(entry_key_cols-j-1);
    // }
    // // for(int i=0; i<len1+len2; i++){
    // //     //因为是大端序，所以要逆序
    // //     //j
    // //     for(int j=entry_key_cols-1; j>=entry_key_cols-key_num; j--){
    // //         entry_key.mShares[0](i, j) = T_c_bool[0].mShares[0](i, entry_key_cols-j-1);
    // //         entry_key.mShares[1](i, j) = T_c_bool[0].mShares[1](i, entry_key_cols-j-1);
    // //     }
    // //     //tid
    // //     entry_key.mShares[0](i, entry_key_cols-key_num-1) = T_c_bool[2].mShares[0](i, 0);
    // //     entry_key.mShares[1](i, entry_key_cols-key_num-1) = T_c_bool[2].mShares[1](i, 0);
    // //     //d
    // //     for(int j=entry_key_cols-key_num-2; j>=0; j--){
    // //         entry_key.mShares[0](i, j) = T_c_bool[1].mShares[0](i, entry_key_cols-key_num-2-j);
    // //         entry_key.mShares[1](i, j) = T_c_bool[1].mShares[1](i, entry_key_cols-key_num-2-j);
    // //     }
   
    // // }


    // sbMatrix entry_key_sorted(len1+len2, entry_key_bitcount);
    // odd_even_merge_sort(entry_key, entry_key_sorted, pIdx, enc, eval, runtime);
    //------odd-even merge sort end------

    //optimize: permutation sort
    auto t1 = std::chrono::high_resolution_clock::now();
    //step 1: tid tag
    si64Matrix tid_0(len1, 1),tid_1(len2, 1);
    set_const_share(pIdx, 0, tid_0, enc, eval, runtime);
    set_const_share(pIdx, 1, tid_1, enc, eval, runtime);

    //step 2: entry_key concat
    int entry_key_cols = key_num+max_other_num+1;
    si64Matrix entry_key_si(len1+len2, entry_key_cols);
    //key列
    for(int j=entry_key_cols-1; j>=entry_key_cols-key_num; j--){
        entry_key_si.mShares[0].col(j).head(len1) = T_1_key[entry_key_cols-1-j].mShares[0].col(0);
        entry_key_si.mShares[1].col(j).head(len1) = T_1_key[entry_key_cols-1-j].mShares[1].col(0);

        entry_key_si.mShares[0].col(j).segment(len1, len2) = T_2_key[entry_key_cols-1-j].mShares[0].col(0);
        entry_key_si.mShares[1].col(j).segment(len1, len2) = T_2_key[entry_key_cols-1-j].mShares[1].col(0);
    }
    //tid列
    entry_key_si.mShares[0].col(entry_key_cols-key_num-1).head(len1) = tid_0.mShares[0].col(0);
    entry_key_si.mShares[1].col(entry_key_cols-key_num-1).head(len1) = tid_0.mShares[1].col(0);
    entry_key_si.mShares[0].col(entry_key_cols-key_num-1).segment(len1, len2) = tid_1.mShares[0].col(0);
    entry_key_si.mShares[1].col(entry_key_cols-key_num-1).segment(len1, len2) = tid_1.mShares[1].col(0);
    //other列（max_other_num列）
    for(int j=entry_key_cols-key_num-2; j>=0; j--){
        if(entry_key_cols-key_num-2-j >= other_1_num){
            entry_key_si.mShares[0].col(j).head(len1).setZero();
            entry_key_si.mShares[1].col(j).head(len1).setZero();
        }else{
            entry_key_si.mShares[0].col(j).head(len1) = T_1_other[entry_key_cols-key_num-2-j].mShares[0].col(0);
            entry_key_si.mShares[1].col(j).head(len1) = T_1_other[entry_key_cols-key_num-2-j].mShares[1].col(0);
        }

        if(entry_key_cols-key_num-2-j >= other_2_num){
            entry_key_si.mShares[0].col(j).segment(len1, len2).setZero();
            entry_key_si.mShares[1].col(j).segment(len1, len2).setZero();
        }else{
            entry_key_si.mShares[0].col(j).segment(len1, len2) = T_2_other[entry_key_cols-key_num-2-j].mShares[0].col(0);
            entry_key_si.mShares[1].col(j).segment(len1, len2) = T_2_other[entry_key_cols-key_num-2-j].mShares[1].col(0);
        }
    }

    //step 3: sort
    si64Matrix perm(len1+len2, 1);
    genPerm(pIdx, entry_key_si, perm, enc, eval, runtime);
    si64Matrix entry_key_sorted_si(len1+len2, entry_key_cols);
    applyPerm(pIdx, perm, entry_key_si, entry_key_sorted_si, enc, eval, runtime);

    int entry_key_bitcount = 64*key_num+64*max_other_num+64;
    sbMatrix entry_key_sorted(len1+len2, entry_key_bitcount);
    arith2bool(pIdx, entry_key_sorted_si, entry_key_sorted, enc, eval, runtime);

    //------permutation sort end------
    std::cout<<"permutation sort end"<<std::endl;


    sbMatrix j_sb(len1+len2, 64*key_num);
    sbMatrix tid_sb(len1+len2, 64);
    sbMatrix d_sb(len1+len2, 64*max_other_num);

    for(size_t i=0; i< len1+len2; i++){
        //j
        for(int j=entry_key_cols-1; j>=entry_key_cols-key_num; j--){
            j_sb.mShares[0](i, entry_key_cols-j-1) = entry_key_sorted.mShares[0](i, j);
            j_sb.mShares[1](i, entry_key_cols-j-1) = entry_key_sorted.mShares[1](i, j);
        }
        //tid
        tid_sb.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, entry_key_cols-key_num-1);
        tid_sb.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, entry_key_cols-key_num-1);
        //d
        for(int j=entry_key_cols-key_num-2; j>=0; j--){
            d_sb.mShares[0](i, entry_key_cols-key_num-2-j) = entry_key_sorted.mShares[0](i, j);
            d_sb.mShares[1](i, entry_key_cols-key_num-2-j) = entry_key_sorted.mShares[1](i, j);
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double step2_decompose_sort_result_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    //timing_results.push_back({"augment_step2_sort_result", step2_decompose_sort_result_time});
    timing_results.push_back({"augment_step2_perm_sort_result", step2_decompose_sort_result_time});

    //step 3: full-dimension
    t1 = std::chrono::high_resolution_clock::now();
    sbMatrix alpha_1(len1+len2, 64), alpha_2(len1+len2, 64);

    //same_attr:
    sbMatrix false_matrix(1,1);
    bool_init_false(pIdx, false_matrix);
    sbMatrix same_attr_partial(len1+len2-1, 1);
    compare_consecutive_rows_bool(pIdx, j_sb, same_attr_partial, enc, eval, runtime);
    sbMatrix same_attr(len1+len2, 1), not_same_attr(len1+len2, 1);
    same_attr.mShares[0](0, 0) = false_matrix.mShares[0](0, 0);
    same_attr.mShares[1](0, 0) = false_matrix.mShares[1](0, 0);
    std::memcpy(same_attr.mShares[0].data() + 1, same_attr_partial.mShares[0].data(), (len1+len2-1) * sizeof(same_attr_partial.mShares[0](0, 0)));
    std::memcpy(same_attr.mShares[1].data() + 1, same_attr_partial.mShares[1].data(), (len1+len2-1) * sizeof(same_attr_partial.mShares[1](0, 0)));
    
    
    sbMatrix one_bit(len1+len2, 1);
    bool_init_true(pIdx, one_bit);
    bool_cipher_sub(pIdx, one_bit, same_attr, not_same_attr, enc, eval, runtime);

    //is_table
    // sbMatrix zero64_vector(len1+len2, 64), one64_vector(len1+len2, 64);
    // bool_init_false(pIdx, zero64_vector);
    // bool_init_true(pIdx, one64_vector);
    // sbMatrix is_table_1(len1+len2, 1), is_table_2(len1+len2, 1);
    // //TODO VECTOR
    // bool_cipher_eq(pIdx, tid_sb, zero64_vector, is_table_1, enc, eval, runtime);
    // bool_cipher_eq(pIdx, tid_sb, one64_vector, is_table_2, enc, eval, runtime);
    sbMatrix zero64_vector(len1+len2, 64);
    bool_init_false(pIdx, zero64_vector);
    sbMatrix is_table_1(len1+len2, 1), is_table_2(len1+len2, 1);
    bool_cipher_eq(pIdx, tid_sb, zero64_vector, is_table_1, enc, eval, runtime);
    //要么是1要么是2
    bool_cipher_sub(pIdx, one_bit, is_table_1, is_table_2, enc, eval, runtime);


    //condition:
    sbMatrix is_table_vector(4*(len1+len2), 1), same_attr_vector(4*(len1+len2), 1);
    std::memcpy(is_table_vector.mShares[0].data(), is_table_1.mShares[0].data(), (len1+len2) * sizeof(is_table_1.mShares[0](0, 0)));
    std::memcpy(is_table_vector.mShares[1].data(), is_table_1.mShares[1].data(), (len1+len2) * sizeof(is_table_1.mShares[1](0, 0)));
    std::memcpy(same_attr_vector.mShares[0].data(), not_same_attr.mShares[0].data(), (len1+len2) * sizeof(not_same_attr.mShares[0](0, 0)));
    std::memcpy(same_attr_vector.mShares[1].data(), not_same_attr.mShares[1].data(), (len1+len2) * sizeof(not_same_attr.mShares[1](0, 0)));
    std::memcpy(is_table_vector.mShares[0].data() + (len1+len2), is_table_1.mShares[0].data(), (len1+len2) * sizeof(is_table_1.mShares[0](0, 0)));
    std::memcpy(is_table_vector.mShares[1].data() + (len1+len2), is_table_1.mShares[1].data(), (len1+len2) * sizeof(is_table_1.mShares[1](0, 0)));
    std::memcpy(same_attr_vector.mShares[0].data() + (len1+len2), same_attr.mShares[0].data(), (len1+len2) * sizeof(same_attr.mShares[0](0, 0)));
    std::memcpy(same_attr_vector.mShares[1].data() + (len1+len2), same_attr.mShares[1].data(), (len1+len2) * sizeof(same_attr.mShares[1](0, 0)));
    std::memcpy(is_table_vector.mShares[0].data() + 2*(len1+len2), is_table_2.mShares[0].data(), (len1+len2) * sizeof(is_table_2.mShares[0](0, 0)));
    std::memcpy(is_table_vector.mShares[1].data() + 2*(len1+len2), is_table_2.mShares[1].data(), (len1+len2) * sizeof(is_table_2.mShares[1](0, 0)));
    std::memcpy(same_attr_vector.mShares[0].data() + 2*(len1+len2), not_same_attr.mShares[0].data(), (len1+len2) * sizeof(not_same_attr.mShares[0](0, 0)));
    std::memcpy(same_attr_vector.mShares[1].data() + 2*(len1+len2), not_same_attr.mShares[1].data(), (len1+len2) * sizeof(not_same_attr.mShares[1](0, 0)));
    std::memcpy(is_table_vector.mShares[0].data() + 3*(len1+len2), is_table_2.mShares[0].data(), (len1+len2) * sizeof(is_table_2.mShares[0](0, 0)));
    std::memcpy(is_table_vector.mShares[1].data() + 3*(len1+len2), is_table_2.mShares[1].data(), (len1+len2) * sizeof(is_table_2.mShares[1](0, 0)));
    std::memcpy(same_attr_vector.mShares[0].data() + 3*(len1+len2), same_attr.mShares[0].data(), (len1+len2) * sizeof(same_attr.mShares[0](0, 0)));
    std::memcpy(same_attr_vector.mShares[1].data() + 3*(len1+len2), same_attr.mShares[1].data(), (len1+len2) * sizeof(same_attr.mShares[1](0, 0)));
    sbMatrix condition(4*(len1+len2), 1);
    bool_cipher_and(pIdx, is_table_vector, same_attr_vector, condition, enc, eval, runtime);
    sbMatrix cond1(len1+len2, 1), cond2(len1+len2, 1), cond3(len1+len2, 1), cond4(len1+len2, 1);
    std::memcpy(cond1.mShares[0].data(), condition.mShares[0].data(), (len1+len2) * sizeof(condition.mShares[0](0, 0)));
    std::memcpy(cond1.mShares[1].data(), condition.mShares[1].data(), (len1+len2) * sizeof(condition.mShares[1](0, 0)));
    std::memcpy(cond2.mShares[0].data(), condition.mShares[0].data() + (len1+len2), (len1+len2) * sizeof(condition.mShares[0](0, 0)));
    std::memcpy(cond2.mShares[1].data(), condition.mShares[1].data() + (len1+len2), (len1+len2) * sizeof(condition.mShares[1](0, 0)));
    std::memcpy(cond3.mShares[0].data(), condition.mShares[0].data() + 2*(len1+len2), (len1+len2) * sizeof(condition.mShares[0](0, 0)));
    std::memcpy(cond3.mShares[1].data(), condition.mShares[1].data() + 2*(len1+len2), (len1+len2) * sizeof(condition.mShares[1](0, 0)));
    std::memcpy(cond4.mShares[0].data(), condition.mShares[0].data() + 3*(len1+len2), (len1+len2) * sizeof(condition.mShares[0](0, 0)));
    std::memcpy(cond4.mShares[1].data(), condition.mShares[1].data() + 3*(len1+len2), (len1+len2) * sizeof(condition.mShares[1](0, 0)));

    //downward scan: 
    for(size_t i=0; i<len1+len2; i++){
        // sbMatrix j_i(1, 64*key_num), tid_i(1, 64);
        // for(int j=0; j<key_num; j++){
        //     j_i.mShares[0](0, j) = j_sb.mShares[0](i, j);
        //     j_i.mShares[1](0, j) = j_sb.mShares[1](i, j);
        // }
        // tid_i.mShares[0](0, 0) = tid_sb.mShares[0](i, 0);
        // tid_i.mShares[1](0, 0) = tid_sb.mShares[1](i, 0);

        // //join_key_same
        // sbMatrix same_attr(1, 1), not_same_attr(1,1);
        // if(i == 0){
        //     bool_init_false(pIdx, same_attr);
        // } else {
        //     sbMatrix j_i_prev(1, 64*key_num);
        //     for(int j=0; j<key_num; j++){
        //         j_i_prev.mShares[0](0, j) = j_sb.mShares[0](i-1, j);
        //         j_i_prev.mShares[1](0, j) = j_sb.mShares[1](i-1, j);
        //     }
        //     bool_cipher_eq(pIdx, j_i, j_i_prev, same_attr, enc, eval, runtime);
        // }
        // sbMatrix one_bit(1, 1);
        // bool_init_true(pIdx, one_bit);
        // bool_cipher_sub(pIdx, one_bit, same_attr, not_same_attr, enc, eval, runtime);

        // //tid
       
        // //TODO VECTOR
        // sbMatrix is_table_1(1, 1);
        // bool_cipher_eq(pIdx, tid_i, zero64, is_table_1, enc, eval, runtime);
        // sbMatrix is_table_2(1, 1);
        // bool_cipher_eq(pIdx, tid_i, one64, is_table_2, enc, eval, runtime);

        // //condition
        // sbMatrix cond1, cond2, cond3, cond4;
        // //TODO VECTOR
        // //TODO condition在for之前算
        // bool_cipher_and(pIdx, is_table_1, not_same_attr, cond1, enc, eval, runtime);
        // bool_cipher_and(pIdx, is_table_1, same_attr, cond2, enc, eval, runtime);
        // bool_cipher_and(pIdx, is_table_2, not_same_attr, cond3, enc, eval, runtime);
        // bool_cipher_and(pIdx, is_table_2, same_attr, cond4, enc, eval, runtime);


        sbMatrix zero64(1, 64);
        sbMatrix one64(1, 64);
        bool_init_false(pIdx, zero64);
        bool_init_true(pIdx, one64);
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
        
        sbMatrix cond_i_vector(4, 64);
        cond_i_vector.mShares[0](0, 0) = (cond1.mShares[0](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[1](0, 0) = (cond1.mShares[1](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[0](1, 0) = (cond2.mShares[0](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[1](1, 0) = (cond2.mShares[1](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[0](2, 0) = (cond3.mShares[0](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[1](2, 0) = (cond3.mShares[1](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[0](3, 0) = (cond4.mShares[0](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[1](3, 0) = (cond4.mShares[1](i, 0) == 1) ? -1 : -0;
        sbMatrix alpha_1_value(4, 64);
        alpha_1_value.mShares[0](0, 0) = one64.mShares[0](0, 0);
        alpha_1_value.mShares[1](0, 0) = one64.mShares[1](0, 0);
        alpha_1_value.mShares[0](1, 0) = alpha_1_plus.mShares[0](0, 0);
        alpha_1_value.mShares[1](1, 0) = alpha_1_plus.mShares[1](0, 0);
        alpha_1_value.mShares[0](2, 0) = zero64.mShares[0](0, 0);
        alpha_1_value.mShares[1](2, 0) = zero64.mShares[1](0, 0);
        alpha_1_value.mShares[0](3, 0) = current_alpha_1.mShares[0](0, 0);
        alpha_1_value.mShares[1](3, 0) = current_alpha_1.mShares[1](0, 0);
        sbMatrix term_1(4, 64);
        bool_cipher_and(pIdx, cond_i_vector, alpha_1_value, term_1, enc, eval, runtime);
        // bool_cipher_and(pIdx, cond1, one64, term1, enc, eval, runtime);
        // bool_cipher_and(pIdx, cond2, alpha_1_plus, term2, enc, eval, runtime);
        // bool_cipher_and(pIdx, zero64, cond3, term3, enc, eval, runtime);
        // bool_cipher_and(pIdx, current_alpha_1, cond4, term4, enc, eval, runtime);
       
        alpha_1.mShares[0](i, 0) = term_1.mShares[0](0, 0) ^ term_1.mShares[0](1, 0)^ term_1.mShares[0](2, 0)^ term_1.mShares[0](3, 0);
        alpha_1.mShares[1](i, 0) = term_1.mShares[1](0, 0) ^ term_1.mShares[1](1, 0)^ term_1.mShares[1](2, 0)^ term_1.mShares[1](3, 0);
        

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

        sbMatrix alpha_2_value(4, 64);
        alpha_2_value.mShares[0](0, 0) = zero64.mShares[0](0, 0);
        alpha_2_value.mShares[1](0, 0) = zero64.mShares[1](0, 0);
        alpha_2_value.mShares[0](1, 0) = zero64.mShares[0](0, 0);
        alpha_2_value.mShares[1](1, 0) = zero64.mShares[1](0, 0);
        alpha_2_value.mShares[0](2, 0) = one64.mShares[0](0, 0);
        alpha_2_value.mShares[1](2, 0) = one64.mShares[1](0, 0);
        alpha_2_value.mShares[0](3, 0) = alpha_2_plus.mShares[0](0, 0);
        alpha_2_value.mShares[1](3, 0) = alpha_2_plus.mShares[1](0, 0);
        sbMatrix term_2(4, 64);
        bool_cipher_and(pIdx, cond_i_vector, alpha_2_value, term_2, enc, eval, runtime);

        // bool_cipher_and(pIdx, zero64, cond1, term5, enc, eval, runtime);
        // bool_cipher_and(pIdx, zero64, cond2, term6, enc, eval, runtime);
        // bool_cipher_and(pIdx, one64, cond3, term7, enc, eval, runtime);
        // bool_cipher_and(pIdx, alpha_2_plus, cond4, term8, enc, eval, runtime);
        
        alpha_2.mShares[0](i, 0) = term_2.mShares[0](0, 0) ^ term_2.mShares[0](1, 0)^ term_2.mShares[0](2, 0)^ term_2.mShares[0](3, 0);
        alpha_2.mShares[1](i, 0) = term_2.mShares[1](0, 0) ^ term_2.mShares[1](1, 0)^ term_2.mShares[1](2, 0)^ term_2.mShares[1](3, 0);
        
    }
    t2 = std::chrono::high_resolution_clock::now();
    double step3_downward_scan_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"augment_step3_downward_scan", step3_downward_scan_time});

    //upward scan:
    t1 = std::chrono::high_resolution_clock::now(); 

    sbMatrix same_inv_attr(len1+len2, 1), not_same_inv_attr(len1+len2, 1);
    std::memcpy(same_inv_attr.mShares[0].data() , same_attr_partial.mShares[0].data(), (len1+len2-1) * sizeof(same_attr_partial.mShares[0](0, 0)));
    std::memcpy(same_inv_attr.mShares[1].data() , same_attr_partial.mShares[1].data(), (len1+len2-1) * sizeof(same_attr_partial.mShares[1](0, 0)));
    same_inv_attr.mShares[0](len1+len2-1, 0) = false_matrix.mShares[0](0, 0);
    same_inv_attr.mShares[1](len1+len2-1, 0) = false_matrix.mShares[1](0, 0);
    bool_cipher_sub(pIdx, one_bit, same_inv_attr, not_same_inv_attr, enc, eval, runtime);

     //condition:
    sbMatrix same_inv_attr_vector(4*(len1+len2), 1);
    std::memcpy(same_inv_attr_vector.mShares[0].data(), not_same_inv_attr.mShares[0].data(), (len1+len2) * sizeof(not_same_inv_attr.mShares[0](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[1].data(), not_same_inv_attr.mShares[1].data(), (len1+len2) * sizeof(not_same_inv_attr.mShares[1](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[0].data() + (len1+len2), same_inv_attr.mShares[0].data(), (len1+len2) * sizeof(same_inv_attr.mShares[0](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[1].data() + (len1+len2), same_inv_attr.mShares[1].data(), (len1+len2) * sizeof(same_inv_attr.mShares[1](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[0].data() + 2*(len1+len2), not_same_inv_attr.mShares[0].data(), (len1+len2) * sizeof(not_same_inv_attr.mShares[0](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[1].data() + 2*(len1+len2), not_same_inv_attr.mShares[1].data(), (len1+len2) * sizeof(not_same_inv_attr.mShares[1](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[0].data() + 3*(len1+len2), same_inv_attr.mShares[0].data(), (len1+len2) * sizeof(same_inv_attr.mShares[0](0, 0)));
    std::memcpy(same_inv_attr_vector.mShares[1].data() + 3*(len1+len2), same_inv_attr.mShares[1].data(), (len1+len2) * sizeof(same_inv_attr.mShares[1](0, 0)));
    

    bool_cipher_and(pIdx, is_table_vector, same_inv_attr_vector, condition, enc, eval, runtime);
    std::memcpy(cond1.mShares[0].data(), condition.mShares[0].data(), (len1+len2) * sizeof(condition.mShares[0](0, 0)));
    std::memcpy(cond1.mShares[1].data(), condition.mShares[1].data(), (len1+len2) * sizeof(condition.mShares[1](0, 0)));
    std::memcpy(cond2.mShares[0].data(), condition.mShares[0].data() + (len1+len2), (len1+len2) * sizeof(condition.mShares[0](0, 0)));
    std::memcpy(cond2.mShares[1].data(), condition.mShares[1].data() + (len1+len2), (len1+len2) * sizeof(condition.mShares[1](0, 0)));
    std::memcpy(cond3.mShares[0].data(), condition.mShares[0].data() + 2*(len1+len2), (len1+len2) * sizeof(condition.mShares[0](0, 0)));
    std::memcpy(cond3.mShares[1].data(), condition.mShares[1].data() + 2*(len1+len2), (len1+len2) * sizeof(condition.mShares[1](0, 0)));
    std::memcpy(cond4.mShares[0].data(), condition.mShares[0].data() + 3*(len1+len2), (len1+len2) * sizeof(condition.mShares[0](0, 0)));
    std::memcpy(cond4.mShares[1].data(), condition.mShares[1].data() + 3*(len1+len2), (len1+len2) * sizeof(condition.mShares[1](0, 0)));
 

    //upward scan:
    for(int i=len1+len2-1; i>=0; i--){
        // sbMatrix j_i(1, 64*key_num), tid_i(1, 64);
        // for(int j=0; j<key_num; j++){
        //     j_i.mShares[0](0, j) = j_sb.mShares[0](i, j);
        //     j_i.mShares[1](0, j) = j_sb.mShares[1](i, j);
        // }
        // tid_i.mShares[0](0, 0) = tid_sb.mShares[0](i, 0);
        // tid_i.mShares[1](0, 0) = tid_sb.mShares[1](i, 0);

        // //join_key_same
        // sbMatrix same_attr(1, 1), not_same_attr(1,1);
        // if(i == len1+len2-1){
        //     bool_init_false(pIdx, same_attr);
        // } else {
        //     sbMatrix j_i_prev(1, 64*key_num);
        //     for(int j=0; j<key_num; j++){
        //         j_i_prev.mShares[0](0, j) = j_sb.mShares[0](i+1, j);
        //         j_i_prev.mShares[1](0, j) = j_sb.mShares[1](i+1, j);
        //     }
        //     bool_cipher_eq(pIdx, j_i, j_i_prev, same_attr, enc, eval, runtime);
        // }
        // sbMatrix one_bit(1, 1);
        // bool_init_true(pIdx, one_bit);
        // bool_cipher_sub(pIdx, one_bit, same_attr, not_same_attr, enc, eval, runtime);
        
        // //tid
        // sbMatrix zero64(1, 64);
        // sbMatrix one64(1, 64);
        // bool_init_false(pIdx, zero64);
        // bool_init_true(pIdx, one64);
        // sbMatrix is_table_1(1, 1);
        // bool_cipher_eq(pIdx, tid_i, zero64, is_table_1, enc, eval, runtime);
        // sbMatrix is_table_2(1, 1);
        // bool_cipher_eq(pIdx, tid_i, one64, is_table_2, enc, eval, runtime);

        //condition
        // sbMatrix cond1, cond2, cond3, cond4;
        // bool_cipher_and(pIdx, is_table_1, not_same_attr, cond1, enc, eval, runtime);
        // bool_cipher_and(pIdx, is_table_1, same_attr, cond2, enc, eval, runtime);
        // bool_cipher_and(pIdx, is_table_2, not_same_attr, cond3, enc, eval, runtime);
        // bool_cipher_and(pIdx, is_table_2, same_attr, cond4, enc, eval, runtime);   
        // cond1.resize(1, 64);
        // cond2.resize(1, 64);
        // cond3.resize(1, 64);
        // cond4.resize(1, 64);
        // cond1.mShares[0](0, 0) = (cond1.mShares[0](0, 0) == 1) ? -1 : -0;
        // cond1.mShares[1](0, 0) = (cond1.mShares[1](0, 0) == 1) ? -1 : -0;
        // cond2.mShares[0](0, 0) = (cond2.mShares[0](0, 0) == 1) ? -1 : -0;
        // cond2.mShares[1](0, 0) = (cond2.mShares[1](0, 0) == 1) ? -1 : -0;
        // cond3.mShares[0](0, 0) = (cond3.mShares[0](0, 0) == 1) ? -1 : -0;
        // cond3.mShares[1](0, 0) = (cond3.mShares[1](0, 0) == 1) ? -1 : -0;
        // cond4.mShares[0](0, 0) = (cond4.mShares[0](0, 0) == 1) ? -1 : -0;
        // cond4.mShares[1](0, 0) = (cond4.mShares[1](0, 0) == 1) ? -1 : -0;

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

        sbMatrix cond_i_vector(4, 64);
        cond_i_vector.mShares[0](0, 0) = (cond1.mShares[0](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[1](0, 0) = (cond1.mShares[1](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[0](1, 0) = (cond2.mShares[0](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[1](1, 0) = (cond2.mShares[1](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[0](2, 0) = (cond3.mShares[0](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[1](2, 0) = (cond3.mShares[1](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[0](3, 0) = (cond4.mShares[0](i, 0) == 1) ? -1 : -0;
        cond_i_vector.mShares[1](3, 0) = (cond4.mShares[1](i, 0) == 1) ? -1 : -0;
        sbMatrix alpha_1_value(4, 64);
        alpha_1_value.mShares[0](0, 0) = current_alpha_1.mShares[0](0, 0);
        alpha_1_value.mShares[1](0, 0) = current_alpha_1.mShares[1](0, 0);
        alpha_1_value.mShares[0](1, 0) = pre_alpha_1.mShares[0](0, 0);
        alpha_1_value.mShares[1](1, 0) = pre_alpha_1.mShares[1](0, 0);
        alpha_1_value.mShares[0](2, 0) = current_alpha_1.mShares[0](0, 0);
        alpha_1_value.mShares[1](2, 0) = current_alpha_1.mShares[1](0, 0);
        alpha_1_value.mShares[0](3, 0) = current_alpha_1.mShares[0](0, 0);
        alpha_1_value.mShares[1](3, 0) = current_alpha_1.mShares[1](0, 0);

        sbMatrix term_1(4, 64);
        bool_cipher_and(pIdx, cond_i_vector, alpha_1_value, term_1, enc, eval, runtime);
        // bool_cipher_and(pIdx, current_alpha_1, cond1, term1, enc, eval, runtime);
        // bool_cipher_and(pIdx, pre_alpha_1, cond2, term2, enc, eval, runtime);
        // bool_cipher_and(pIdx, current_alpha_1, cond3, term3, enc, eval, runtime);
        // bool_cipher_and(pIdx, current_alpha_1, cond4, term4, enc, eval, runtime);

        alpha_1.mShares[0](i, 0) = term_1.mShares[0](0, 0) ^ term_1.mShares[0](1, 0)^ term_1.mShares[0](2, 0)^ term_1.mShares[0](3, 0);
        alpha_1.mShares[1](i, 0) = term_1.mShares[1](0, 0) ^ term_1.mShares[1](1, 0)^ term_1.mShares[1](2, 0)^ term_1.mShares[1](3, 0);
        
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

        sbMatrix alpha_2_value(4, 64);
        alpha_2_value.mShares[0](0, 0) = current_alpha_2.mShares[0](0, 0);
        alpha_2_value.mShares[1](0, 0) = current_alpha_2.mShares[1](0, 0);
        alpha_2_value.mShares[0](1, 0) = pre_alpha_2.mShares[0](0, 0);
        alpha_2_value.mShares[1](1, 0) = pre_alpha_2.mShares[1](0, 0);
        alpha_2_value.mShares[0](2, 0) = current_alpha_2.mShares[0](0, 0);
        alpha_2_value.mShares[1](2, 0) = current_alpha_2.mShares[1](0, 0);
        alpha_2_value.mShares[0](3, 0) = pre_alpha_2.mShares[0](0, 0);
        alpha_2_value.mShares[1](3, 0) = pre_alpha_2.mShares[1](0, 0);
        sbMatrix term_2(4, 64);
        bool_cipher_and(pIdx, cond_i_vector, alpha_2_value, term_2, enc, eval, runtime);

        // sbMatrix term5(1, 64), term6(1, 64), term7(1, 64), term8(1, 64);
        // bool_cipher_and(pIdx, current_alpha_2, cond1, term5, enc, eval, runtime);
        // bool_cipher_and(pIdx, pre_alpha_2, cond2, term6, enc, eval, runtime);
        // bool_cipher_and(pIdx, current_alpha_2, cond3, term7, enc, eval, runtime);
        // bool_cipher_and(pIdx, pre_alpha_2, cond4, term8, enc, eval, runtime);

        alpha_2.mShares[0](i, 0) = term_2.mShares[0](0, 0) ^ term_2.mShares[0](1, 0)^ term_2.mShares[0](2, 0)^ term_2.mShares[0](3, 0);
        alpha_2.mShares[1](i, 0) = term_2.mShares[1](0, 0) ^ term_2.mShares[1](1, 0)^ term_2.mShares[1](2, 0)^ term_2.mShares[1](3, 0);
        
    }
    t2 = std::chrono::high_resolution_clock::now();
    double step3_upward_scan_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"augment_step3_upward_scan", step3_upward_scan_time});

    //step 4: sort by: tid
    t1 = std::chrono::high_resolution_clock::now();
    //还要加上alpha1,2
    sbMatrix entry_key(len1+len2, entry_key_bitcount+64+64);
    entry_key_sorted.resize(len1+len2,entry_key_bitcount+64+64);
    for(int i=0; i<len1+len2; i++){
        //因为是大端序，所以要逆序
        //tid
        entry_key.mShares[0](i, entry_key_cols+1) = tid_sb.mShares[0](i, 0);
        entry_key.mShares[1](i, entry_key_cols+1) = tid_sb.mShares[1](i, 0);
        //j
        for(int j=entry_key_cols; j>=entry_key_cols-key_num+1; j--){
            entry_key.mShares[0](i, j) = j_sb.mShares[0](i, entry_key_cols-j);
            entry_key.mShares[1](i, j) = j_sb.mShares[1](i, entry_key_cols-j);
        }
        //d
        for(int j=entry_key_cols-key_num; j>=2; j--){
            entry_key.mShares[0](i, j) = d_sb.mShares[0](i, entry_key_cols-key_num-j);
            entry_key.mShares[1](i, j) = d_sb.mShares[1](i, entry_key_cols-key_num-j);
        }
        //alpha1
        entry_key.mShares[0](i, 1) = alpha_1.mShares[0](i, 0);
        entry_key.mShares[1](i, 1) = alpha_1.mShares[1](i, 0);
        //alpha2
        entry_key.mShares[0](i, 0) = alpha_2.mShares[0](i, 0);
        entry_key.mShares[1](i, 0) = alpha_2.mShares[1](i, 0);
    }
    odd_even_merge_sort(entry_key, entry_key_sorted, pIdx, enc, eval, runtime);

    sbMatrix j_sorted(len1+len2, 64*key_num), d_sorted(len1+len2, 64*max_other_num),tid_sorted(len1+len2, 64);
    sbMatrix alpha_1_sorted(len1+len2, 64), alpha_2_sorted(len1+len2, 64);
 
    for(size_t i=0; i< len1+len2; i++){
        //tid
        tid_sorted.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, entry_key_cols+1);
        tid_sorted.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, entry_key_cols+1);
        //j
        for(int j=entry_key_cols; j>=entry_key_cols-key_num+1; j--){
            j_sorted.mShares[0](i, entry_key_cols-j) = entry_key_sorted.mShares[0](i, j);
            j_sorted.mShares[1](i, entry_key_cols-j) = entry_key_sorted.mShares[1](i, j);
        }
        //d
        for(int j=entry_key_cols-key_num; j>=2; j--){
            d_sorted.mShares[0](i, entry_key_cols-key_num-j) = entry_key_sorted.mShares[0](i, j);
            d_sorted.mShares[1](i, entry_key_cols-key_num-j) = entry_key_sorted.mShares[1](i, j);
        }
        //alpha1
        alpha_1_sorted.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, 1);
        alpha_1_sorted.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, 1);
        //alpha2
        alpha_2_sorted.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, 0);
        alpha_2_sorted.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, 0);
    }
    t2 = std::chrono::high_resolution_clock::now();
    double step4_decompose_sort_result_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"augment_step4_sort_result", step4_decompose_sort_result_time});
   
    //step 5: abstract from T_c
    t1 = std::chrono::high_resolution_clock::now();
    T_1_auged.resize(4);
    T_2_auged.resize(4);
    T_1_auged[0].resize(len1, 64*key_num);
    T_1_auged[1].resize(len1, 64*max_other_num);
    T_2_auged[0].resize(len2, 64*key_num);
    T_2_auged[1].resize(len2, 64*max_other_num);
    for(size_t i=2; i<4; i++){
        T_1_auged[i].resize(len1, 64);
        T_2_auged[i].resize(len2, 64);
    }

    std::memcpy(T_1_auged[0].mShares[0].data(), j_sorted.mShares[0].data(), len1 * key_num * sizeof(j_sorted.mShares[0](0, 0)));
    std::memcpy(T_1_auged[0].mShares[1].data(), j_sorted.mShares[1].data(), len1 * key_num * sizeof(j_sorted.mShares[1](0, 0)));
    std::memcpy(T_1_auged[1].mShares[0].data(), d_sorted.mShares[0].data(), len1 * max_other_num * sizeof(d_sorted.mShares[0](0, 0)));
    std::memcpy(T_1_auged[1].mShares[1].data(), d_sorted.mShares[1].data(), len1 * max_other_num * sizeof(d_sorted.mShares[1](0, 0)));
    std::memcpy(T_1_auged[2].mShares[0].data(), alpha_1_sorted.mShares[0].data(), len1 * sizeof(alpha_1_sorted.mShares[0](0, 0)));
    std::memcpy(T_1_auged[2].mShares[1].data(), alpha_1_sorted.mShares[1].data(), len1 * sizeof(alpha_1_sorted.mShares[1](0, 0)));
    std::memcpy(T_1_auged[3].mShares[0].data(), alpha_2_sorted.mShares[0].data(), len1 * sizeof(alpha_2_sorted.mShares[0](0, 0)));
    std::memcpy(T_1_auged[3].mShares[1].data(), alpha_2_sorted.mShares[1].data(), len1 * sizeof(alpha_2_sorted.mShares[1](0, 0)));

    // for(size_t i=0; i<len1; i++){
    //     for(size_t j=0; j<key_num; j++){
    //         T_1_auged[0].mShares[0](i, j) = j_sorted.mShares[0](i, j);
    //         T_1_auged[0].mShares[1](i, j) = j_sorted.mShares[1](i, j);
    //     }
    //     for(size_t j=0; j<max_other_num; j++){
    //         T_1_auged[1].mShares[0](i, j) = d_sorted.mShares[0](i, j);
    //         T_1_auged[1].mShares[1](i, j) = d_sorted.mShares[1](i, j);
    //     }
    //     T_1_auged[2].mShares[0](i, 0) = alpha_1_sorted.mShares[0](i, 0);
    //     T_1_auged[2].mShares[1](i, 0) = alpha_1_sorted.mShares[1](i, 0);
    //     T_1_auged[3].mShares[0](i, 0) = alpha_2_sorted.mShares[0](i, 0);
    //     T_1_auged[3].mShares[1](i, 0) = alpha_2_sorted.mShares[1](i, 0);
    // }

    std::memcpy(T_2_auged[0].mShares[0].data(), j_sorted.mShares[0].data() + len1 * key_num, len2 * key_num * sizeof(j_sorted.mShares[0](0, 0)));
    std::memcpy(T_2_auged[0].mShares[1].data(), j_sorted.mShares[1].data() + len1 * key_num, len2 * key_num * sizeof(j_sorted.mShares[1](0, 0)));
    std::memcpy(T_2_auged[1].mShares[0].data(), d_sorted.mShares[0].data() + len1 * max_other_num, len2 * max_other_num * sizeof(d_sorted.mShares[0](0, 0)));
    std::memcpy(T_2_auged[1].mShares[1].data(), d_sorted.mShares[1].data() + len1 * max_other_num, len2 * max_other_num * sizeof(d_sorted.mShares[1](0, 0)));
    std::memcpy(T_2_auged[2].mShares[0].data(), alpha_1_sorted.mShares[0].data() + len1, len2 * sizeof(alpha_1_sorted.mShares[0](0, 0)));
    std::memcpy(T_2_auged[2].mShares[1].data(), alpha_1_sorted.mShares[1].data() + len1, len2 * sizeof(alpha_1_sorted.mShares[1](0, 0)));
    std::memcpy(T_2_auged[3].mShares[0].data(), alpha_2_sorted.mShares[0].data() + len1, len2 * sizeof(alpha_2_sorted.mShares[0](0, 0)));
    std::memcpy(T_2_auged[3].mShares[1].data(), alpha_2_sorted.mShares[1].data() + len1, len2 * sizeof(alpha_2_sorted.mShares[1](0, 0)));

    // for(size_t i=0; i<len2; i++){
    //     for(size_t j=0; j<key_num; j++){
    //         T_2_auged[0].mShares[0](i, j) = j_sorted.mShares[0](i+len1, j);
    //         T_2_auged[0].mShares[1](i, j) = j_sorted.mShares[1](i+len1, j);
    //     }
    //     for(size_t j=0; j<max_other_num; j++){
    //         T_2_auged[1].mShares[0](i, j) = d_sorted.mShares[0](i+len1, j);
    //         T_2_auged[1].mShares[1](i, j) = d_sorted.mShares[1](i+len1, j);
    //     }
    //     T_2_auged[2].mShares[0](i, 0) = alpha_1_sorted.mShares[0](i+len1, 0);
    //     T_2_auged[2].mShares[1](i, 0) = alpha_1_sorted.mShares[1](i+len1, 0);
    //     T_2_auged[3].mShares[0](i, 0) = alpha_2_sorted.mShares[0](i+len1, 0);
    //     T_2_auged[3].mShares[1](i, 0) = alpha_2_sorted.mShares[1](i+len1, 0);
    // }
    t2 = std::chrono::high_resolution_clock::now();
    double step5_abstract_from_Tc_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"augment_step5_abstract_from_Tc", step5_abstract_from_Tc_time});

    // 计算总时间
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    timing_results.push_back({"augment_total", total_time});

    // 将时间测量结果写入文件（仅 role 0 写入）
    if (pIdx == 0) {
        std::string filename = "./join_timing_results_role0.txt";
        std::ofstream outFile(filename, std::ios::app);
        if (outFile.is_open()) {
            outFile << "--- augment_table Internal Timing (Role " << pIdx << ") ---" << std::endl;
            outFile << std::fixed << std::setprecision(3);
            for(const auto& result : timing_results) {
                outFile << "  " << result.first << ": " << result.second << " ms" << std::endl;
            }
            outFile << "----------------------------------------" << std::endl;
            outFile << std::endl;
            outFile.close();
        }
    }
    
    // //DEBUG
    // //T_1_auged  
    // for(int l=0; l<key_num ; l++){
    //     i64Matrix j_sb_plain(len1, 1);
    //     sbMatrix j_sb_l(len1, 64);
    //     for(int i=0; i<len1; i++){
    //         j_sb_l.mShares[0](i, 0) = T_1_auged[0].mShares[0](i, l);
    //         j_sb_l.mShares[1](i, 0) = T_1_auged[0].mShares[1](i, l);
    //     }
    //     enc.revealAll(runtime, j_sb_l, j_sb_plain).get();
    //     if(pIdx ==0){
    //         std::cout<< "j_sb_"<<l<<": "<<std::endl;
    //         for(size_t i=0; i<len1; i++){
    //             std::cout << j_sb_plain(i, 0) << " ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }
    // for(int l=0; l<max_other_num ; l++){
    //     i64Matrix d_sb_plain(len1, 1);
    //     sbMatrix d_sb_l(len1, 64);
    //     for(int i=0; i<len1; i++){
    //         d_sb_l.mShares[0](i, 0) = T_1_auged[1].mShares[0](i, l);
    //         d_sb_l.mShares[1](i, 0) = T_1_auged[1].mShares[1](i, l);
    //     }
    //     enc.revealAll(runtime, d_sb_l, d_sb_plain).get();
    //     if(pIdx ==0){
    //         std::cout<< "d_sb_"<<l<<": "<<std::endl;
    //         for(size_t i=0; i<len1; i++){
    //             std::cout << d_sb_plain(i, 0) << " ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }

    
    // i64Matrix t1_alpha_1_plain(len1, 1);
    // sbMatrix t1_alpha_1(len1, 64);
    // for(int i=0; i<len1; i++){
    //     t1_alpha_1.mShares[0](i, 0) = T_1_auged[2].mShares[0](i, 0);
    //     t1_alpha_1.mShares[1](i, 0) = T_1_auged[2].mShares[1](i, 0);
    // }
    // enc.revealAll(runtime, t1_alpha_1, t1_alpha_1_plain).get();
    // if(pIdx ==0){
    //     std::cout << "t1_alpha_1_plain: " << std::endl;
    //     for(size_t i=0; i<len1; i++){
    //         std::cout << t1_alpha_1_plain(i, 0) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // i64Matrix t1_alpha_2_plain(len1, 1);
    // sbMatrix t1_alpha_2(len1, 64);
    // for(int i=0; i<len1; i++){
    //     t1_alpha_2.mShares[0](i, 0) = T_1_auged[3].mShares[0](i, 0);
    //     t1_alpha_2.mShares[1](i, 0) = T_1_auged[3].mShares[1](i, 0);
    // }
    // enc.revealAll(runtime, t1_alpha_2, t1_alpha_2_plain).get();
    // if(pIdx ==0){
    //     std::cout << "t1_alpha_2_plain: " << std::endl;
    //     for(size_t i=0; i<len1; i++){
    //         std::cout << t1_alpha_2_plain(i, 0) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // //T_2_auged
    // for(int l=0; l<key_num ; l++){
    //     i64Matrix j_sb_plain(len2, 1);
    //     sbMatrix j_sb_l(len2, 64);
    //     for(int i=0; i<len2; i++){
    //         j_sb_l.mShares[0](i, 0) = T_2_auged[0].mShares[0](i, l);
    //         j_sb_l.mShares[1](i, 0) = T_2_auged[0].mShares[1](i, l);
    //     }
    //     enc.revealAll(runtime, j_sb_l, j_sb_plain).get();
    //     if(pIdx ==0){
    //         std::cout<< "j_sb_"<<l<<": "<<std::endl;
    //         for(size_t i=0; i<len2; i++){
    //             std::cout << j_sb_plain(i, 0) << " ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }
    // for(int l=0; l<max_other_num ; l++){
    //     i64Matrix d_sb_plain(len2, 1);
    //     sbMatrix d_sb_l(len2, 64);
    //     for(int i=0; i<len2; i++){
    //         d_sb_l.mShares[0](i, 0) = T_2_auged[1].mShares[0](i, l);
    //         d_sb_l.mShares[1](i, 0) = T_2_auged[1].mShares[1](i, l);
    //     }
    //     enc.revealAll(runtime, d_sb_l, d_sb_plain).get();
    //     if(pIdx ==0){
    //         std::cout<< "d_sb_"<<l<<": "<<std::endl;
    //         for(size_t i=0; i<len2; i++){
    //             std::cout << d_sb_plain(i, 0) << " ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }

    
    // i64Matrix t2_alpha_1_plain(len2, 1);
    // sbMatrix t2_alpha_1(len2, 64);
    // for(int i=0; i<len2; i++){
    //     t2_alpha_1.mShares[0](i, 0) = T_2_auged[2].mShares[0](i, 0);
    //     t2_alpha_1.mShares[1](i, 0) = T_2_auged[2].mShares[1](i, 0);
    // }
    // enc.revealAll(runtime, t2_alpha_1, t2_alpha_1_plain).get();
    // if(pIdx ==0){
    //     std::cout << "t2_alpha_1_plain: " << std::endl;
    //     for(size_t i=0; i<len2; i++){
    //         std::cout << t2_alpha_1_plain(i, 0) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // i64Matrix t2_alpha_2_plain(len2, 1);
    // sbMatrix t2_alpha_2(len2, 64);
    // for(int i=0; i<len2; i++){
    //     t2_alpha_2.mShares[0](i, 0) = T_2_auged[3].mShares[0](i, 0);
    //     t2_alpha_2.mShares[1](i, 0) = T_2_auged[3].mShares[1](i, 0);
    // }
    // enc.revealAll(runtime, t2_alpha_2, t2_alpha_2_plain).get();
    // if(pIdx ==0){
    //     std::cout << "t2_alpha_2_plain: " << std::endl;
    //     for(size_t i=0; i<len2; i++){
    //         std::cout << t2_alpha_2_plain(i, 0) << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // //----correct

    return ;
}
*/


//optimized:prefixsum(si) + selector
void oblivious_expand(int pIdx, std::vector<si64Matrix> &T, std::vector<sbMatrix> &A, i64 tid,
    i64Matrix &s_plain,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    // // 时间测量变量（按需保留，当前禁用）
    // auto total_start = std::chrono::high_resolution_clock::now();
    // std::vector<std::pair<std::string, double>> timing_results;
    
    int n = T[0].rows();
    int key_num = T[0].cols();
    int max_other_num = T[1].cols();

    //step 1: 计算f(x)和s
    // auto t1 = std::chrono::high_resolution_clock::now();
    //alpha的值64bit
    si64Matrix g(n, 1);
    if(tid == 0){
        g = T[3];
    } else {
        g = T[2];
    }

    si64Matrix zero(n, 1);
    sbMatrix flag(n,1), not_flag(n,1);
    set_const_share(pIdx, 0, zero, enc, eval, runtime);
    cipher_eq(pIdx, g, zero, flag, eval, runtime);
    bool_cipher_not(pIdx, flag, not_flag);
    // bool_cipher_not 对 i64 做 ~ 会把 0->-1、1->-2，asyncMul 要求 share 为 0/1，需取 LSB
    for(size_t i=0; i<n; i++){
        not_flag.mShares[0](i, 0) &= 1;
        not_flag.mShares[1](i, 0) &= 1;
    }

    si64Matrix g_sum(n, 1), s(1, 1), fx(n, 1);  
    i64 sum_0=0,sum_1=0;
    g_sum.mShares[0](0, 0) = sum_0;
    g_sum.mShares[1](0, 0) = sum_1;
    for(size_t i=1; i<n; i++){
        sum_0 += g.mShares[0](i-1, 0);
        sum_1 += g.mShares[1](i-1, 0);
        g_sum.mShares[0](i, 0) = sum_0;
        g_sum.mShares[1](i, 0) = sum_1;
    }

    //s=g_sum(n-1)
    s.mShares[0](0, 0) = g_sum.mShares[0](n-1, 0) + g.mShares[0](n-1, 0);
    s.mShares[1](0, 0) = g_sum.mShares[1](n-1, 0) + g.mShares[1](n-1, 0);
    s_plain.resize(1, 1);
    enc.revealAll(runtime, s, s_plain).get();
    int m = s_plain(0,0);
      

    //fx=(not flag) × g_sum；cipher_mul 要求 sbMatrix 的 bitCount==1
    cipher_mul(pIdx, g_sum, not_flag, fx, eval, enc, runtime);
    // auto t2 = std::chrono::high_resolution_clock::now();
    // double step1_compute_fx_s_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"oblivious_expand_step1_compute_fx_s", step1_compute_fx_s_time});


    //step 2 obli_distribute
    // t1 = std::chrono::high_resolution_clock::now();
    sbMatrix A_vector(m, 64*key_num+64*max_other_num+64+64);
    sbMatrix flag_sorted_auged(m, 1);
    oblivious_distribute(pIdx, T, flag, fx, s_plain, A_vector, flag_sorted_auged, enc, eval, runtime);
    int A_vector_i64cols = A_vector.i64Cols();
    // t2 = std::chrono::high_resolution_clock::now();
    // double step2_obli_distribute_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"oblivious_expand_step2_obli_distribute", step2_obli_distribute_time});


    //step 3 fill down
    // t1 = std::chrono::high_resolution_clock::now();
    // sbMatrix 第二参数是 bitCount，需 64*A_vector_i64cols 才能容纳 A_vector_i64cols 个 i64
    int A_vector_bitcount = 64 * A_vector_i64cols;
    sbMatrix px_vector(1, A_vector_bitcount);
    sbMatrix px(1, 64);
    bool_init_i64(pIdx, std::numeric_limits<i64>::max(), px, enc, eval, runtime);
    bool_cols_expand(px, px_vector);

    sbMatrix one_share(m, 1);
    bool_init_i64(pIdx, 1, one_share, enc, eval, runtime);
    sbMatrix cond(m, 1);
    bool_cipher_eq(pIdx, flag_sorted_auged, one_share, cond, enc, eval, runtime);


    // 3. 倍增循环：将 O(m) 降为 O(log m)
    for (int step = 1; step < m; step <<= 1) {
        sbMatrix prev_A(m, A_vector_bitcount);
        sbMatrix prev_cond(m, 1);

        // 策略：先整体拷贝 A_vector 到 prev_A，使 i < step 的部分自动“保持原样”
        std::memcpy(prev_A.mShares[0].data(), A_vector.mShares[0].data(), m * A_vector_i64cols * sizeof(i64));
        std::memcpy(prev_A.mShares[1].data(), A_vector.mShares[1].data(), m * A_vector_i64cols * sizeof(i64));
        
        // 同理拷贝 cond
        std::memcpy(prev_cond.mShares[0].data(), cond.mShares[0].data(), m * sizeof(i64));
        std::memcpy(prev_cond.mShares[1].data(), cond.mShares[1].data(), m * sizeof(i64));

        // 只更新 i >= step 的部分，将其指向“上方跳跃点”
        for (int i = step; i < m; ++i) {
            int target_idx = i - step;
            
            // 覆盖 prev_A 的第 i 行数据为 A_vector 的第 target_idx 行
            std::memcpy(prev_A.mShares[0].data() + i * A_vector_i64cols, 
                        A_vector.mShares[0].data() + target_idx * A_vector_i64cols, 
                        A_vector_i64cols * sizeof(i64));
            std::memcpy(prev_A.mShares[1].data() + i * A_vector_i64cols, 
                        A_vector.mShares[1].data() + target_idx * A_vector_i64cols, 
                        A_vector_i64cols * sizeof(i64));
                        
            // 覆盖 prev_cond 的第 i 行
            prev_cond.mShares[0](i, 0) = cond.mShares[0](target_idx, 0);
            prev_cond.mShares[1](i, 0) = cond.mShares[1](target_idx, 0);
        }


        sbMatrix next_A(m, A_vector_bitcount);
        bool_cipher_selector(pIdx, cond, prev_A, A_vector, next_A, enc, eval, runtime);   
        bool_cipher_and(pIdx, cond, prev_cond, cond, enc, eval, runtime);
        
        A_vector = next_A; 
    }

    //就这样得到A——vector不行吗，为啥要拆？不然不知道每一部分有多少啊
    //将A_vector拆分成A[k]
    A.resize(4);
    A[0].resize(m, 64*key_num);
    A[1].resize(m, 64*max_other_num);
    A[2].resize(m, 64);
    A[3].resize(m, 64);

    for(size_t i=0; i<m; i++){
        A[3].mShares[0](i, 0) = A_vector.mShares[0](i, 0);
        A[3].mShares[1](i, 0) = A_vector.mShares[1](i, 0);
        A[2].mShares[0](i, 0) = A_vector.mShares[0](i, 1);
        A[2].mShares[1](i, 0) = A_vector.mShares[1](i, 1);
        std::memcpy(A[1].mShares[0].data() + i*max_other_num, A_vector.mShares[0].data() + i*A_vector_i64cols + 2, max_other_num*sizeof(A_vector.mShares[0](i, 0)));
        std::memcpy(A[1].mShares[1].data() + i*max_other_num, A_vector.mShares[1].data() + i*A_vector_i64cols + 2, max_other_num*sizeof(A_vector.mShares[1](i, 0)));
        std::memcpy(A[0].mShares[0].data() + i*key_num, A_vector.mShares[0].data() + i*A_vector_i64cols + 2 + max_other_num, key_num*sizeof(A_vector.mShares[0](i, 0)));
        std::memcpy(A[0].mShares[1].data() + i*key_num, A_vector.mShares[1].data() + i*A_vector_i64cols + 2 + max_other_num, key_num*sizeof(A_vector.mShares[1](i, 0)));
    }
    // t2 = std::chrono::high_resolution_clock::now();
    // double step3_fill_down_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"oblivious_expand_step3_fill_down", step3_fill_down_time});
    // auto total_end = std::chrono::high_resolution_clock::now();
    // double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    // timing_results.push_back({"oblivious_expand_total", total_time});
    // if (pIdx == 0) {
    //     std::string filename = "./join_timing_results_role0.txt";
    //     std::ofstream outFile(filename, std::ios::app);
    //     if (outFile.is_open()) {
    //         outFile << "--- oblivious_expand Internal Timing (Role " << pIdx << ") ---" << std::endl;
    //         outFile << std::fixed << std::setprecision(3);
    //         for(const auto& result : timing_results) {
    //             outFile << "  " << result.first << ": " << result.second << " ms" << std::endl;
    //         }
    //         outFile << "----------------------------------------" << std::endl;
    //         outFile << std::endl;
    //         outFile.close();
    //     }
    // }

    return ;
}

/*
void oblivious_expand(int pIdx, std::vector<sbMatrix> &T, std::vector<sbMatrix> &A, i64 tid,
    i64Matrix &s_plain,
    aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    
    // 时间测量变量
    auto total_start = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<std::string, double>> timing_results;
    
    int n = T[0].rows();
    int key_num = T[0].i64Cols();
    int max_other_num = T[1].i64Cols();

    //step 1: 计算f(x)和s
    auto t1 = std::chrono::high_resolution_clock::now();
    //alpha的值64bit
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

    // //f(x)
    // sbMatrix fx(n, 64);
    // sbMatrix s(1, 64);
    // bool_init_false(pIdx, s);
    // for(size_t i=0; i<n; i++){
    //     sbMatrix f_i(1, 64), g_i(1, 64), zero_i(1, 64);
    //     sbMatrix flag_i(1, 1);
    //     bool_init_false(pIdx, zero_i);
    //     flag_i.mShares[0](0, 0) = flag.mShares[0](i, 0);
    //     flag_i.mShares[1](0, 0) = flag.mShares[1](i, 0);
       
    //     bool_cipher_selector(pIdx, flag_i, zero_i, s, f_i, enc, eval, runtime);

    //     fx.mShares[0](i, 0) = f_i.mShares[0](0, 0);
    //     fx.mShares[1](i, 0) = f_i.mShares[1](0, 0);
    //     g_i.mShares[0](0, 0) = g.mShares[0](i, 0);
    //     g_i.mShares[1](0, 0) = g.mShares[1](i, 0);
        
    //     bool_cipher_add(pIdx, s, g_i, s, enc, eval, runtime);
    // }

    //optimized:先求前缀和载selector
    

    s_plain.resize(1, 1);
    enc.revealAll(runtime, s, s_plain).get();
    int m = s_plain(0,0);
    auto t2 = std::chrono::high_resolution_clock::now();
    double step1_compute_fx_s_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"oblivious_expand_step1_compute_fx_s", step1_compute_fx_s_time});

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
    // i64Matrix flag_plain(n, 1);
    // enc.revealAll(runtime, flag, flag_plain).get();
    // if(pIdx == 0){
    //     std::cout << "flag_plain: " << std::endl;
    //     for(size_t i=0; i<n; i++){
    //         std::cout << flag_plain(i, 0) << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // //-----correct

    //step 2: obli_distribute
    t1 = std::chrono::high_resolution_clock::now();
    sbMatrix A_vector(m, 64*key_num+64*max_other_num+64+64);
    sbMatrix flag_sorted_auged(m, 1);
    oblivious_distribute(pIdx, T, flag, fx, s_plain, A_vector, flag_sorted_auged, enc, eval, runtime);
    int A_vector_i64cols = A_vector.i64Cols();
    t2 = std::chrono::high_resolution_clock::now();
    double step2_obli_distribute_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"oblivious_expand_step2_obli_distribute", step2_obli_distribute_time});

    //step 3:fill down
    t1 = std::chrono::high_resolution_clock::now();
    sbMatrix px_vector(1, 64*key_num+64*max_other_num+64+64);
    sbMatrix px(1, 64);
    bool_init_i64(pIdx, std::numeric_limits<i64>::max(), px, enc, eval, runtime);
    bool_cols_expand(px, px_vector);

    sbMatrix one_share(m, 1);
    bool_init_i64(pIdx, 1, one_share, enc, eval, runtime);
    sbMatrix cond(m, 1);
    bool_cipher_eq(pIdx, flag_sorted_auged, one_share, cond, enc, eval, runtime);
    

    for(size_t i=0; i<m; i++){
        sbMatrix cond_i(1, 1);
        cond_i.mShares[0](0, 0) = cond.mShares[0](i, 0);
        cond_i.mShares[1](0, 0) = cond.mShares[1](i, 0);

        //TODO
        // sbMatrix one_share(1, 1);
        // bool_init_i64(pIdx, 1, one_share, enc, eval, runtime);

        // sbMatrix cond(1, 1);
        // bool_cipher_eq(pIdx, flag_i, one_share, cond, enc, eval, runtime);
        //----放到外面

        sbMatrix A_k_i(1, 64*key_num+64*max_other_num+64+64), A_k_i_new(1, 64*key_num+64*max_other_num+64+64);
        for(size_t l=0; l<A_vector_i64cols; l++){
            A_k_i.mShares[0](0, l) = A_vector.mShares[0](i, l);
            A_k_i.mShares[1](0, l) = A_vector.mShares[1](i, l);
        }
        bool_cipher_selector(pIdx, cond_i, px_vector, A_k_i, A_k_i_new, enc, eval, runtime);
        for(size_t l=0; l<A_vector_i64cols; l++){
            A_vector.mShares[0](i, l) = A_k_i_new.mShares[0](0, l);
            A_vector.mShares[1](i, l) = A_k_i_new.mShares[1](0, l);
            px_vector.mShares[0](0, l) = A_k_i_new.mShares[0](0, l);
            px_vector.mShares[1](0, l) = A_k_i_new.mShares[1](0, l);
        }

    }

    //将A_vector拆分成A[k]
    A.resize(4);
    A[0].resize(m, 64*key_num);
    A[1].resize(m, 64*max_other_num);
    A[2].resize(m, 64);
    A[3].resize(m, 64);

    for(size_t i=0; i<m; i++){
        for(size_t l=0; l<key_num; l++){
            A[0].mShares[0](i, l) = A_vector.mShares[0](i, l);
            A[0].mShares[1](i, l) = A_vector.mShares[1](i, l);
        }
        for(size_t l=key_num; l<key_num + max_other_num; l++){
            A[1].mShares[0](i, l-key_num) = A_vector.mShares[0](i, l);
            A[1].mShares[1](i, l-key_num) = A_vector.mShares[1](i, l);
        }
        A[2].mShares[0](i, 0) = A_vector.mShares[0](i, key_num + max_other_num);
        A[2].mShares[1](i, 0) = A_vector.mShares[1](i, key_num + max_other_num);
        A[3].mShares[0](i, 0) = A_vector.mShares[0](i, key_num + max_other_num + 1);
        A[3].mShares[1](i, 0) = A_vector.mShares[1](i, key_num + max_other_num + 1);
    }
    t2 = std::chrono::high_resolution_clock::now();
    double step3_fill_down_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"oblivious_expand_step3_fill_down", step3_fill_down_time});

    // 计算总时间
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    timing_results.push_back({"oblivious_expand_total", total_time});

    // 将时间测量结果写入文件（仅 role 0 写入）
    if (pIdx == 0) {
        std::string filename = "./join_timing_results_role0.txt";
        std::ofstream outFile(filename, std::ios::app);
        if (outFile.is_open()) {
            outFile << "--- oblivious_expand Internal Timing (Role " << pIdx << ") ---" << std::endl;
            outFile << std::fixed << std::setprecision(3);
            for(const auto& result : timing_results) {
                outFile << "  " << result.first << ": " << result.second << " ms" << std::endl;
            }
            outFile << "----------------------------------------" << std::endl;
            outFile << std::endl;
            outFile.close();
        }
    }

    // //DEBUG
    // for(int l=0; l< key_num  ; l++){
    //     i64Matrix j_plain(m, 1);
    //     sbMatrix j_l(m, 64);
    //     for(int i=0; i<m; i++){
    //         j_l.mShares[0](i, 0) = A[0].mShares[0](i, l);
    //         j_l.mShares[1](i, 0) = A[0].mShares[1](i, l);
    //     }
    //     enc.revealAll(runtime, j_l, j_plain).get();
    //     if(pIdx ==0){
    //         std::cout<< "j_"<<l<<": "<<std::endl;
    //         for(size_t i=0; i<m; i++){
    //             std::cout <<j_plain(i, 0) << " ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }
    // for(int l=0; l< max_other_num  ; l++){
    //     i64Matrix d_plain(m, 1);
    //     sbMatrix d_l(m, 64);
    //     for(int i=0; i<m; i++){
    //         d_l.mShares[0](i, 0) = A[1].mShares[0](i, l);
    //         d_l.mShares[1](i, 0) = A[1].mShares[1](i, l);
    //     }
    //     enc.revealAll(runtime, d_l, d_plain).get();
    //     if(pIdx ==0){
    //         std::cout<< "d_"<<l<<": "<<std::endl;
    //         for(size_t i=0; i<m; i++){
    //             std::cout <<d_plain(i, 0) << " ";
    //         }
    //         std::cout<<std::endl;
    //     }
    // }
    // i64Matrix alpha_1_plain(m, 1);
    // sbMatrix alpha_1_f(m, 64);
    // for(int i=0; i<m; i++){
    //     alpha_1_f.mShares[0](i, 0) = A[2].mShares[0](i, 0);
    //     alpha_1_f.mShares[1](i, 0) = A[2].mShares[1](i, 0);
    // }
    // enc.revealAll(runtime, alpha_1_f, alpha_1_plain).get();
    // if(pIdx ==0){
    //     std::cout<< "alpha_1_f: "<<std::endl;
    //     for(size_t i=0; i<m; i++){
    //         std::cout <<alpha_1_plain(i, 0) << " ";
    //     }
    //     std::cout<<std::endl;
    // }
    // i64Matrix alpha_2_plain(m, 1);
    // sbMatrix alpha_2_f(m, 64);
    // for(int i=0; i<m; i++){
    //     alpha_2_f.mShares[0](i, 0) = A[3].mShares[0](i, 0);
    //     alpha_2_f.mShares[1](i, 0) = A[3].mShares[1](i, 0);
    // }
    // enc.revealAll(runtime, alpha_2_f, alpha_2_plain).get();
    // if(pIdx ==0){
    //     std::cout<< "alpha_2_f: "<<std::endl;
    //     for(size_t i=0; i<m; i++){
    //         std::cout <<alpha_2_plain(i, 0) << " ";
    //     }
    //     std::cout<<std::endl;
    // }
    // //-----correct: after fill down


    return ;
}
*/

//optimized:permsort(si) 
void oblivious_distribute(int pIdx, std::vector<si64Matrix> &T_prime, sbMatrix &flag, si64Matrix &fx, i64Matrix &s_plain, 
    sbMatrix &A_vector, sbMatrix &flag_sorted_auged,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    // // 时间测量变量（按需保留，当前禁用）
    // auto total_start = std::chrono::high_resolution_clock::now();
    // std::vector<std::pair<std::string, double>> timing_results;
    
    int n = T_prime[0].rows();
    i64 m = s_plain(0, 0);

    int key_num = T_prime[0].cols();
    int max_other_num = T_prime[1].cols();

    //step 1: sort T_prime
    // auto t1 = std::chrono::high_resolution_clock::now();

    //flag->si64
    si64Matrix flag_si64(n, 1);
    bool2arith(pIdx, flag, flag_si64, enc, eval, runtime);

    //sort by:flag, fx, j, d, alpha1, alpha2 
    int entry_key_cols = 1+1+key_num+max_other_num+1+1;
    si64Matrix entry_key_si(n, entry_key_cols);
    //flag列
    entry_key_si.mShares[0].col(entry_key_cols-1) = flag_si64.mShares[0].col(0);
    entry_key_si.mShares[1].col(entry_key_cols-1) = flag_si64.mShares[1].col(0);
    //fx列
    entry_key_si.mShares[0].col(entry_key_cols-2) = fx.mShares[0].col(0);
    entry_key_si.mShares[1].col(entry_key_cols-2) = fx.mShares[1].col(0);
    //j列
    entry_key_si.mShares[0].block(0, 2+max_other_num, n, key_num) = T_prime[0].mShares[0].block(0, 0, n, key_num);
    entry_key_si.mShares[1].block(0, 2+max_other_num, n, key_num) = T_prime[0].mShares[1].block(0, 0, n, key_num);
    //d列
    entry_key_si.mShares[0].block(0, 2, n, max_other_num) = T_prime[1].mShares[0].block(0, 0, n, max_other_num);
    entry_key_si.mShares[1].block(0, 2, n, max_other_num) = T_prime[1].mShares[1].block(0, 0, n, max_other_num);
    //alpha1列
    entry_key_si.mShares[0].col(1) = T_prime[2].mShares[0].col(0);
    entry_key_si.mShares[1].col(1) = T_prime[2].mShares[1].col(0);
    //alpha2列
    entry_key_si.mShares[0].col(0) = T_prime[3].mShares[0].col(0);
    entry_key_si.mShares[1].col(0) = T_prime[3].mShares[1].col(0);

    
    si64Matrix perm(n, 1);
    //perm according flag,fx
    si64Matrix entry_key_perm_si(n, 2);
    entry_key_perm_si.mShares[0].col(0) = entry_key_si.mShares[0].block(0, entry_key_cols-2, n, 1);
    entry_key_perm_si.mShares[0].col(1) = entry_key_si.mShares[0].block(0, entry_key_cols-1, n, 1);
    entry_key_perm_si.mShares[1].col(0) = entry_key_si.mShares[1].block(0, entry_key_cols-2, n, 1);
    entry_key_perm_si.mShares[1].col(1) = entry_key_si.mShares[1].block(0, entry_key_cols-1, n, 1);

    //genPerm(pIdx, entry_key_si, perm, enc, eval, runtime);
    genPerm(pIdx, entry_key_perm_si, perm, enc, eval, runtime);
    si64Matrix entry_key_sorted_si(n, entry_key_cols);
    applyPerm(pIdx, perm, entry_key_si, entry_key_sorted_si, enc, eval, runtime);
    // auto t2 = std::chrono::high_resolution_clock::now();
    // double step1_sort_T_prime_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"oblivious_distribute_step1_sort_T_prime", step1_sort_T_prime_time});



    //step 2 填充空值得到A_vector, fx_sorted_auged, flag_sorted_auged
    //flag_sorted_auged_si(m, 1), fx_sorted_auged_si(m, 1);
    // t1 = std::chrono::high_resolution_clock::now();
    sbMatrix entry_key_sorted_sb(n, 64*entry_key_cols);
    arith2bool(pIdx, entry_key_sorted_si, entry_key_sorted_sb, enc, eval, runtime);
    
    //A_vector.resize(m, 64*entry_key_cols);
    sbMatrix A_vector_concat(m, 64*entry_key_cols);
    if(m >= n){
        sbMatrix null_share_m(m-n, (key_num+max_other_num+2)*64), zero_m(m-n, 64), one_share_m(m-n, 64);
        set_const_share_bool(pIdx, 0, zero_m, enc, eval, runtime);
        set_const_share_bool(pIdx, 1, one_share_m, enc, eval, runtime);
        set_const_share_bool(pIdx, std::numeric_limits<i64>::max(), null_share_m, enc, eval, runtime);
        for(size_t s=0; s<2; s++){
            // copy first n rows (same stride: entry_key_cols)
            std::memcpy(A_vector_concat.mShares[s].data(), entry_key_sorted_sb.mShares[s].data(), n * entry_key_cols * sizeof(i64));
            for(size_t i=n; i<m; i++){
                for(size_t j=0; j<(key_num+max_other_num+2); j++){
                    A_vector_concat.mShares[s](i, j) = null_share_m.mShares[s](i-n, j);
                }
                A_vector_concat.mShares[s](i, key_num+max_other_num+2) = zero_m.mShares[s](i-n, 0);
                A_vector_concat.mShares[s](i, key_num+max_other_num+3) = one_share_m.mShares[s](i-n, 0);
            }
        }
    }
    else{
        for(size_t s=0; s<2; s++){
            // copy first m rows (same stride: entry_key_cols)
            std::memcpy(A_vector_concat.mShares[s].data(), entry_key_sorted_sb.mShares[s].data(), m * entry_key_cols * sizeof(i64));
        }
    }
    // t2 = std::chrono::high_resolution_clock::now();
    // double step2_fill_null_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"oblivious_distribute_step2_fill_null", step2_fill_null_time});

    //step 3: distribute loop
    // t1 = std::chrono::high_resolution_clock::now();
    // 计算 j = 2^⌈log₂(m)⌉ - 1
    i64 step_size;
    i64 power = static_cast<i64>(std::ceil(std::log2(m)))-1;
    step_size = 1ULL << power;

    i64Matrix idx_plain(m, 1);
    for(size_t i=0; i<m; i++){
        idx_plain(i, 0) = i;
    }
    sbMatrix idx(m, 64);
    if(pIdx == 0){
        enc.localBinMatrix(runtime, idx_plain, idx).get();
    }else{
        enc.remoteBinMatrix(runtime, idx).get();
    }

    while(step_size >= 1){

        sbMatrix i_plus_j(m-step_size, 64), fx_vector(m-step_size, 64);
     
        for(size_t s=0; s<2; s++){
            // idx is single-column, contiguous: shift-copy by step_size-1
            std::memcpy(i_plus_j.mShares[s].data(), idx.mShares[s].data() + (step_size-1), (m-step_size) * sizeof(i64));
            for(size_t i=0; i<m-step_size; i++){
                fx_vector.mShares[s](i,0) = A_vector_concat.mShares[s](i, key_num+max_other_num+2);
            }
        }
        sbMatrix cond1(m-step_size, 1), cond1_concat(2, 1);

        bool_cipher_lt(pIdx, i_plus_j, fx_vector, cond1, enc, eval, runtime);

        // sbMatrix A_true_value(2,64*entry_key_cols);
        // sbMatrix A_false_value(2,64*entry_key_cols);
        // sbMatrix A_new_value(2,64*entry_key_cols);

        //for(i64 i = m-1-step_size; i >= 0; i--){
        
        //     cond1_concat.mShares[0](0, 0) = cond1.mShares[0](i, 0);
        //     cond1_concat.mShares[1](0, 0) = cond1.mShares[1](i, 0);
        //     cond1_concat.mShares[0](1, 0) = cond1.mShares[0](i, 0);
        //     cond1_concat.mShares[1](1, 0) = cond1.mShares[1](i, 0);

        //     for(size_t s=0; s<2; s++){
        //         for(size_t l=0; l<entry_key_cols; l++){
        //             // cond=1(swap): row0->A[i] 取 A[i+step_size], row1->A[i+step_size] 取 A[i]
        //             A_true_value.mShares[s](0, l) = A_vector_concat.mShares[s](i+step_size, l);
        //             A_true_value.mShares[s](1, l) = A_vector_concat.mShares[s](i, l);
        //             // cond=0(no swap): row0->A[i] 取 A[i], row1->A[i+step_size] 取 A[i+step_size]
        //             A_false_value.mShares[s](0, l) = A_vector_concat.mShares[s](i, l);
        //             A_false_value.mShares[s](1, l) = A_vector_concat.mShares[s](i+step_size, l);
        //         }
        //     }

        //     bool_cipher_selector(pIdx, cond1_concat, A_true_value, A_false_value, A_new_value, enc, eval, runtime);

         
        //     //更新值
        //     for(size_t s=0; s<2; s++){
        //         for(size_t l=0; l<entry_key_cols; l++){
        //             A_vector_concat.mShares[s](i, l) = A_new_value.mShares[s](0, l);
        //         }
        //     }
        //     for(size_t s=0; s<2; s++){
        //         for(size_t l=0; l<entry_key_cols; l++){
        //             A_vector_concat.mShares[s](i+step_size, l) = A_new_value.mShares[s](1, l);
        //         }
        //     }
        // }

        const i64 j = step_size;
        const i64 max_t = (m - 1) / j - 1;
        // 每个 t 对应一个“无冲突子层”
        // 同一子层里要交换的pair: (r + t*j, r + (t+1)*j), r = 0,1,...,pair_cnt-1
        for (i64 t = max_t; t >= 0; --t) {
            //r + (t+1)*j < m
            const i64 pair_cnt = std::min<i64>(j, m - (t + 1) * j);
            if (pair_cnt <= 0) continue;
            // 一个 pair 需要 2 行：
            // row 2*k     -> 目标写回 A[i]
            // row 2*k + 1 -> 目标写回 A[i + j]
            sbMatrix cond_batch(2 * pair_cnt, 1);
            sbMatrix A_true_value(2 * pair_cnt, 64 * entry_key_cols);
            sbMatrix A_false_value(2 * pair_cnt, 64 * entry_key_cols);
            sbMatrix A_new_value(2 * pair_cnt, 64 * entry_key_cols);
            
            for (i64 r = 0; r < pair_cnt; ++r) {
                const i64 i = r + t * j;
                const i64 ip = i + j;
                
                cond_batch.mShares[0](2 * r, 0)     = cond1.mShares[0](i, 0);
                cond_batch.mShares[1](2 * r, 0)     = cond1.mShares[1](i, 0);
                cond_batch.mShares[0](2 * r + 1, 0) = cond1.mShares[0](i, 0);
                cond_batch.mShares[1](2 * r + 1, 0) = cond1.mShares[1](i, 0);
                for (size_t s = 0; s < 2; ++s) {
                    // cond = 1 -> swap
                    std::memcpy(
                        A_true_value.mShares[s].data() + (2 * r) * entry_key_cols,
                        A_vector_concat.mShares[s].data() + ip * entry_key_cols,
                        entry_key_cols * sizeof(i64));
                    std::memcpy(
                        A_true_value.mShares[s].data() + (2 * r + 1) * entry_key_cols,
                        A_vector_concat.mShares[s].data() + i * entry_key_cols,
                        entry_key_cols * sizeof(i64));
                    // cond = 0 -> no swap
                    std::memcpy(
                        A_false_value.mShares[s].data() + (2 * r) * entry_key_cols,
                        A_vector_concat.mShares[s].data() + i * entry_key_cols,
                        entry_key_cols * sizeof(i64));
                    std::memcpy(
                        A_false_value.mShares[s].data() + (2 * r + 1) * entry_key_cols,
                        A_vector_concat.mShares[s].data() + ip * entry_key_cols,
                        entry_key_cols * sizeof(i64));
                }
            }
            // 当前无冲突子层，一次 selector
            bool_cipher_selector(pIdx, cond_batch, A_true_value, A_false_value, A_new_value, enc, eval, runtime);
           
            for (i64 r = 0; r < pair_cnt; ++r) {
                const i64 i = r + t * j;
                const i64 ip = i + j;
                for (size_t s = 0; s < 2; ++s) {
                    std::memcpy(A_vector_concat.mShares[s].data() + i  * entry_key_cols, A_new_value.mShares[s].data() + (2*r)   * entry_key_cols, entry_key_cols * sizeof(i64));
                    std::memcpy(A_vector_concat.mShares[s].data() + ip * entry_key_cols, A_new_value.mShares[s].data() + (2*r+1) * entry_key_cols, entry_key_cols * sizeof(i64));
                }
            }
        }
        
        step_size = step_size / 2;
    }

    A_vector.resize(m, 64*(2+key_num+max_other_num));
    flag_sorted_auged.resize(m, 1);
    for(size_t s=0; s<2; s++){
        for(size_t i=0; i<m; i++){
            for(size_t l=0; l<(2+key_num+max_other_num); l++){
                A_vector.mShares[s](i, l) = A_vector_concat.mShares[s](i, l);
            }
            flag_sorted_auged.mShares[s](i, 0) = A_vector_concat.mShares[s](i, key_num+max_other_num+3);
        }
    }
    // t2 = std::chrono::high_resolution_clock::now();
    // double step3_distribute_loop_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"oblivious_distribute_step3_distribute_loop", step3_distribute_loop_time});
    // auto total_end = std::chrono::high_resolution_clock::now();
    // double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    // timing_results.push_back({"oblivious_distribute_total", total_time});
    // if (pIdx == 0) {
    //     std::string filename = "./join_timing_results_role0.txt";
    //     std::ofstream outFile(filename, std::ios::app);
    //     if (outFile.is_open()) {
    //         outFile << "--- oblivious_distribute Internal Timing (Role " << pIdx << ") ---" << std::endl;
    //         outFile << std::fixed << std::setprecision(3);
    //         for(const auto& result : timing_results) {
    //             outFile << "  " << result.first << ": " << result.second << " ms" << std::endl;
    //         }
    //         outFile << "----------------------------------------" << std::endl;
    //         outFile << std::endl;
    //         outFile.close();
    //     }
    // }
        
        return;
    }
/*
void oblivious_distribute(int pIdx, std::vector<sbMatrix> &T_prime, sbMatrix &flag, sbMatrix &fx, i64Matrix &s_plain, 
    sbMatrix &A_vector, sbMatrix &flag_sorted_auged,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    // 时间测量变量
    auto total_start = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<std::string, double>> timing_results;
    
    int n = T_prime[0].rows();
    i64 m = s_plain(0, 0);

    int key_num = T_prime[0].i64Cols();
    int max_other_num = T_prime[1].i64Cols();

    //step 1: sort T_prime(转成了bool进行拼接排序)
    auto t1 = std::chrono::high_resolution_clock::now();

    //flag(1bit),j(key_num*64bit),d(max_other_num*64bit),fx(64bit),alpha1(64bit),alpha2(64bit)
    int entry_key_bitcount = 64+64*key_num+64*max_other_num+64+64+64;
    int entry_key_cols = 1+key_num+max_other_num+1+1+1;
    sbMatrix entry_key(n, entry_key_bitcount);
    for(size_t i=0; i<n; i++){
        //flag
        entry_key.mShares[0](i, entry_key_cols-1) = flag.mShares[0](i, 0);
        entry_key.mShares[1](i, entry_key_cols-1) = flag.mShares[1](i, 0);
        //j
        for(int l=entry_key_cols-2; l>=3+max_other_num; l--){
            entry_key.mShares[0](i, l) = T_prime[0].mShares[0](i, entry_key_cols-2-l);
            entry_key.mShares[1](i, l) = T_prime[0].mShares[1](i, entry_key_cols-2-l);
        }
        //d
        for(int l=2+max_other_num ; l>=3; l--){
            entry_key.mShares[0](i, l) = T_prime[1].mShares[0](i, 2+max_other_num-l);
            entry_key.mShares[1](i, l) = T_prime[1].mShares[1](i, 2+max_other_num-l);
        }
        //fx
        entry_key.mShares[0](i, 2) = fx.mShares[0](i, 0);
        entry_key.mShares[1](i, 2) = fx.mShares[1](i, 0);
        //alpha1
        entry_key.mShares[0](i, 1) = T_prime[2].mShares[0](i, 0);
        entry_key.mShares[1](i, 1) = T_prime[2].mShares[1](i, 0);
        //alpha2
        entry_key.mShares[0](i, 0) = T_prime[3].mShares[0](i, 0);
        entry_key.mShares[1](i, 0) = T_prime[3].mShares[1](i, 0);
    }

    sbMatrix entry_key_sorted(n, entry_key_bitcount);
    odd_even_merge_sort(entry_key, entry_key_sorted, pIdx, enc, eval, runtime);

    sbMatrix j(n, 64*key_num), d(n, 64*max_other_num), fx_prime(n, 64), flag_prime(n, 1);
    sbMatrix alpha_1(n, 64), alpha_2(n, 64);
    for(size_t i=0; i< n; i++){
        //flag
        flag_prime.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, entry_key_cols-1);
        flag_prime.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, entry_key_cols-1);
        //j
        for(int l=entry_key_cols-2; l>=3+max_other_num; l--){
            j.mShares[0](i, entry_key_cols-2-l) = entry_key_sorted.mShares[0](i, l);
            j.mShares[1](i, entry_key_cols-2-l) = entry_key_sorted.mShares[1](i, l);
        }
        //d
        for(int l=2+max_other_num ; l>=3; l--){
            d.mShares[0](i, 2+max_other_num-l) = entry_key_sorted.mShares[0](i, l);
            d.mShares[1](i, 2+max_other_num-l) = entry_key_sorted.mShares[1](i, l);
        }
        //fx
        fx_prime.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, 2);
        fx_prime.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, 2);
        //alpha1
        alpha_1.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, 1);
        alpha_1.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, 1);
        //alpha2
        alpha_2.mShares[0](i, 0) = entry_key_sorted.mShares[0](i, 0);
        alpha_2.mShares[1](i, 0) = entry_key_sorted.mShares[1](i, 0);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double step1_sort_T_prime_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"oblivious_distribute_step1_sort_T_prime", step1_sort_T_prime_time});

    std::vector<sbMatrix> T_sorted(4);
    T_sorted[0] = j ;
    T_sorted[1] = d;
    T_sorted[2] = alpha_1;
    T_sorted[3] = alpha_2;

    sbMatrix fx_sorted_auged(m, 64);
    flag_sorted_auged.resize(m, 1);
    std::vector<sbMatrix> A(4);
    A[0].resize(m, 64*key_num);
    A[1].resize(m, 64*max_other_num);
    A[2].resize(m, 64);
    A[3].resize(m, 64);
    
    //step 2 填充空值得到A, fx_sorted_auged
    t1 = std::chrono::high_resolution_clock::now();
    if(m >= n){
        // for(size_t i=0; i<n; i++){
        //     for(size_t j=0; j<key_num; j++){
        //         A[0].mShares[0](i, j) = T_sorted[0].mShares[0](i, j);
        //         A[0].mShares[1](i, j) = T_sorted[0].mShares[1](i, j);
        //     }
        //     for(size_t j=0; j<max_other_num; j++){
        //         A[1].mShares[0](i, j) = T_sorted[1].mShares[0](i, j);
        //         A[1].mShares[1](i, j) = T_sorted[1].mShares[1](i, j);
        //     }
        // }
        std::memcpy(A[0].mShares[0].data(), T_sorted[0].mShares[0].data(), n * key_num * sizeof(T_sorted[0].mShares[0](0, 0)));
        std::memcpy(A[0].mShares[1].data(), T_sorted[0].mShares[1].data(), n * key_num * sizeof(T_sorted[0].mShares[1](0, 0)));
        std::memcpy(A[1].mShares[0].data(), T_sorted[1].mShares[0].data(), n * max_other_num * sizeof(T_sorted[1].mShares[0](0, 0)));
        std::memcpy(A[1].mShares[1].data(), T_sorted[1].mShares[1].data(), n * max_other_num * sizeof(T_sorted[1].mShares[1](0, 0)));
        for(size_t i=2; i<4; i++){
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
        
     
        // for(size_t i=n; i<m; i++){
        //     for(size_t j=0; j<key_num; j++){
        //         A[0].mShares[0](i, j) = null_share_m.mShares[0](i-n, 0);
        //         A[0].mShares[1](i, j) = null_share_m.mShares[1](i-n, 0);
        //     }
        //     for(size_t j=0; j<max_other_num; j++){
        //         A[1].mShares[0](i, j) = null_share_m.mShares[0](i-n, 0);
        //         A[1].mShares[1](i, j) = null_share_m.mShares[1](i-n, 0);
        //     }
        // }
        std::memcpy(A[0].mShares[0].data() + n * key_num, null_share_m.mShares[0].data(), (m-n) * key_num * sizeof(null_share_m.mShares[0](0, 0)));
        std::memcpy(A[0].mShares[1].data() + n * key_num, null_share_m.mShares[1].data(), (m-n) * key_num * sizeof(null_share_m.mShares[1](0, 0)));
        std::memcpy(A[1].mShares[0].data() + n * max_other_num, null_share_m.mShares[0].data(), (m-n) * max_other_num * sizeof(null_share_m.mShares[0](0, 0)));
        std::memcpy(A[1].mShares[1].data() + n * max_other_num, null_share_m.mShares[1].data(), (m-n) * max_other_num * sizeof(null_share_m.mShares[1](0, 0)));

        for(size_t i=2; i<4; i++){
            std::memcpy(A[i].mShares[0].data() + n, null_share_m.mShares[0].data(), (m-n) * sizeof(null_share_m.mShares[0](0, 0)));
            std::memcpy(A[i].mShares[1].data() + n, null_share_m.mShares[1].data(), (m-n) * sizeof(null_share_m.mShares[1](0, 0)));
        }
        std::memcpy(fx_sorted_auged.mShares[0].data() + n, zero_m.mShares[0].data(), (m-n) * sizeof(zero_m.mShares[0](0, 0)));
        std::memcpy(fx_sorted_auged.mShares[1].data() + n, zero_m.mShares[1].data(), (m-n) * sizeof(zero_m.mShares[1](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[0].data() + n, one_share_m.mShares[0].data(), (m-n) * sizeof(one_share_m.mShares[0](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[1].data() + n, one_share_m.mShares[1].data(), (m-n) * sizeof(one_share_m.mShares[1](0, 0)));
    }
    else{
        // for(size_t i=0; i<m; i++){
        //     for(size_t j=0; j<key_num; j++){
        //         A[0].mShares[0](i, j) = T_sorted[0].mShares[0](i, j);
        //         A[0].mShares[1](i, j) = T_sorted[0].mShares[1](i, j);
        //     }
        //     for(size_t j=0; j<max_other_num; j++){
        //         A[1].mShares[0](i, j) = T_sorted[1].mShares[0](i, j);
        //         A[1].mShares[1](i, j) = T_sorted[1].mShares[1](i, j);
        //     }
        // }
        std::memcpy(A[0].mShares[0].data(), T_sorted[0].mShares[0].data(), m * key_num * sizeof(T_sorted[0].mShares[0](0, 0)));
        std::memcpy(A[0].mShares[1].data(), T_sorted[0].mShares[1].data(), m * key_num * sizeof(T_sorted[0].mShares[1](0, 0)));
        std::memcpy(A[1].mShares[0].data(), T_sorted[1].mShares[0].data(), m * max_other_num * sizeof(T_sorted[1].mShares[0](0, 0)));
        std::memcpy(A[1].mShares[1].data(), T_sorted[1].mShares[1].data(), m * max_other_num * sizeof(T_sorted[1].mShares[1](0, 0)));
        for(size_t i=2; i<4; i++){
            std::memcpy(A[i].mShares[0].data(), T_sorted[i].mShares[0].data(), m * sizeof(T_sorted[i].mShares[0](0, 0)));
            std::memcpy(A[i].mShares[1].data(), T_sorted[i].mShares[1].data(), m * sizeof(T_sorted[i].mShares[1](0, 0)));
        }
        std::memcpy(fx_sorted_auged.mShares[0].data(), fx_prime.mShares[0].data(), m * sizeof(fx_prime.mShares[0](0, 0)));
        std::memcpy(fx_sorted_auged.mShares[1].data(), fx_prime.mShares[1].data(), m * sizeof(fx_prime.mShares[1](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[0].data(), flag_prime.mShares[0].data(), m * sizeof(flag_prime.mShares[0](0, 0)));
        std::memcpy(flag_sorted_auged.mShares[1].data(), flag_prime.mShares[1].data(), m * sizeof(flag_prime.mShares[1](0, 0)));
    }
    t2 = std::chrono::high_resolution_clock::now();
    double step2_fill_null_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"oblivious_distribute_step2_fill_null", step2_fill_null_time});

    //step 3: distribute loop
    t1 = std::chrono::high_resolution_clock::now();
    // 计算 j = 2^⌈log₂(m)⌉ - 1
    i64 step_size;
    i64 power = static_cast<i64>(std::ceil(std::log2(m)))-1;
    step_size = 1ULL << power;

    //拼接成一个大的vector, 为了后面批量的bool_cipher_selector
    A_vector.resize(m, 64*key_num+64*max_other_num+64+64);
    for(size_t i=0; i<m ; i++){
        for(size_t l=0; l<key_num; l++){
            A_vector.mShares[0](i, l) = A[0].mShares[0](i, l);
            A_vector.mShares[1](i, l) = A[0].mShares[1](i, l);
        }
        for(size_t l=key_num; l<key_num + max_other_num; l++){
            A_vector.mShares[0](i, l) = A[1].mShares[0](i, l-key_num);
            A_vector.mShares[1](i, l) = A[1].mShares[1](i, l-key_num);
        }
        A_vector.mShares[0](i, key_num + max_other_num) = A[2].mShares[0](i, 0);
        A_vector.mShares[1](i, key_num + max_other_num) = A[2].mShares[1](i, 0);
        A_vector.mShares[0](i, key_num + max_other_num + 1) = A[3].mShares[0](i, 0);
        A_vector.mShares[1](i, key_num + max_other_num + 1) = A[3].mShares[1](i, 0);    
    }

    while(step_size >= 1){
        //i+j放在这里？bool_init?
        sbMatrix i_plus_j(m-step_size, 64);
        i64Matrix i_plus_j_plain(m-step_size, 1);
        for(size_t i=0; i<m-step_size; i++){
            i_plus_j_plain(i, 0) = i + step_size - 1;
        }
        if(pIdx == 0){
            enc.localBinMatrix(runtime, i_plus_j_plain, i_plus_j).get();
        }else{
            enc.remoteBinMatrix(runtime, i_plus_j).get();
        }
        
        for(i64 i = m-1-step_size; i >= 0; i--){
            sbMatrix flag_i(1, 1), fx_i(1, 64);
            sbMatrix flag_i_plus_j(1, 1), fx_i_plus_j(1, 64);
            fx_i.mShares[0](0, 0) = fx_sorted_auged.mShares[0](i, 0);
            fx_i.mShares[1](0, 0) = fx_sorted_auged.mShares[1](i, 0);
            // flag_i.mShares[0](0, 0) = flag_sorted_auged.mShares[0](i, 0);
            // flag_i.mShares[1](0, 0) = flag_sorted_auged.mShares[1](i, 0);

            fx_i_plus_j.mShares[0](0, 0) = fx_sorted_auged.mShares[0](i+step_size, 0);
            fx_i_plus_j.mShares[1](0, 0) = fx_sorted_auged.mShares[1](i+step_size, 0);
            // flag_i_plus_j.mShares[0](0, 0) = flag_sorted_auged.mShares[0](i+step_size, 0);
            // flag_i_plus_j.mShares[1](0, 0) = flag_sorted_auged.mShares[1](i+step_size, 0);

            // sbMatrix i_plus_j(1, 64);
            // bool_init_i64(pIdx, i+step_size-1, i_plus_j, enc, eval, runtime);
            sbMatrix i_plus_j_cond(1, 64);
            i_plus_j_cond.mShares[0](0, 0) = i_plus_j.mShares[0](i, 0);
            i_plus_j_cond.mShares[1](0, 0) = i_plus_j.mShares[1](i, 0);

            //目标idx超过i+j则交换，否则值不变
            sbMatrix cond1(1, 1),cond1_concat(2,1);
            bool_cipher_lt(pIdx, i_plus_j_cond, fx_i, cond1, enc, eval, runtime);
            cond1_concat.mShares[0](0, 0) = cond1.mShares[0](0, 0);
            cond1_concat.mShares[1](0, 0) = cond1.mShares[1](0, 0);
            cond1_concat.mShares[0](1, 0) = cond1.mShares[0](0, 0);
            cond1_concat.mShares[1](1, 0) = cond1.mShares[1](0, 0);

            //A_vector + fx + flag
            sbMatrix A_true_value(2,64*key_num+64*max_other_num+64+64+64+64);
            sbMatrix A_false_value(2,64*key_num+64*max_other_num+64+64+64+64);
            sbMatrix A_new_value(2,64*key_num+64*max_other_num+64+64+64+64);
            int A_vector_i64cols = A_vector.i64Cols();
            for(size_t l=0; l<A_vector_i64cols; l++){
                A_true_value.mShares[0](0, l) = A_vector.mShares[0](i+step_size, l);
                A_true_value.mShares[1](0, l) = A_vector.mShares[1](i+step_size, l);
                A_true_value.mShares[0](1, l) = A_vector.mShares[0](i, l);
                A_true_value.mShares[1](1, l) = A_vector.mShares[1](i, l);
                A_false_value.mShares[0](0, l) = A_vector.mShares[0](i, l);
                A_false_value.mShares[1](0, l) = A_vector.mShares[1](i, l);
                A_false_value.mShares[0](1, l) = A_vector.mShares[0](i+step_size, l);
                A_false_value.mShares[1](1, l) = A_vector.mShares[1](i+step_size, l);
            }
            //fx
            A_true_value.mShares[0](0, A_vector_i64cols) = fx_i_plus_j.mShares[0](0, 0);
            A_true_value.mShares[1](0, A_vector_i64cols) = fx_i_plus_j.mShares[1](0, 0);
            A_true_value.mShares[0](1, A_vector_i64cols) = fx_i.mShares[0](0, 0);
            A_true_value.mShares[1](1, A_vector_i64cols) = fx_i.mShares[1](0, 0);
            A_false_value.mShares[0](0, A_vector_i64cols) = fx_i.mShares[0](0, 0);
            A_false_value.mShares[1](0, A_vector_i64cols) = fx_i.mShares[1](0, 0);
            A_false_value.mShares[0](1, A_vector_i64cols) = fx_i_plus_j.mShares[0](0, 0);
            A_false_value.mShares[1](1, A_vector_i64cols) = fx_i_plus_j.mShares[1](0, 0);
            //flag
            A_true_value.mShares[0](0, A_vector_i64cols + 1) = flag_sorted_auged.mShares[0](i+step_size, 0);
            A_true_value.mShares[1](0, A_vector_i64cols + 1) = flag_sorted_auged.mShares[1](i+step_size, 0);
            A_true_value.mShares[0](1, A_vector_i64cols + 1) = flag_sorted_auged.mShares[0](i, 0);
            A_true_value.mShares[1](1, A_vector_i64cols + 1) = flag_sorted_auged.mShares[1](i, 0);
            A_false_value.mShares[0](0, A_vector_i64cols + 1) = flag_sorted_auged.mShares[0](i, 0);
            A_false_value.mShares[1](0, A_vector_i64cols + 1) = flag_sorted_auged.mShares[1](i, 0);
            A_false_value.mShares[0](1, A_vector_i64cols + 1) = flag_sorted_auged.mShares[0](i+step_size, 0);
            A_false_value.mShares[1](1, A_vector_i64cols + 1) = flag_sorted_auged.mShares[1](i+step_size, 0);

            bool_cipher_selector(pIdx, cond1_concat, A_true_value, A_false_value, A_new_value, enc, eval, runtime);

            // sbMatrix A_k_i(1, 64*key_num+64*max_other_num+64+64), A_k_i_plus_j(1, 64*key_num+64*max_other_num+64+64);
            // sbMatrix A_k_i_new(1, 64*key_num+64*max_other_num+64+64), A_k_i_plus_new(1, 64*key_num+64*max_other_num+64+64);
            // int A_vector_i64cols = A_vector.i64Cols();
            // for(size_t l=0; l<A_vector_i64cols; l++){
            //     A_k_i.mShares[0](0, l) = A_vector.mShares[0](i, l);
            //     A_k_i.mShares[1](0, l) = A_vector.mShares[1](i, l);
            // }
            // for(size_t l=0; l<A_vector_i64cols; l++){
            //     A_k_i_plus_j.mShares[0](0, l) = A_vector.mShares[0](i+step_size, l);
            //     A_k_i_plus_j.mShares[1](0, l) = A_vector.mShares[1](i+step_size, l);
            // }
            
            // bool_cipher_selector(pIdx, cond1, A_k_i_plus_j, A_k_i, A_k_i_new, enc, eval, runtime);
            // bool_cipher_selector(pIdx, cond1, A_k_i, A_k_i_plus_j, A_k_i_plus_new, enc, eval, runtime);

            for(size_t l=0; l<A_vector_i64cols; l++){
                A_vector.mShares[0](i, l) = A_new_value.mShares[0](0, l);
                A_vector.mShares[1](i, l) = A_new_value.mShares[1](0, l);
            }
            for(size_t l=0; l<A_vector_i64cols; l++){
                A_vector.mShares[0](i+step_size, l) = A_new_value.mShares[0](1, l);
                A_vector.mShares[1](i+step_size, l) = A_new_value.mShares[1](1, l);
            }
            
            //---应该把下面的flag和fx都拼到A_vector后面再加两列
            // sbMatrix flag_i_new(1, 1), fx_i_new(1, 64);
            // bool_cipher_selector(pIdx, cond1, fx_i_plus_j, fx_i, fx_i_new, enc, eval, runtime);
            // bool_cipher_selector(pIdx, cond1, flag_i_plus_j, flag_i, flag_i_new, enc, eval, runtime);
            
            fx_sorted_auged.mShares[0](i, 0) = A_new_value.mShares[0](0, A_vector_i64cols);
            fx_sorted_auged.mShares[1](i, 0) = A_new_value.mShares[1](0, A_vector_i64cols);
            flag_sorted_auged.mShares[0](i, 0) = A_new_value.mShares[0](0, A_vector_i64cols + 1);
            flag_sorted_auged.mShares[1](i, 0) = A_new_value.mShares[1](0, A_vector_i64cols + 1);

            // sbMatrix flag_i_plus_new(1, 1), fx_i_plus_new(1, 64);
            // bool_cipher_selector(pIdx, cond1, fx_i, fx_i_plus_j, fx_i_plus_new, enc, eval, runtime);
            // bool_cipher_selector(pIdx, cond1, flag_i, flag_i_plus_j, flag_i_plus_new, enc, eval, runtime);
         
            fx_sorted_auged.mShares[0](i+step_size, 0) =  A_new_value.mShares[0](1, A_vector_i64cols);
            fx_sorted_auged.mShares[1](i+step_size, 0) = A_new_value.mShares[1](1, A_vector_i64cols);
            flag_sorted_auged.mShares[0](i+step_size, 0) = A_new_value.mShares[0](1, A_vector_i64cols + 1);
            flag_sorted_auged.mShares[1](i+step_size, 0) = A_new_value.mShares[1](1, A_vector_i64cols + 1);
        }
        step_size = step_size / 2;
    }
    t2 = std::chrono::high_resolution_clock::now();
    double step3_distribute_loop_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"oblivious_distribute_step3_distribute_loop", step3_distribute_loop_time});

    // 计算总时间
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    timing_results.push_back({"oblivious_distribute_total", total_time});

    // 将时间测量结果写入文件（仅 role 0 写入）
    if (pIdx == 0) {
        std::string filename = "./join_timing_results_role0.txt";
        std::ofstream outFile(filename, std::ios::app);
        if (outFile.is_open()) {
            outFile << "--- oblivious_distribute Internal Timing (Role " << pIdx << ") ---" << std::endl;
            outFile << std::fixed << std::setprecision(3);
            for(const auto& result : timing_results) {
                outFile << "  " << result.first << ": " << result.second << " ms" << std::endl;
            }
            outFile << "----------------------------------------" << std::endl;
            outFile << std::endl;
            outFile.close();
        }
    }

    return ;
}
*/

//optimized:permsort(si)
void align_table(int pIdx, std::vector<sbMatrix> &T, std::vector<si64Matrix> &T_aligned,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    // // 时间测量变量（按需保留，当前禁用）
    // auto total_start = std::chrono::high_resolution_clock::now();
    // std::vector<std::pair<std::string, double>> timing_results;
    
    int m=T[0].rows();
    int key_num = T[0].i64Cols();
    int max_other_num = T[1].i64Cols();

    //step 1: set e.ii
    // auto t1 = std::chrono::high_resolution_clock::now();

        //step 1.1 same_attr
    sbMatrix false_matrix(1,1);
    bool_init_false(pIdx, false_matrix);
    sbMatrix same_attr_partial(m-1, 1);
    compare_consecutive_rows_bool(pIdx, T[0], same_attr_partial, enc, eval, runtime);
    sbMatrix same_attr(m, 1);
    same_attr.mShares[0](0, 0) = false_matrix.mShares[0](0, 0);
    same_attr.mShares[1](0, 0) = false_matrix.mShares[1](0, 0);
    std::memcpy(same_attr.mShares[0].data() + 1, same_attr_partial.mShares[0].data(), (m-1) * sizeof(same_attr_partial.mShares[0](0, 0)));
    std::memcpy(same_attr.mShares[1].data() + 1, same_attr_partial.mShares[1].data(), (m-1) * sizeof(same_attr_partial.mShares[1](0, 0)));
    
    //DEBUG print alpha_plain
    // i64Matrix alpha1_plain(m, 1);
    // enc.revealAll(runtime, T[2], alpha1_plain).get();
    // if (pIdx == 0) {
    //     std::ofstream debug_out("./debug_alpha1_plain.txt", std::ios::app);
    //     if (debug_out.is_open()) {
    //         debug_out << "--- align_table alpha1_plain (rows=" << m << ") ---" << std::endl;
    //         for (int i = 0; i < m; ++i) {
    //             debug_out << alpha1_plain(i, 0) << std::endl;
    //         }
    //         debug_out << std::endl;
    //         debug_out.close();
    //     }
    // }

    // Parallel segmented prefix scan (O(log m) rounds):
    // q[i] = same_attr[i] ? (q[i-1] + 1) : 0
    // Initialize:
    //   q as 64-bit values converted from same_attr bits (0/1)
    //   continue_flag = same_attr (whether this position can connect leftward)
    sbMatrix q(m, 64);
    std::memcpy(q.mShares[0].data(), same_attr.mShares[0].data(), m * sizeof(same_attr.mShares[0](0, 0)));
    std::memcpy(q.mShares[1].data(), same_attr.mShares[1].data(), m * sizeof(same_attr.mShares[1](0, 0)));
    sbMatrix continue_flag = same_attr;

    for (size_t d = 1; d < (size_t)m; d <<= 1) {
        size_t tail = (size_t)m - d;
        if (tail == 0) break;

        sbMatrix q_prev = q;
        sbMatrix flag_prev = continue_flag;

        sbMatrix q_right(tail, 64), q_left(tail, 64);
        sbMatrix flag_right(tail, 1), flag_left(tail, 1);

        std::memcpy(q_right.mShares[0].data(), q_prev.mShares[0].data() + d , tail * sizeof(q_prev.mShares[0](0, 0)));
        std::memcpy(q_right.mShares[1].data(), q_prev.mShares[1].data() + d , tail * sizeof(q_prev.mShares[1](0, 0)));
        std::memcpy(q_left.mShares[0].data(), q_prev.mShares[0].data(), tail * sizeof(q_prev.mShares[0](0, 0)));
        std::memcpy(q_left.mShares[1].data(), q_prev.mShares[1].data(), tail * sizeof(q_prev.mShares[1](0, 0)));
        std::memcpy(flag_right.mShares[0].data(), flag_prev.mShares[0].data() + d , tail * sizeof(flag_prev.mShares[0](0, 0)));
        std::memcpy(flag_right.mShares[1].data(), flag_prev.mShares[1].data() + d , tail * sizeof(flag_prev.mShares[1](0, 0)));
        std::memcpy(flag_left.mShares[0].data(), flag_prev.mShares[0].data(), tail * sizeof(flag_prev.mShares[0](0, 0)));
        std::memcpy(flag_left.mShares[1].data(), flag_prev.mShares[1].data(), tail * sizeof(flag_prev.mShares[1](0, 0)));

        sbMatrix zero_tail(tail, 64);
        bool_init_false(pIdx, zero_tail);

        sbMatrix q_left_masked(tail, 64), q_right_new(tail, 64), flag_right_new(tail, 1);

        // q_right_new = q_right + (flag_right ? q_left : 0)
        bool_cipher_selector(pIdx, flag_right, q_left, zero_tail, q_left_masked, enc, eval, runtime);
        bool_cipher_add(pIdx, q_right, q_left_masked, q_right_new, enc, eval, runtime);

        // New linkability for next doubling step.
        bool_cipher_and(pIdx, flag_right, flag_left, flag_right_new, enc, eval, runtime);

        std::memcpy(q.mShares[0].data() + d, q_right_new.mShares[0].data(), tail * sizeof(q_right_new.mShares[0](0, 0)));
        std::memcpy(q.mShares[1].data() + d, q_right_new.mShares[1].data(), tail * sizeof(q_right_new.mShares[1](0, 0)));
        std::memcpy(continue_flag.mShares[0].data() + d, flag_right_new.mShares[0].data(), tail * sizeof(flag_right_new.mShares[0](0, 0)));
        std::memcpy(continue_flag.mShares[1].data() + d, flag_right_new.mShares[1].data(), tail * sizeof(flag_right_new.mShares[1](0, 0)));
    }

    //step 1.2 bool2arith(q,alpha_2, T[0], T[1])
    si64Matrix q_si(m, 1);
    si64Matrix ii_si(m, 1);
    si64Matrix alpha1_si(m, 1), alpha_2_si(m, 1);

    sbMatrix concat_all(m, 64*3+64*max_other_num+64*key_num);
    for(size_t i=0; i< m; i++){
        concat_all.mShares[0](i, 0) = q.mShares[0](i, 0);
        concat_all.mShares[1](i, 0) = q.mShares[1](i, 0);
        concat_all.mShares[0](i, 1) = T[2].mShares[0](i, 0);
        concat_all.mShares[1](i, 1) = T[2].mShares[1](i, 0);
        concat_all.mShares[0](i, 2) = T[3].mShares[0](i, 0);
        concat_all.mShares[1](i, 2) = T[3].mShares[1](i, 0);
        for(size_t j=0; j<key_num; j++){
            concat_all.mShares[0](i, 3+j) = T[0].mShares[0](i, j);
            concat_all.mShares[1](i, 3+j) = T[0].mShares[1](i, j);
        }
        for(size_t j=0; j<max_other_num; j++){
            concat_all.mShares[0](i, 3+key_num+j) = T[1].mShares[0](i, j);
            concat_all.mShares[1](i, 3+key_num+j) = T[1].mShares[1](i, j);
        }
    }

    si64Matrix concat_all_si(m, 3+max_other_num+key_num);
    bool2arith(pIdx, concat_all, concat_all_si, enc, eval, runtime);

    q_si.mShares[0] = concat_all_si.mShares[0].block(0, 0, m, 1);
    q_si.mShares[1] = concat_all_si.mShares[1].block(0, 0, m, 1);
    alpha1_si.mShares[0] = concat_all_si.mShares[0].block(0, 1, m, 1);
    alpha1_si.mShares[1] = concat_all_si.mShares[1].block(0, 1, m, 1);
    alpha_2_si.mShares[0] = concat_all_si.mShares[0].block(0, 2, m, 1);
    alpha_2_si.mShares[1] = concat_all_si.mShares[1].block(0, 2, m, 1);
  
    // cipher_mul(pIdx, q_si, alpha_2_si, ii_si, eval, enc, runtime);
    // ii_si = ii_si + k_si;
    sbMatrix q_div_alpha1(m, 64);
    bool_cipher_div_pow2(pIdx, q, T[2], q_div_alpha1, eval, runtime);
    si64Matrix q_mod_alpha1(m, 1), q_div_alpha1_si(m, 1), q_prime(m, 1);
    bool2arith(pIdx, q_div_alpha1, q_div_alpha1_si, enc, eval, runtime);
    cipher_mul(pIdx, q_div_alpha1_si, alpha1_si, q_prime, eval, enc, runtime);
    q_mod_alpha1 = q_si - q_prime;

    si64Matrix q_mod_mul_alpha2(m, 1);
    cipher_mul(pIdx, q_mod_alpha1, alpha_2_si, q_mod_mul_alpha2, eval, enc, runtime);
    ii_si = q_div_alpha1_si + q_mod_mul_alpha2;
    // auto t2 = std::chrono::high_resolution_clock::now();
    // double step1_set_e_ii_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"align_table_step1_set_e_ii", step1_set_e_ii_time});

    //step 2:permsort(si)
    // t1 = std::chrono::high_resolution_clock::now();
    //sort by T[0] ii_si T[1]
    si64Matrix entry_key(m, 1+key_num+max_other_num);
    //T[1]
    entry_key.mShares[0].block(0, 0, m, max_other_num) = concat_all_si.mShares[0].block(0, 3+key_num, m, max_other_num);
    entry_key.mShares[1].block(0, 0, m, max_other_num) = concat_all_si.mShares[1].block(0, 3+key_num, m, max_other_num);
    //ii_si
    entry_key.mShares[0].block(0, max_other_num, m, 1) = ii_si.mShares[0].block(0, 0, m, 1);
    entry_key.mShares[1].block(0, max_other_num, m, 1) = ii_si.mShares[1].block(0, 0, m, 1);
    //T[0]
    entry_key.mShares[0].block(0, 1+max_other_num, m, key_num) = concat_all_si.mShares[0].block(0, 3, m, key_num);
    entry_key.mShares[1].block(0, 1+max_other_num, m, key_num) = concat_all_si.mShares[1].block(0, 3, m, key_num);

    si64Matrix entry_key_sorted(m, 1+key_num+max_other_num), perm(m, 1);
    si64Matrix entry_key_perm(m, key_num+1);
    entry_key_perm.mShares[0].block(0, 0, m, 1) = entry_key.mShares[0].block(0, max_other_num, m, 1);
    entry_key_perm.mShares[1].block(0, 0, m, 1) = entry_key.mShares[1].block(0, max_other_num, m, 1);
    entry_key_perm.mShares[0].block(0, 1, m, key_num) = entry_key.mShares[0].block(0, 1+max_other_num, m, key_num);
    entry_key_perm.mShares[1].block(0, 1, m, key_num) = entry_key.mShares[1].block(0, 1+max_other_num, m, key_num);
    
    genPerm(pIdx, entry_key_perm, perm, enc, eval, runtime);
    applyPerm(pIdx, perm, entry_key, entry_key_sorted, enc, eval, runtime);
    
    T_aligned.resize(2);
    T_aligned[0].resize(m, key_num);
    T_aligned[1].resize(m, max_other_num);
    T_aligned[0].mShares[0].block(0, 0, m, key_num) = entry_key_sorted.mShares[0].block(0, 1+max_other_num, m, key_num);
    T_aligned[0].mShares[1].block(0, 0, m, key_num) = entry_key_sorted.mShares[1].block(0, 1+max_other_num, m, key_num);
    T_aligned[1].mShares[0].block(0, 0, m, max_other_num) = entry_key_sorted.mShares[0].block(0, 0, m, max_other_num);
    T_aligned[1].mShares[1].block(0, 0, m, max_other_num) = entry_key_sorted.mShares[1].block(0, 0, m, max_other_num);
    // t2 = std::chrono::high_resolution_clock::now();
    // double step2_sort_j_ii_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"align_table_step2_sort_j_ii", step2_sort_j_ii_time});
    // auto total_end = std::chrono::high_resolution_clock::now();
    // double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    // timing_results.push_back({"align_table_total", total_time});
    // if (pIdx == 0) {
    //     std::string filename = "./join_timing_results_role0.txt";
    //     std::ofstream outFile(filename, std::ios::app);
    //     if (outFile.is_open()) {
    //         outFile << "--- align_table Internal Timing (Role " << pIdx << ") ---" << std::endl;
    //         outFile << std::fixed << std::setprecision(3);
    //         for(const auto& result : timing_results) {
    //             outFile << "  " << result.first << ": " << result.second << " ms" << std::endl;
    //         }
    //         outFile << "----------------------------------------" << std::endl;
    //         outFile << std::endl;
    //         outFile.close();
    //     }
    // }

        return ;
    }
/*
void align_table(int pIdx, std::vector<sbMatrix> &T, std::vector<sbMatrix> &T_aligned,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    // 时间测量变量
    auto total_start = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<std::string, double>> timing_results;
    
    int m=T[0].rows();
    int key_num = T[0].i64Cols();
    int max_other_num = T[1].i64Cols();


    //sbMatrix q(m, 64), k(m, 64);
    sbMatrix q_k_concat(2*m, 64);
    bool_init_false(pIdx, q_k_concat);
    //bool_init_false(pIdx, k);

    sbMatrix false_matrix(1,1);
    bool_init_false(pIdx, false_matrix);
    sbMatrix same_attr_partial(m-1, 1);
    compare_consecutive_rows_bool(pIdx, T[0], same_attr_partial, enc, eval, runtime);
    sbMatrix same_attr(m, 1);
    same_attr.mShares[0](0, 0) = false_matrix.mShares[0](0, 0);
    same_attr.mShares[1](0, 0) = false_matrix.mShares[1](0, 0);
    std::memcpy(same_attr.mShares[0].data() + 1, same_attr_partial.mShares[0].data(), (m-1) * sizeof(same_attr_partial.mShares[0](0, 0)));
    std::memcpy(same_attr.mShares[1].data() + 1, same_attr_partial.mShares[1].data(), (m-1) * sizeof(same_attr_partial.mShares[1](0, 0)));
    


    //step 1: set e.ii
    auto t1 = std::chrono::high_resolution_clock::now();
    for(size_t i=0; i<m; i++){
        //j是否相同: same_attr
        // sbMatrix j_i(1, 64*key_num);
        // for(size_t l=0; l<key_num; l++){
        //     j_i.mShares[0](0, l) = T[0].mShares[0](i, l);
        //     j_i.mShares[1](0, l) = T[0].mShares[1](i, l);
        // }

        // sbMatrix same_attr(1, 1), not_same_attr(1,1);
        // if(i == 0){
        //     bool_init_false(pIdx, same_attr);
        // } else {
        //     sbMatrix j_i_prev(1, 64*key_num);
        //     for(size_t l=0; l<key_num; l++){
        //         j_i_prev.mShares[0](0, l) = T[0].mShares[0](i-1, l);
        //         j_i_prev.mShares[1](0, l) = T[0].mShares[1](i-1, l);
        //     }
        //     bool_cipher_eq(pIdx, j_i, j_i_prev, same_attr, enc, eval, runtime);
        // }

        sbMatrix same_attr_i(2, 1);
        same_attr_i.mShares[0](0, 0) = same_attr.mShares[0](i, 0);
        same_attr_i.mShares[1](0, 0) = same_attr.mShares[1](i, 0);
        same_attr_i.mShares[0](1, 0) = same_attr.mShares[0](i, 0);
        same_attr_i.mShares[1](1, 0) = same_attr.mShares[1](i, 0);

        sbMatrix zero(2, 64);
        bool_init_false(pIdx, zero);
        sbMatrix one(2, 64);
        bool_init_true(pIdx, one);

        //sbMatrix current_q(1, 64),current_k(1, 64);
        sbMatrix current_q_k(2, 64);
        if(i == 0) {
            current_q_k = zero;
            // current_q = zero;
            // current_k = zero;
        } else {
            current_q_k.mShares[0](0, 0) = q_k_concat.mShares[0](i-1, 0);
            current_q_k.mShares[1](0, 0) = q_k_concat.mShares[1](i-1, 0);
            current_q_k.mShares[0](1, 0) = q_k_concat.mShares[0](i+m-1, 0);
            current_q_k.mShares[1](1, 0) = q_k_concat.mShares[1](i+m-1, 0);
            // current_q.mShares[0](0, 0) = q.mShares[0](i-1, 0);
            // current_q.mShares[1](0, 0) = q.mShares[1](i-1, 0);
            // current_k.mShares[0](0, 0) = k.mShares[0](i-1, 0);
            // current_k.mShares[1](0, 0) = k.mShares[1](i-1, 0);
        }

        sbMatrix q_k_plus(2, 64);
        bool_cipher_add(pIdx, current_q_k, one, q_k_plus, enc, eval, runtime);
        sbMatrix true_value_1(2, 64);
        true_value_1.mShares[0](0, 0) = q_k_plus.mShares[0](0, 0);
        true_value_1.mShares[1](0, 0) = q_k_plus.mShares[1](0, 0);
        true_value_1.mShares[0](1, 0) = current_q_k.mShares[0](1, 0);
        true_value_1.mShares[1](1, 0) = current_q_k.mShares[1](1, 0);
        // sbMatrix q_plus(1, 64);
        // bool_cipher_add(pIdx, current_q, one, q_plus, enc, eval, runtime);
        // sbMatrix k_plus(1, 64);
        // bool_cipher_add(pIdx, current_k, one, k_plus, enc, eval, runtime);

        //if same :q++, k不变 else :q=0,k=0
        //sbMatrix q_mid(1, 64), k_mid(1, 64);
        sbMatrix q_k_mid(2, 64);
        bool_cipher_selector(pIdx, same_attr_i, true_value_1, zero, q_k_mid, enc, eval, runtime);
        // bool_cipher_selector(pIdx, same_attr, q_plus, zero, q_mid, enc, eval, runtime);
        // bool_cipher_selector(pIdx, same_attr, current_k, zero, k_mid, enc, eval, runtime);
   
        //if  q_new>alpha_1-1 : q=0 ,k++ ; else q=0,k不变
        sbMatrix cond(1, 1);
        sbMatrix q_mid(1, 64);
        q_mid.mShares[0](0, 0) = q_k_mid.mShares[0](0, 0);
        q_mid.mShares[1](0, 0) = q_k_mid.mShares[1](0, 0);

        sbMatrix alpha_1_i(1, 64),alpha_1_minus_one(1, 64),one_one(1, 64);
        bool_init_true(pIdx, one_one);
        alpha_1_i.mShares[0](0, 0) = T[2].mShares[0](i, 0);
        alpha_1_i.mShares[1](0, 0) = T[2].mShares[1](i, 0);
        bool_cipher_sub(pIdx, alpha_1_i, one_one, alpha_1_minus_one, enc, eval, runtime);
        bool_cipher_lt(pIdx, alpha_1_minus_one, q_mid, cond, enc, eval, runtime);

        sbMatrix cond_concat(2, 1);
        cond_concat.mShares[0](0, 0) = cond.mShares[0](0, 0);
        cond_concat.mShares[1](0, 0) = cond.mShares[1](0, 0);
        cond_concat.mShares[0](1, 0) = cond.mShares[0](0, 0);
        cond_concat.mShares[1](1, 0) = cond.mShares[1](0, 0);
        sbMatrix true_value_2(2, 64);
        true_value_2.mShares[0](0, 0) = zero.mShares[0](0, 0);
        true_value_2.mShares[1](0, 0) = zero.mShares[1](0, 0);
        true_value_2.mShares[0](1, 0) = q_k_plus.mShares[0](1, 0);
        true_value_2.mShares[1](1, 0) = q_k_plus.mShares[1](1, 0);
        sbMatrix q_k_new(2, 64);
        bool_cipher_selector(pIdx, cond_concat, true_value_2, q_k_mid, q_k_new, enc, eval, runtime);
        // sbMatrix q_new(1, 64), k_new(1, 64);
        // bool_cipher_selector(pIdx, cond, zero, q_mid, q_new, enc, eval, runtime);
        // bool_cipher_selector(pIdx, cond, k_plus, k_mid, k_new, enc, eval, runtime);


        q_k_concat.mShares[0](i, 0) = q_k_new.mShares[0](0, 0);
        q_k_concat.mShares[1](i, 0) = q_k_new.mShares[1](0, 0);
        q_k_concat.mShares[0](i+m, 0) = q_k_new.mShares[0](1, 0);
        q_k_concat.mShares[1](i+m, 0) = q_k_new.mShares[1](1, 0);

    }

    //这个好像必须得转成si64才能算算术乘法
    si64Matrix q_si(m, 1), k_si(m, 1);
    si64Matrix ii_si(m, 1);
    si64Matrix alpha_2_si(m, 1);
    sbMatrix q_k_alpha_concat(3*m, 64);
    std::memcpy(q_k_alpha_concat.mShares[0].data(), q_k_concat.mShares[0].data(), 2*m * sizeof(q_k_concat.mShares[0](0, 0)));
    std::memcpy(q_k_alpha_concat.mShares[1].data(), q_k_concat.mShares[1].data(), 2*m * sizeof(q_k_concat.mShares[1](0, 0)));
    std::memcpy(q_k_alpha_concat.mShares[0].data() + 2*m, T[3].mShares[0].data(), m * sizeof(T[3].mShares[0](0, 0)));
    std::memcpy(q_k_alpha_concat.mShares[1].data() + 2*m, T[3].mShares[1].data(), m * sizeof(T[3].mShares[1](0, 0)));
    
    si64Matrix q_k_alpha_si(3*m, 1);
    bool2arith(pIdx, q_k_alpha_concat, q_k_alpha_si, enc, eval, runtime);
    std::memcpy(q_si.mShares[0].data(), q_k_alpha_si.mShares[0].data(), m * sizeof(q_k_alpha_si.mShares[0](0, 0)));
    std::memcpy(q_si.mShares[1].data(), q_k_alpha_si.mShares[1].data(), m * sizeof(q_k_alpha_si.mShares[1](0, 0)));
    std::memcpy(k_si.mShares[0].data(), q_k_alpha_si.mShares[0].data()+ m, m * sizeof(q_k_alpha_si.mShares[0](0, 0)));
    std::memcpy(k_si.mShares[1].data(), q_k_alpha_si.mShares[1].data()+ m, m * sizeof(q_k_alpha_si.mShares[1](0, 0)));
    std::memcpy(alpha_2_si.mShares[0].data(), q_k_alpha_si.mShares[0].data() + 2*m, m * sizeof(q_k_alpha_si.mShares[0](0, 0)));
    std::memcpy(alpha_2_si.mShares[1].data(), q_k_alpha_si.mShares[1].data() + 2*m, m * sizeof(q_k_alpha_si.mShares[1](0, 0)));

    // bool2arith(pIdx, q, q_si, enc, eval, runtime);
    // bool2arith(pIdx, k, k_si, enc, eval, runtime);
    // bool2arith(pIdx, T[3], alpha_2_si, enc, eval, runtime);
    cipher_mul(pIdx, q_si, alpha_2_si, ii_si, eval, enc, runtime);
    ii_si = ii_si + k_si;
    auto t2 = std::chrono::high_resolution_clock::now();
    double step1_set_e_ii_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"align_table_step1_set_e_ii", step1_set_e_ii_time});

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
    t1 = std::chrono::high_resolution_clock::now();
    sbMatrix  ii_sb(m, 64);
    arith2bool(pIdx, ii_si, ii_sb, enc, eval, runtime);


    sbMatrix entry_key(m, 64*key_num+64*max_other_num+64);
    //合并顺序为：T[1],ii_sb,T[0]
    for(size_t i=0; i<m; i++){
        //T[0]
        for(int l=key_num + max_other_num; l>=max_other_num+1 ; l--){
            entry_key.mShares[0](i, l) = T[0].mShares[0](i, key_num + max_other_num -l);
            entry_key.mShares[1](i, l) = T[0].mShares[1](i, key_num + max_other_num -l);
        }
        //ii_sb
        entry_key.mShares[0](i, max_other_num) = ii_sb.mShares[0](i, 0);
        entry_key.mShares[1](i, max_other_num) = ii_sb.mShares[1](i, 0);
        //T[1]
        for(int l=max_other_num-1; l>=0; l--){
            entry_key.mShares[0](i, l) = T[1].mShares[0](i, max_other_num-1-l);
            entry_key.mShares[1](i, l) = T[1].mShares[1](i, max_other_num-1-l);
        }
    }

    sbMatrix entry_key_sorted(m, 64*key_num+64*max_other_num+64);
    odd_even_merge_sort(entry_key, entry_key_sorted, pIdx, enc, eval, runtime);

    T_aligned.resize(2);
    T_aligned[0].resize(m, 64*key_num);
    T_aligned[1].resize(m, 64*max_other_num);

    for(size_t i=0; i< m; i++){
        //T[0]
        for(int l=key_num + max_other_num; l>=max_other_num+1 ; l--){
            T_aligned[0].mShares[0](i, key_num + max_other_num -l) = entry_key_sorted.mShares[0](i, l);
            T_aligned[0].mShares[1](i, key_num + max_other_num -l) = entry_key_sorted.mShares[1](i, l);
        }
        //T[1]
        for(int l=max_other_num-1; l>=0; l--){
            T_aligned[1].mShares[0](i, max_other_num-1-l) = entry_key_sorted.mShares[0](i, l);
            T_aligned[1].mShares[1](i, max_other_num-1-l) = entry_key_sorted.mShares[1](i, l);
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    double step2_sort_j_ii_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    timing_results.push_back({"align_table_step2_sort_j_ii", step2_sort_j_ii_time});

    // 计算总时间
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    timing_results.push_back({"align_table_total", total_time});

    // 将时间测量结果写入文件（仅 role 0 写入）
    if (pIdx == 0) {
        std::string filename = "./join_timing_results_role0.txt";
        std::ofstream outFile(filename, std::ios::app);
        if (outFile.is_open()) {
            outFile << "--- align_table Internal Timing (Role " << pIdx << ") ---" << std::endl;
            outFile << std::fixed << std::setprecision(3);
            for(const auto& result : timing_results) {
                outFile << "  " << result.first << ": " << result.second << " ms" << std::endl;
            }
            outFile << "----------------------------------------" << std::endl;
            outFile << std::endl;
            outFile.close();
        }
    }

    return ;
    
}
*/

void join(int pIdx, std::vector<si64Matrix> &T_1_key, std::vector<si64Matrix> &T_1_other,
     std::vector<si64Matrix> &T_2_key, std::vector<si64Matrix> &T_2_other,
     std::vector<si64Matrix> &T_joined,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    // // 时间测量变量（按需保留，当前禁用）
    // auto total_start = std::chrono::high_resolution_clock::now();
    // std::vector<std::pair<std::string, double>> timing_results;
    
    int key_num = T_1_key.size();
    int other_1_num = T_1_other.size();
    int other_2_num = T_2_other.size();
    int max_other_num = std::max(other_1_num, other_2_num);
    
    // // 测量 augment_table 时间
    // auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<si64Matrix> T_1_auged(4),T_2_auged(4);
    augment_table(pIdx, T_1_key, T_1_other, T_2_key, T_2_other, T_1_auged, T_2_auged, enc, eval, runtime);
    // auto t2 = std::chrono::high_resolution_clock::now();
    // double augment_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"augment_table", augment_time});
 
    // // 测量 oblivious_expand (T_1) 时间
    // t1 = std::chrono::high_resolution_clock::now();
    std::vector<sbMatrix> T_1_expanded(4);
    i64Matrix output_size(1, 1);
    oblivious_expand(pIdx, T_1_auged,  T_1_expanded, 0, output_size, enc, eval, runtime);
    // t2 = std::chrono::high_resolution_clock::now();
    // double expand_t1_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"oblivious_expand_T1", expand_t1_time});
    
    // // 测量 oblivious_expand (T_2) 时间
    // t1 = std::chrono::high_resolution_clock::now();
    std::vector<sbMatrix> T_2_expanded(4);
    oblivious_expand(pIdx, T_2_auged,  T_2_expanded, 1, output_size, enc, eval, runtime);
    // t2 = std::chrono::high_resolution_clock::now();
    // double expand_t2_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"oblivious_expand_T2", expand_t2_time});
    
    i64 m = output_size(0, 0);

    // // 测量 align_table 时间
    // t1 = std::chrono::high_resolution_clock::now();
    std::vector<si64Matrix> T_2_aligned(2);
    align_table(pIdx,T_2_expanded, T_2_aligned, enc, eval, runtime);
    // t2 = std::chrono::high_resolution_clock::now();
    // double align_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"align_table", align_time});

    // // 测量 concat join 时间
    // t1 = std::chrono::high_resolution_clock::now();
    T_joined.resize(key_num+other_1_num+other_2_num);
    for(size_t i=0; i<key_num+other_1_num+other_2_num; i++){
        T_joined[i].resize(m, 1);
    }

    sbMatrix T1_joined_sb_concat(m, 64*(key_num+other_1_num));
    si64Matrix T1_joined_si_concat(m, key_num+other_1_num);

    for(size_t i=0; i<m; i++){
        for(size_t l=0; l<key_num; l++){
            T1_joined_sb_concat.mShares[0](i, l) = T_1_expanded[0].mShares[0](i, l);
            T1_joined_sb_concat.mShares[1](i, l) = T_1_expanded[0].mShares[1](i, l);
        }
        //if(pIdx == 0){std::cout<<"0.1"<<std::endl;}
        for(size_t l=0; l<other_1_num; l++){
            T1_joined_sb_concat.mShares[0](i, l+key_num) = T_1_expanded[1].mShares[0](i, l);
            T1_joined_sb_concat.mShares[1](i, l+key_num) = T_1_expanded[1].mShares[1](i, l);
        }
        
    }
    
    bool2arith(pIdx, T1_joined_sb_concat, T1_joined_si_concat, enc, eval, runtime);


    for(size_t i=0; i<key_num; i++){
        T_joined[i].mShares[0].block(0, 0, m, 1) = T1_joined_si_concat.mShares[0].block(0, key_num-i-1, m, 1);
        T_joined[i].mShares[1].block(0, 0, m, 1) = T1_joined_si_concat.mShares[1].block(0, key_num-i-1, m, 1);
    }
    for(size_t i=0; i<other_1_num; i++){
        T_joined[i+key_num].mShares[0].block(0, 0, m, 1) = T1_joined_si_concat.mShares[0].block(0, key_num+other_1_num-i-1, m, 1);
        T_joined[i+key_num].mShares[1].block(0, 0, m, 1) = T1_joined_si_concat.mShares[1].block(0, key_num+other_1_num-i-1, m, 1);
    }
    for(size_t i=0; i<other_2_num; i++){
        T_joined[i+key_num+other_1_num].mShares[0].block(0, 0, m, 1) = T_2_aligned[1].mShares[0].block(0, max_other_num-i-1, m, 1);
        T_joined[i+key_num+other_1_num].mShares[1].block(0, 0, m, 1) = T_2_aligned[1].mShares[1].block(0, max_other_num-i-1, m, 1);
    }
    // t2 = std::chrono::high_resolution_clock::now();
    // double concat_join_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    // timing_results.push_back({"concat_join", concat_join_time});
    // auto total_end = std::chrono::high_resolution_clock::now();
    // double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    // timing_results.push_back({"total", total_time});
    // if (pIdx == 0) {
    //     std::string filename = "./join_timing_results_role0.txt";
    //     std::ofstream outFile(filename, std::ios::app);
    //     if (outFile.is_open()) {
    //         outFile << "========== Join Function Timing (Role " << pIdx << ") ==========" << std::endl;
    //         outFile << std::fixed << std::setprecision(3);
    //         for(const auto& result : timing_results) {
    //             outFile << result.first << ": " << result.second << " ms" << std::endl;
    //         }
    //         outFile << "================================================" << std::endl;
    //         outFile << std::endl;
    //         outFile.close();
    //     }
    // }

    return ;
}


void filter(int pIdx, std::vector<si64Matrix> &T, int filColIdx, int value, bool is_scalar, std::string op_str,
    std::vector<si64Matrix> &T_filtered,
    Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    int n = T[0].rows();
    int m = T.size();
    
    si64Matrix filCol = T[filColIdx];
    si64Matrix valCol(n, 1);
    if(is_scalar){
        set_const_share(pIdx, value, valCol, enc, eval, runtime);
    }
    else{
        valCol = T[value];
    }
    

    sbMatrix flag(n,1);
    if(op_str == "="){
        circuit_cipher_eq(pIdx, filCol, valCol, flag, eval, runtime);
    }
    else if(op_str == ">"){
        cipher_gt(pIdx, filCol, valCol, flag, eval, runtime);
    }
    else if(op_str == ">="){
        cipher_ge(pIdx, filCol, valCol, flag, eval, enc, runtime);
    }
    else if(op_str == "<"){
        cipher_gt(pIdx, valCol, filCol, flag, eval, runtime);
    }
    else if(op_str == "<="){
        cipher_ge(pIdx, valCol, filCol, flag, eval, enc, runtime);
    }


    std::vector<si64Matrix> res(m);
    si64Matrix T_vector(m*n, 1), res_vector(m*n, 1);
    sbMatrix flag_vector(m*n, 1);
    for(size_t i=0; i<m; i++){
        // for(size_t j=0; j<n; j++){
        //     T_vector.mShares[0](i*n+j, 0) = T[i].mShares[0](j, 0);
        //     T_vector.mShares[1](i*n+j, 0) = T[i].mShares[1](j, 0);
        // }
        std::memcpy(T_vector.mShares[0].data() + i*n, T[i].mShares[0].data(), n * sizeof(T[i].mShares[0](0, 0)));
        std::memcpy(T_vector.mShares[1].data() + i*n, T[i].mShares[1].data(), n * sizeof(T[i].mShares[1](0, 0)));
        std::memcpy(flag_vector.mShares[0].data() + i*n, flag.mShares[0].data(), n * sizeof(flag.mShares[0](0, 0)));
        std::memcpy(flag_vector.mShares[1].data() + i*n, flag.mShares[1].data(), n * sizeof(flag.mShares[1](0, 0)));
    }
    cipher_mul(pIdx, T_vector, flag_vector, res_vector, eval, enc, runtime);
    for(size_t i=0; i<m; i++){
        res[i].resize(n, 1);
        //cipher_mul(pIdx, T[i], flag, res[i], eval, enc, runtime);
        std::memcpy(res[i].mShares[0].data(), res_vector.mShares[0].data() + i*n, n * sizeof(res_vector.mShares[0](0, 0)));
        std::memcpy(res[i].mShares[1].data(), res_vector.mShares[1].data() + i*n, n * sizeof(res_vector.mShares[1](0, 0)));
    }

    T_filtered.resize(m);
    for(size_t i=0; i<m; i++){
        T_filtered[i].resize(n, 1);
    }


    shuffle(pIdx, res, T_filtered, enc, eval, runtime);

    return ;
}

void semi_join(int pIdx, std::vector<si64Matrix> &T_1_key, std::vector<si64Matrix> &T_1_other,
    std::vector<si64Matrix> &T_2_key, std::vector<si64Matrix> &T_2_other,
    std::vector<si64Matrix> &T_semi_joined,
   Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

    int n = T_1_key[0].rows();
    int key_num = T_1_key.size();
    int other_1_num = T_1_other.size();
    int total_cols = key_num + other_1_num;

    std::vector<si64Matrix> T_1_auged(4),T_2_auged(4);
    augment_table(pIdx, T_1_key, T_1_other, T_2_key, T_2_other, T_1_auged, T_2_auged, enc, eval, runtime);

    //T_1_auged: [0]=key(n,key_num), [1]=other(n,max_other_num), [2]=alpha1, [3]=alpha2
    // alpha2 != 0 表示该行的key在T_2中存在匹配
    si64Matrix alpha2(n, 1), zero(n, 1);
    alpha2 = T_1_auged[3];
    set_const_share(pIdx, 0, zero, enc, eval, runtime);

    // eq_flag=1 当 alpha2==0（无匹配），取反得到 flag=1 当有匹配
    sbMatrix eq_flag(n, 1), flag(n, 1);
    cipher_eq(pIdx, alpha2, zero, eq_flag, eval, runtime);
    bool_cipher_not(pIdx, eq_flag, flag);
    for(size_t i = 0; i < (size_t)n; i++){
        flag.mShares[0](i, 0) &= 1;
        flag.mShares[1](i, 0) &= 1;
    }


    // 用 flag 乘以各列，无匹配的行归零
    // 将所有列 pack 成一个大向量，一次 cipher_mul 完成
    si64Matrix T_1_vector(n * total_cols, 1);
    sbMatrix flag_vector(n * total_cols, 1);

    // T_1_auged[0] 是 (n, key_num) 的列主序矩阵，逐列 memcpy
    for(int i = 0; i < key_num; i++){
        T_1_vector.mShares[0].block(i * n, 0, n, 1) = T_1_auged[0].mShares[0].col(i);
        T_1_vector.mShares[1].block(i * n, 0, n, 1) = T_1_auged[0].mShares[1].col(i);
    }
    for(int i = 0; i < other_1_num; i++){
        T_1_vector.mShares[0].block((key_num + i) * n, 0, n, 1) = T_1_auged[1].mShares[0].col(i);
        T_1_vector.mShares[1].block((key_num + i) * n, 0, n, 1) = T_1_auged[1].mShares[1].col(i);
    }
    for(int i = 0; i < total_cols; i++){
        std::memcpy(flag_vector.mShares[0].data() + i * n, flag.mShares[0].data(), n * sizeof(i64));
        std::memcpy(flag_vector.mShares[1].data() + i * n, flag.mShares[1].data(), n * sizeof(i64));
    }

    //
    si64Matrix res_vector(n * total_cols, 1);
    cipher_mul(pIdx, T_1_vector, flag_vector, res_vector, eval, enc, runtime);

    // 拆分结果，注意 auged 中列序是反的：key[0]在最高位列(key_num-1)，需要反转回来
    T_semi_joined.resize(total_cols);
    for(int i = 0; i < key_num; i++){
        T_semi_joined[i].resize(n, 1);
        std::memcpy(T_semi_joined[i].mShares[0].data(), res_vector.mShares[0].data() + (key_num - i - 1) * n, n * sizeof(i64));
        std::memcpy(T_semi_joined[i].mShares[1].data(), res_vector.mShares[1].data() + (key_num - i - 1) * n, n * sizeof(i64));
    }
    for(int i = 0; i < other_1_num; i++){
        T_semi_joined[i + key_num].resize(n, 1);
        std::memcpy(T_semi_joined[i + key_num].mShares[0].data(), res_vector.mShares[0].data() + (key_num + other_1_num - i - 1) * n, n * sizeof(i64));
        std::memcpy(T_semi_joined[i + key_num].mShares[1].data(), res_vector.mShares[1].data() + (key_num + other_1_num - i - 1) * n, n * sizeof(i64));
    }

   }



