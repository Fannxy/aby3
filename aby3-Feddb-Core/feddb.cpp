#include "feddb.h"
#include "utils.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3-GORAM-Core/Sort.h"
#include "../aby3-GORAM-Core/Basics.h"

using namespace oc;
using namespace aby3;

// void shuffle(int pIdx, std::vector<si64Matrix>& T, std::vector<si64Matrix> &Tres, 
//     Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
//         size_t len = T.size();
//         std::vector<sbMatrix> T_sb(len),Tres_sb(len);
//         for(size_t i=0; i<len; i++){
//             //不行，arith2bool转成的bitcount是64bit的，但是efficient_shuffle只支持1bit的
//            // T_sb[i].resize(T[i].rows(), 64);
//             arith2bool(pIdx, T[i], T_sb[i], enc, eval, runtime);
//             //std::cout<<"T_sb[i].bitCount()"<<T_sb[i].bitCount()<<std::endl;
//         }

//         efficient_shuffle(T_sb, pIdx, Tres_sb, enc, eval, runtime);
        
//         for(size_t i=0; i<len; i++){
//             Tres[i].resize(T[i].rows(), 1);
//             bool2arith(pIdx, Tres_sb[i], Tres[i], enc, eval, runtime);
//         }
//         return ;
//     }

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
    std::string filename = "/root/GORAM-ABY3/aby3/aby3-Feddb-tmpfile/plain/" + table_name + "_" + std::to_string(pIdx) + ".txt";
    
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing" << std::endl;
        return;
    }

    int rows = T[0].rows();
    outFile << "Matrix dimensions: " << rows << std::endl;

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
    std::string filename = "/root/GORAM-ABY3/aby3/aby3-Feddb-tmpfile/plain/" + table_name + "_" + std::to_string(pIdx) + ".txt";
    
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for reading" << std::endl;
        return;
    }

    std::string dummy;
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






