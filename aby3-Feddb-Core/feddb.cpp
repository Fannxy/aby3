#include "feddb.h"
#include "utils.h"
#include "../aby3-GORAM-Core/Shuffle.h"
#include "../aby3-GORAM-Core/Sort.h"
#include "../aby3-GORAM-Core/Basics.h"

using namespace oc;
using namespace aby3;

//sbMatrix：一列
// void concatRows_sbMatrix(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB,
//     sbMatrix &res, Sh3Encryptor &enc, Sh3Evaluator &eval,
//     Sh3Runtime &runtime){
//         assert(sharedA.bitCount() == sharedB.bitCount());
//         //将sharedA和sharedB按行拼接成res
//         int rowA = sharedA.rows();
//         int rowB = sharedB.rows();
//         int rowRes = rowA + rowB;
//         res.resize(rowRes, sharedA.bitCount());
//         res.mShares[0].block(0, 0, rowA, sharedA.i64Cols()) = sharedA.mShares[0];
//         res.mShares[1].block(0, 0, rowA, sharedA.i64Cols()) = sharedA.mShares[1];
//         res.mShares[0].block(rowA, 0, rowB, sharedB.i64Cols()) = sharedB.mShares[0];
//         res.mShares[1].block(rowA, 0, rowB, sharedB.i64Cols()) = sharedB.mShares[1];
//         return;
//     }

// void concatRows_i64Matrix(int pIdx, i64Matrix &sharedA, i64Matrix &sharedB,
//     i64Matrix &res, Sh3Encryptor &enc, Sh3Evaluator &eval,
//     Sh3Runtime &runtime){
//         assert(sharedA.bitCount() == sharedB.bitCount());
//         //将sharedA和sharedB按行拼接成res
//         int rowA = sharedA.rows();
//         int rowB = sharedB.rows();
//         int rowRes = rowA + rowB;
//         res.resize(rowRes, sharedA.bitCount());
//         res.block(0, 0, rowA, sharedA.i64Cols()) = sharedA;
//         res.block(rowA, 0, rowB, sharedB.i64Cols()) = sharedB;
//         return;
//     }

// void shuffle_sbMatrix(int pIdx, sbMatrix &T, sbMatrix &Tres, 
//     Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

//         efficient_shuffle_full_matrix(T, pIdx, Tres, enc, eval, runtime);
//         return ;
//     }

// void project_sbMatrix(int pIdx, sbMatrix &T, std::vector<int> &cols, sbMatrix &Tres, 
//     Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
        
//     size_t rows = T.rows();
//     size_t selected_cols = cols.size();
    
//     Tres.resize(rows, T.bitCount());
    
//     // Copy the selected columns from input to output
//     for (size_t i = 0; i < rows; i++) {
//         for (size_t j = 0; j < selected_cols; j++) {
//             // Copy both shares for the selected column
//             Tres.mShares[0](i, j) = T.mShares[0](i, cols[j]);
//             Tres.mShares[1](i, j) = T.mShares[1](i, cols[j]);
//         }
//     }
        
//     return ;
// }

// void project_i64Matrix(int pIdx, i64Matrix &T, std::vector<int> &cols, i64Matrix &Tres, 
//     Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
        
//     size_t rows = T.rows();
//     size_t selected_cols = cols.size();
    
//     Tres.resize(rows, T.bitCount());

//     for (size_t i = 0; i < rows; i++) {
//         for (size_t j = 0; j < selected_cols; j++) {
//             Tres(i, j) = T(i, cols[j]);
//         }
//     }
        
//     return ;
// }

// 原始左表；标志表；行数信息
// void flag_join_sbMatrix(int pIdx, sbMatrix &left_table, sbMatrix &flags_table, si64Matrix &num_rows_info,sbMatrix &result,              
//         Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime) {
    
//     size_t left_rows = left_table.rows();
//     size_t left_cols = left_table.cols();
//     size_t left_bits = left_table.bitCount();
    
//     size_t flags_rows = flags_table.rows();
    
//     // 恢复期望行数的明文值
//     i64Matrix expected_result_rows_plaintext(1, 1);
//     enc.revealAll(runtime, num_rows_info, expected_result_rows_plaintext).get();
//     size_t expected_result_rows = static_cast<size_t>(expected_result_rows_plaintext(0, 0));
    
//     // 验证输入维度
//     if (left_rows != flags_rows) {
//         throw std::runtime_error("FlagJoin: left table rows must match flags table rows");
//     }
    
//     // 初始化结果矩阵为期望的大小
//     result.resize(expected_result_rows, left_bits);
    
//     size_t write_pos = 0;
//     for (size_t i = 0; i < left_rows && write_pos < expected_result_rows; i++) {
//         // 获取每行的标志
//         sbMatrix flag_val(1, 1);
//         flag_val.mShares[0](0, 0) = flags_table.mShares[0](i, 0);
//         flag_val.mShares[1](0, 0) = flags_table.mShares[1](i, 0);
        
//         // 检查标志是否为零
//         i64Matrix plain_zero(1, 1);
//         plain_zero(0, 0) = 0;
//         sbMatrix is_flag_zero(1, 1);
//         bool_cipher_eq(pIdx, flag_val, plain_zero, is_flag_zero, enc, eval, runtime);
//         i64Matrix plain_is_zero(1, 1);
//         enc.revealAll(runtime, is_flag_zero, plain_is_zero).get();
        
//         // 如果标志不为零，则执行条件复制并写入结果
//         if (plain_is_zero(0, 0) == 0) {
//             // 对于左表的每一列，执行条件复制
//             for (size_t j = 0; j < left_cols; j++) {
//                 // 获取左表的元素
//                 sbMatrix left_element(1, left_bits);
//                 left_element.mShares[0](0, 0) = left_table.mShares[0](i, j);
//                 left_element.mShares[1](0, 0) = left_table.mShares[1](i, j);
                
//                 sbMatrix masked_element(1, left_bits);
//                 bool_cipher_and(pIdx, flag_val, left_element, masked_element, enc, eval, runtime);
                
//                 result.mShares[0](write_pos, j) = masked_element.mShares[0](0, 0);
//                 result.mShares[1](write_pos, j) = masked_element.mShares[1](0, 0);
//             }
//             write_pos++;
//         }
//     }
    
//     return;
// }

//TODO：join

// template<typename MatrixType>
// void persist_cipher(int pIdx, string &table_name, MatrixType &T){
//     string filename = "GORAM-ABY3/aby3/aby3-Feddb-tmpfile/cipher" + table_name + "_" + to_string(pIdx) + ".txt";
    
//     ofstream outFile(filename);
//     if (!outFile.is_open()) {
//         cerr << "Error: Unable to open file " << filename << " for writing" << endl;
//         return;
//     }
    
//     int rows = T.rows();
//     int cols = T.cols();
    
//     outFile << "Matrix dimensions: " << rows << " x " << cols << endl;

//     outFile << "Matrix shares[0]:" << endl;
//     //T.mShares[0]和T.mShares[1]分别写入
//     for (int i = 0; i < rows; i++) {
//         for (int j = 0; j < cols; j++) {
//             outFile << T.mShares[0](i, j) << " ";
//         }
//         outFile << endl;
//     }
//     outFile << "Matrix shares[1]:" << endl;
//     for (int i = 0; i < rows; i++) {
//         for (int j = 0; j < cols; j++) {
//             outFile << T.mShares[1](i, j) << " ";
//         }
//         outFile << endl;
//     }

//     outFile.close();

//     return;
// }

// void persist_plain(int pIdx, string &table_name, i64Matrix &T){
//     string filename = "GORAM-ABY3/aby3/aby3-Feddb-tmpfile/i64Matrix" + table_name + "_" + to_string(pIdx) + ".txt";
    
//     ofstream outFile(filename);
//     if (!outFile.is_open()) {
//         cerr << "Error: Unable to open file " << filename << " for writing" << endl;
//         return;
//     }

//     int rows = T.rows();
//     int cols = T.cols();

//     outFile << "Matrix dimensions: " << rows << " x " << cols << endl;

//     outFile << "Matrix values:" << endl;
//     for (int i = 0; i < rows; i++) {
//         for (int j = 0; j < cols; j++) {
//             outFile << T(i, j) << " ";
//         }
//         outFile << endl;
//     }

//     outFile.close();

//     return;
// }



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





