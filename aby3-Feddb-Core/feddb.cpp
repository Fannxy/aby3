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

    //DEBUG
    i64Matrix t_plain(t_len, idx.cols());
    enc.revealAll(runtime, t, t_plain).get();
    std::cout << "t_plain: " << std::endl;
    for(size_t i=0; i<t_len; i++){  
        std::cout << t_plain(i, 0) << " ";
    }
    std::cout << std::endl;  
    std::cout.flush();       
    //----t_plain correct

    //step one: 对t进行argsort得到sigma
    si64Matrix sigma_si(t_len, 1);

    //TODO cipher_argsort or rtr_cipher_argsort
    fed_argsort(pIdx,t,sigma_si, enc, eval, runtime);
    //DEBUG
    i64Matrix temp(t_len, 1);
    enc.revealAll(runtime, sigma_si, temp).get();
    std::cout << "sigma_plain: " << std::endl;
    for(size_t i=0; i<t_len; i++){  
        std::cout << temp(i, 0) << " ";
    }
    std::cout << std::endl;  
    std::cout.flush();      

    //step two: prefixsum_{-1}(v) && u_si
    si64Matrix prefix_inv(v_len, v.cols());
    prefixsum_inv(pIdx, v, prefix_inv);
    si64Matrix u_si(t_len, v.cols());
    for(size_t i=0; i<v_len; i++){
        u_si.mShares[0](i, 0) = prefix_inv.mShares[0](i, 0);
        u_si.mShares[1](i, 0) = prefix_inv.mShares[1](i, 0);
    }
    for(size_t i=v_len; i<t_len; i++){
        u_si.mShares[0](i, 0) = 0;
        u_si.mShares[1](i, 0) = 0;
    }

    //---------------------------------apply permutation start
    //step three: 对u进行sigma^{-1}的permutation
        //step0: si64matrix利用Sh3Converter转成sbmatrix
    Sh3Converter convt;
    convt.init(runtime, enc.mShareGen);
    Sh3Task task = runtime.noDependencies();

    sbMatrix sigma,u;
    convt.toBinaryMatrix(task, sigma_si, sigma).get();
    convt.toBinaryMatrix(task, u_si, u).get();

    //sigma_len=u_len=len
    size_t len = sigma.rows();
    size_t bitsize = sigma.bitCount();
    size_t unit_size = (sigma.bitCount() + 63) / 64;

    sbMatrix rsigma(len, bitsize);

        //step1: each party get pi_{i} and pi_{i+1} ;pi_{i}^{-1} and pi_{i+1}^{-1}
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();

    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);
    std::vector<size_t> prev_inverse_permutation;
    std::vector<size_t> next_inverse_permutation;
    get_inverse_permutation(prev_permutation, prev_inverse_permutation);
    get_inverse_permutation(next_permutation, next_inverse_permutation);
    
        //Step 2:对 [[sigma]]进行pi=p2*p0*p1序变换
    //  generate the random masks Z.
    i64Matrix prev_maskZ(len, unit_size);
    i64Matrix next_maskZ(len, unit_size);
    get_random_mask(pIdx, prev_maskZ, prevSeed);
    get_random_mask(pIdx, next_maskZ, nextSeed);
    
 
    if (pIdx == 0) {
        i64Matrix maskB(len, unit_size);
        i64Matrix maskA(len, unit_size);
        get_random_mask(pIdx, maskB, nextSeed);
        get_random_mask(pIdx, maskA, prevSeed);
        
        // 计算sharedX1 = sigma ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedX1(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX1(i, 0) = sigma.mShares[0](i, 0) ^ sigma.mShares[1](i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedX1);
        
        // 计算sharedX2 = sharedX1 ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedX2(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX2(i, 0) = sharedX1(i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedX2);
        
        // 发送sharedX2给P1
        large_data_sending(pIdx, sharedX2, runtime, true);
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            rsigma.mShares[1](i, 0) = maskA(i, 0);
            rsigma.mShares[0](i, 0) = maskB(i, 0);
        }
    }
    else if (pIdx == 1) {
        
        i64Matrix maskB(len, unit_size);
        get_random_mask(pIdx, maskB, prevSeed);
        
        // 计算sharedY1 = sigma ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedY1(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY1(i, 0) = sigma.mShares[0](i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedY1);
        
        // 发送sharedY1给P2
        large_data_sending(pIdx, sharedY1, runtime, true);
        // 接收来自P0的sharedX2
        i64Matrix sharedX2(len, unit_size);
        large_data_receiving(pIdx, sharedX2, runtime, true);
        
        // 计算sharedX3 = sharedX2 ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedX3(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX3(i, 0) = sharedX2(i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedX3);
        
        // 计算maskedC1 = sharedX3 ⊕ maskB
        i64Matrix maskedC1(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskedC1(i, 0) = sharedX3(i, 0) ^ maskB(i, 0);
        }
        large_data_sending(pIdx, maskedC1, runtime, true);
        
        // 接收来自P2的maskedC2
        i64Matrix maskedC2(len, unit_size);
        large_data_receiving(pIdx, maskedC2, runtime, true);
        
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            rsigma.mShares[1](i, 0) = maskB(i, 0);
            rsigma.mShares[0](i, 0) = maskedC1(i, 0) ^ maskedC2(i, 0);
        }
    }
    else if (pIdx == 2) {
        // maskA=maskz0
        i64Matrix maskA(len, unit_size);
        get_random_mask(pIdx, maskA, nextSeed);
        
        // 接收来自P1的sharedY1
        i64Matrix sharedY1(len, unit_size);
        large_data_receiving(pIdx, sharedY1, runtime, true);
        
        // 计算sharedY2 = sharedY1 ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedY2(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY2(i, 0) = sharedY1(i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedY2);
        
        // 计算sharedY3 = sharedY2 ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedY3(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY3(i, 0) = sharedY2(i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedY3);
        
        // 计算maskedC2 = sharedY3 ⊕ maskA
        i64Matrix maskedC2(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskedC2(i, 0) = sharedY3(i, 0) ^ maskA(i, 0);
        }
        large_data_sending(pIdx, maskedC2, runtime, true);
        
        // 接收来自P1的maskedC1
        i64Matrix maskedC1(len, unit_size);
        large_data_receiving(pIdx, maskedC1, runtime, true);
        
        // 计算maskC = maskedC1 ⊕ maskedC2
        i64Matrix maskC(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskC(i, 0) = maskedC1(i, 0) ^ maskedC2(i, 0);
        }
        
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            rsigma.mShares[1](i, 0) = maskC(i, 0);
            rsigma.mShares[0](i, 0) = maskA(i, 0);
        }
    }
    
        // Step 3-0: 恢复rsigma至明文
        // Step 3-1: 得到rsigma^{-1}
    i64Matrix rsigma_plain(len, unit_size);
    for (size_t i = 0; i < len; i++) {
        enc.revealAll(runtime, rsigma, rsigma_plain).get();
    }
    i64Matrix rsigma_inverse_plain(len, unit_size);
    permutation_inverse(rsigma_plain, rsigma_inverse_plain);
    
    //Step 4: 三方协同进行u1=pi([[u]])
    sbMatrix u1(len, bitsize);
    if (pIdx == 0) {
        i64Matrix maskB(len, unit_size);
        i64Matrix maskA(len, unit_size);
        get_random_mask(pIdx, maskB, nextSeed);
        get_random_mask(pIdx, maskA, prevSeed);
        
        // 计算sharedX1 = u[0] ⊕ u[1] ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedX1(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX1(i, 0) = u.mShares[0](i, 0) ^ u.mShares[1](i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedX1);
        
        // 计算sharedX2 = sharedX1 ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedX2(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX2(i, 0) = sharedX1(i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedX2);
        
        // 发送sharedX2给P1
        large_data_sending(pIdx, sharedX2, runtime, true);
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            u1.mShares[1](i, 0) = maskA(i, 0);
            u1.mShares[0](i, 0) = maskB(i, 0);
        }
    }
    else if (pIdx == 1) {
        
        i64Matrix maskB(len, unit_size);
        get_random_mask(pIdx, maskB, prevSeed);
        
        // 计算sharedY1 = u[2] ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedY1(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY1(i, 0) = u.mShares[0](i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedY1);
        
        // 发送sharedY1给P2
        large_data_sending(pIdx, sharedY1, runtime, true);
        // 接收来自P0的sharedX2
        i64Matrix sharedX2(len, unit_size);
        large_data_receiving(pIdx, sharedX2, runtime, true);
        
        // 计算sharedX3 = sharedX2 ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedX3(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX3(i, 0) = sharedX2(i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedX3);
        
        // 计算maskedC1 = sharedX3 ⊕ maskB
        i64Matrix maskedC1(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskedC1(i, 0) = sharedX3(i, 0) ^ maskB(i, 0);
        }
        large_data_sending(pIdx, maskedC1, runtime, true);
        
        // 接收来自P2的maskedC2
        i64Matrix maskedC2(len, unit_size);
        large_data_receiving(pIdx, maskedC2, runtime, true);
        
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            u1.mShares[1](i, 0) = maskB(i, 0);
            u1.mShares[0](i, 0) = maskedC1(i, 0) ^ maskedC2(i, 0);
        }
    }
    else if (pIdx == 2) {
        // maskA=maskz0
        i64Matrix maskA(len, unit_size);
        get_random_mask(pIdx, maskA, nextSeed);
        
        // 接收来自P1的sharedY1
        i64Matrix sharedY1(len, unit_size);
        large_data_receiving(pIdx, sharedY1, runtime, true);
        
        // 计算sharedY2 = sharedY1 ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedY2(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY2(i, 0) = sharedY1(i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedY2);
        
        // 计算sharedY3 = sharedY2 ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedY3(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY3(i, 0) = sharedY2(i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedY3);
        
        // 计算maskedC2 = sharedY3 ⊕ maskA
        i64Matrix maskedC2(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskedC2(i, 0) = sharedY3(i, 0) ^ maskA(i, 0);
        }
        large_data_sending(pIdx, maskedC2, runtime, true);
        
        // 接收来自P1的maskedC1
        i64Matrix maskedC1(len, unit_size);
        large_data_receiving(pIdx, maskedC1, runtime, true);
        
        // 计算maskC = maskedC1 ⊕ maskedC2
        i64Matrix maskC(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskC(i, 0) = maskedC1(i, 0) ^ maskedC2(i, 0);
        }
        
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            u1.mShares[1](i, 0) = maskC(i, 0);
            u1.mShares[0](i, 0) = maskA(i, 0);
        }
    }

        // Step 5: 根据rsigma^{-1}对v1进行置换得到v_prime = rsigma^{-1} (v1)=sigma^{-1} * pi^{-1} * pi([[v]]) = sigma^{-1}(v)
    sbMatrix u_prime(len, bitsize);
    for (size_t i = 0; i < len; i++) {
        size_t new_pos = rsigma_inverse_plain(i, 0);
        if (new_pos < len) {
            u_prime.mShares[0](new_pos, 0) = u1.mShares[0](i, 0);
            u_prime.mShares[1](new_pos, 0) = u1.mShares[1](i, 0);
        }
    } 
    si64Matrix u_prime_si(len, v.cols());
    convt.bitInjection(task, u_prime, u_prime_si).get();
    //---------------------------------apply permutation finished

    //step four: prefixsum_(u_prime)
    si64Matrix u_prime_prefix_si(len, v.cols());
    prefixsum(pIdx, u_prime_si, u_prime_prefix_si);

    //step five: unapply permutation
    //---------------------------------unapply permutation start
        //step 0:u_prime_prefix转成sbmatrix
    sbMatrix u_prime_prefix(len, bitsize);
    convt.toBinaryMatrix(task, u_prime_prefix_si, u_prime_prefix).get();

        //step 1:[[u1]]=rsigma(u')
    sbMatrix u_1_(len, bitsize);
    for (size_t i = 0; i < len; i++) {
        size_t new_pos = rsigma_plain(i, 0);
        if (new_pos < len) {
            u_1_.mShares[0](new_pos, 0) = u_prime_prefix.mShares[0](i, 0);
            u_1_.mShares[1](new_pos, 0) = u_prime_prefix.mShares[1](i, 0);
        }
    } 
       //step 2: [[u']]=pi^{-1}([[u1]])
    sbMatrix u_prime_(len, bitsize);
    //pi^{-1}=p1^{-1}*p0^{-1}*p2^{-1}，故pidx=0和2互换
    if (pIdx == 2) {
        i64Matrix maskB(len, unit_size);
        i64Matrix maskA(len, unit_size);
        get_random_mask(pIdx, maskB, nextSeed);
        get_random_mask(pIdx, maskA, prevSeed);
        
        // 计算sharedX1 = u_1 ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedX1(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX1(i, 0) = u_1_.mShares[0](i, 0) ^ u_1_.mShares[1](i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(prev_inverse_permutation, sharedX1);
        
        // 计算sharedX2 = sharedX1 ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedX2(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX2(i, 0) = sharedX1(i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(next_inverse_permutation, sharedX2);
        
        // 发送sharedX2给P1
        large_data_sending(pIdx, sharedX2, runtime, true);
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            u_prime_.mShares[1](i, 0) = maskA(i, 0);
            u_prime_.mShares[0](i, 0) = maskB(i, 0);
        }
    }
    else if (pIdx == 1) {
        
        i64Matrix maskB(len, unit_size);
        get_random_mask(pIdx, maskB, prevSeed);
        
        // 计算sharedY1 = u_1 ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedY1(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY1(i, 0) = u_1_.mShares[0](i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(next_inverse_permutation, sharedY1);
        
        // 发送sharedY1给P2
        large_data_sending(pIdx, sharedY1, runtime, true);
        // 接收来自P0的sharedX2
        i64Matrix sharedX2(len, unit_size);
        large_data_receiving(pIdx, sharedX2, runtime, true);
        
        // 计算sharedX3 = sharedX2 ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedX3(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX3(i, 0) = sharedX2(i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(prev_inverse_permutation, sharedX3);
        
        // 计算maskedC1 = sharedX3 ⊕ maskB
        i64Matrix maskedC1(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskedC1(i, 0) = sharedX3(i, 0) ^ maskB(i, 0);
        }
        large_data_sending(pIdx, maskedC1, runtime, true);
        
        // 接收来自P2的maskedC2
        i64Matrix maskedC2(len, unit_size);
        large_data_receiving(pIdx, maskedC2, runtime, true);
        
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            u_prime_.mShares[1](i, 0) = maskB(i, 0);
            u_prime_.mShares[0](i, 0) = maskedC1(i, 0) ^ maskedC2(i, 0);
        }
    }
    else if (pIdx == 0) {
        // maskA=maskz0
        i64Matrix maskA(len, unit_size);
        get_random_mask(pIdx, maskA, nextSeed);
        
        // 接收来自P1的sharedY1
        i64Matrix sharedY1(len, unit_size);
        large_data_receiving(pIdx, sharedY1, runtime, true);
        
        // 计算sharedY2 = sharedY1 ⊕ next_maskZ，然后应用next_permutation
        i64Matrix sharedY2(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY2(i, 0) = sharedY1(i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(prev_inverse_permutation, sharedY2);
        
        // 计算sharedY3 = sharedY2 ⊕ prev_maskZ，然后应用prev_permutation
        i64Matrix sharedY3(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY3(i, 0) = sharedY2(i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(next_inverse_permutation, sharedY3);
        
        // 计算maskedC2 = sharedY3 ⊕ maskA
        i64Matrix maskedC2(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskedC2(i, 0) = sharedY3(i, 0) ^ maskA(i, 0);
        }
        large_data_sending(pIdx, maskedC2, runtime, true);
        
        // 接收来自P1的maskedC1
        i64Matrix maskedC1(len, unit_size);
        large_data_receiving(pIdx, maskedC1, runtime, true);
        
        // 计算maskC = maskedC1 ⊕ maskedC2
        i64Matrix maskC(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskC(i, 0) = maskedC1(i, 0) ^ maskedC2(i, 0);
        }
        
        // 计算最终份额
        for (size_t i = 0; i < len; i++) {
            u_prime_.mShares[1](i, 0) = maskC(i, 0);
            u_prime_.mShares[0](i, 0) = maskA(i, 0);
        }
    }
    si64Matrix u_prime_si_(len, v.cols());
    convt.bitInjection(task,u_prime_, u_prime_si_).get();
    //---------------------------------unapply permutation finished

    //step six:get result
    result.resize(idx_len, v.cols());
    for (size_t i = 0; i < idx_len; i++) {
        result.mShares[0](i, 0) = u_prime_si_.mShares[0](i+v_len, 0);
        result.mShares[1](i, 0) = u_prime_si_.mShares[1](i+v_len, 0);
    }

    return; 
}





