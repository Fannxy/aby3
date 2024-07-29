#include "Matrix.h"

using namespace oc;
using namespace aby3;

int fixed_matrix_mult(sf64Matrix<D8>& A, sf64Matrix<D8>& B, sf64Matrix<D8>& C, int pIdx, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime) {
    if (A.cols() != B.rows() || A.rows() != C.rows() || B.cols() != C.cols()) {
        THROW_RUNTIME_ERROR("Matrix dimensions do not match.");
    }

    std::vector<sf64Matrix<D8>> sharedRow(A.rows()),
                                sharedCol(B.cols());

    for(u64 i = 0; i < A.rows(); ++i) {
        sharedRow[i] = sf64Matrix<D8>(A.cols(), 1);
        for(u64 j = 0; j < A.cols(); ++j)
            sharedRow[i](j, 0, A(i, j));
    }

    for(u64 i = 0; i < B.cols(); ++i) {
        sharedCol[i] = sf64Matrix<D8>(A.cols(), 1);
        for(u64 j = 0; j < A.cols(); ++j)
            sharedCol[i](j, 0, B(j, i));
    }

    for(u64 i = 0; i < A.rows(); ++i)
        for(u64 j = 0; j < B.cols(); ++j) {
            sf64Matrix<D8> sharedProd;
            eval.asyncMul(runtime, sharedRow[i], sharedCol[j], sharedProd).get();
            sf64<D8> sum = sharedProd(0, 0);
            for(u64 k = 1; k < A.cols(); ++k)
                sum = sum + sharedProd(k, 0);
            C(i, j, sum);
        }
    
    return 0;
}

// int fixed_matrix_mult(sf64Matrix<D8>& A, sf64Matrix<D8>& B, sf64Matrix<D8>& C, int pIdx, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime) {
//     if (A.cols() != B.rows() || A.rows() != C.rows() || B.cols() != C.cols()) {
//         THROW_RUNTIME_ERROR("Matrix dimensions do not match.");
//     }

//     sf64Matrix<D8> repeatCol(A.rows(), A.cols());

//     for(u64 i = 0; i < B.cols(); ++i) {
//         for(u64 j = 0; j < A.rows(); ++j)
//             for(u64 k = 0; k < A.cols(); ++k)
//                 repeatCol(j, k, B(k, i));
//         sf64Matrix<D8> sharedProd;
//         eval.asyncMul(runtime, A, repeatCol, sharedProd).get();
//         for(u64 j = 0; j < A.rows(); ++j) {
//             sf64<D8> sum = sharedProd(j, 0);
//             for(u64 k = 1; k < A.cols(); ++k)
//                 sum = sum + sharedProd(j, k);
//             C(j, i, sum);
//         }
//     }
    
//     return 0;
// }