#include "Matrix.h"

using namespace oc;
using namespace aby3;

int fixed_matrix_mult(sf64Matrix<D8>& A, sf64Matrix<D8>& B, sf64Matrix<D8>& C, int pIdx, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime) {
    if (A.cols() != B.rows() || A.rows() != C.rows() || B.cols() != C.cols()) {
        THROW_RUNTIME_ERROR("Matrix dimensions do not match.");
    }

    sf64Matrix<D8> paddedA(A.rows() * A.cols(), 1);
    for(u64 i = 0; i < B.cols(); ++i) {
        sf64Matrix<D8> repeatCol(A.rows() * A.cols(), 1);
        for(u64 j = 0; j < A.rows(); ++j)
            for(u64 k = 0; k < A.cols(); ++k)
                repeatCol(j * A.cols() + k, 0, B(k, i));
        sf64Matrix<D8> sharedProd;
        eval.asyncMul(runtime, paddedA, repeatCol, sharedProd).get();
        for(u64 j = 0; j < A.rows(); ++j) {
            sf64<D8> sum = sharedProd(j * A.cols(), 0);
            for(u64 k = 1; k < A.cols(); ++k)
                sum = sum + sharedProd(j * A.cols() + k, 0);
            C(j, i, sum);
        }
    }
    
    return 0;
}

int fixed_matrix_mult_shift(sf64Matrix<D8>& A, sf64Matrix<D8>& B, sf64Matrix<D8>& C, int shift, int pIdx, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime) {
    if (A.cols() != B.rows() || A.rows() != C.rows() || B.cols() != C.cols()) {
        THROW_RUNTIME_ERROR("Matrix dimensions do not match.");
    }

    sf64Matrix<D8> paddedA(A.rows() * A.cols(), 1);
    for(u64 i = 0; i < B.cols(); ++i) {
        sf64Matrix<D8> repeatCol(A.rows() * A.cols(), 1);
        for(u64 j = 0; j < A.rows(); ++j)
            for(u64 k = 0; k < A.cols(); ++k)
                repeatCol(j * A.cols() + k, 0, B(k, i));
        sf64Matrix<D8> sharedProd;
        eval.asyncMul(runtime, paddedA, repeatCol, sharedProd, shift).get();
        for(u64 j = 0; j < A.rows(); ++j) {
            sf64<D8> sum = sharedProd(j * A.cols(), 0);
            for(u64 k = 1; k < A.cols(); ++k)
                sum = sum + sharedProd(j * A.cols() + k, 0);
            C(j, i, sum);
        }
    }
    
    return 0;
}

int fixed_matrix_sum(sf64Matrix<D8>& A, sf64Matrix<D8>& B, sf64Matrix<D8>& C, int pIdx, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime) {
    if (A.rows() != B.rows() || A.cols() != B.cols() || A.rows() != C.rows() || A.cols() != C.cols()) {
        THROW_RUNTIME_ERROR("Matrix dimensions do not match.");
    }

    C = A + B;
    
    return 0;
}
