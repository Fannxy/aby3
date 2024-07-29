#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-RTR/BuildingBlocks.h"


using namespace oc;
using namespace aby3;

int int_matrix_multiplication_test(oc::CLP& cmd) {

    BASIC_TEST_INIT

    u64 sizeX = 40,
        sizeY = 60,
        sizeZ = 70;
    
    eMatrix<i64> plainMatrix(sizeX, sizeY);
    for (u64 i = 0; i < sizeX; ++i)
        for (u64 j = 0; j < sizeY; ++j)
            plainMatrix(i, j) = i + j;
    
    eMatrix<i64> plainMatrix2(sizeY, sizeZ);
    for (u64 i = 0; i < sizeY; ++i)
        for (u64 j = 0; j < sizeZ; ++j)
            plainMatrix2(i, j) = i + j;

    eMatrix<i64> prodMtx(sizeX, sizeZ);
    for(u64 i = 0; i < sizeX; ++i)
        for(u64 j = 0; j < sizeZ; ++j) {
            eMatrix<i64> row(sizeY, 1),
                        col(sizeY, 1);
            for(u64 k = 0; k < sizeY; ++k) {
                row(k, 0) = plainMatrix(i, k);
                col(k, 0) = plainMatrix2(k, j);
            }

            si64Matrix sharedRow(sizeY, 1),
                            sharedCol(sizeY, 1);
            if(role == 0) {
                enc.localIntMatrix(runtime, row, sharedRow).get();
                enc.localIntMatrix(runtime, col, sharedCol).get();
            } else {
                enc.remoteIntMatrix(runtime, sharedRow).get();
                enc.remoteIntMatrix(runtime, sharedCol).get();
            }

            si64Matrix sharedProd;
            eval.asyncMul(runtime, sharedRow, sharedCol, sharedProd).get();

            si64 sum = sharedProd(0, 0);
            for(u64 k = 1; k < sizeY; ++k)
                sum = sum + sharedProd(k, 0);
            enc.revealAll(runtime, sum, prodMtx(i, j)).get();
        }

    return 0;
}

// int fixed_matrix_multiplication_test(oc::CLP& cmd) {

//     BASIC_TEST_INIT

//     u64 sizeX = cmd.getMany<int>("sizeX")[0],
//         sizeY = cmd.getMany<int>("sizeY")[0],
//         sizeZ = cmd.getMany<int>("sizeZ")[0];
    
//     f64Matrix<D8> fixedMatrix(sizeX, sizeY);
//     for (u64 i = 0; i < sizeX; ++i)
//         for (u64 j = 0; j < sizeY; ++j)
//             fixedMatrix(i, j) = double(i + 1) / (i + j + 1);
    
//     f64Matrix<D8> fixedMatrix2(sizeY, sizeZ);
//     for (u64 i = 0; i < sizeY; ++i)
//         for (u64 j = 0; j < sizeZ; ++j)
//             fixedMatrix2(i, j) = double(i + 1) / (i + j + 1);
    
//     std::vector<sf64Matrix<D8>> sharedRow(sizeX),
//                                 sharedCol(sizeZ);

//     for(u64 i = 0; i < sizeX; ++i) {
//         f64Matrix<D8> row(sizeY, 1);
//         for(u64 j = 0; j < sizeY; ++j)
//             row(j, 0) = fixedMatrix(i, j);
//         sharedRow[i] = sf64Matrix<D8>(sizeY, 1);
//         if(role == 0)
//             enc.localFixedMatrix(runtime, row, sharedRow[i]).get();
//         else
//             enc.remoteFixedMatrix(runtime, sharedRow[i]).get();
//     }

//     for(u64 i = 0; i < sizeZ; ++i) {
//         f64Matrix<D8> col(sizeY, 1);
//         for(u64 j = 0; j < sizeY; ++j)
//             col(j, 0) = fixedMatrix2(j, i);
//         sharedCol[i] = sf64Matrix<D8>(sizeY, 1);
//         if(role == 0)
//             enc.localFixedMatrix(runtime, col, sharedCol[i]).get();
//         else
//             enc.remoteFixedMatrix(runtime, sharedCol[i]).get();
//     }

//     f64Matrix<D8> prodMtx(sizeX, sizeZ);
//     for(u64 i = 0; i < sizeX; ++i)
//         for(u64 j = 0; j < sizeZ; ++j) {
//             sf64Matrix<D8> sharedProd;
//             eval.asyncMul(runtime, sharedRow[i], sharedCol[j], sharedProd).get();
//             sf64<D8> sum = sharedProd(0, 0);
//             for(u64 k = 1; k < sizeY; ++k)
//                 sum = sum + sharedProd(k, 0);
//             enc.revealAll(runtime, sum, prodMtx(i, j)).get();
//         }
    
//     return 0;
// }

int fixed_matrix_multiplication_test(oc::CLP& cmd) {

    BASIC_TEST_INIT

    u64 sizeX = cmd.getMany<int>("sizeX")[0],
        sizeY = cmd.getMany<int>("sizeY")[0],
        sizeZ = cmd.getMany<int>("sizeZ")[0];
    
    f64Matrix<D8> fixedMatrix(sizeX, sizeY);
    for (u64 i = 0; i < sizeX; ++i)
        for (u64 j = 0; j < sizeY; ++j)
            fixedMatrix(i, j) = double(i + 1) / (i + j + 1);
    
    f64Matrix<D8> fixedMatrix2(sizeZ, sizeY);
    for (u64 i = 0; i < sizeZ; ++i)
        for (u64 j = 0; j < sizeY; ++j)
            fixedMatrix2(i, j) = double(j + 1) / (i + j + 1);
    
    sf64Matrix<D8> sharedRow(sizeX, sizeY);
    if(role == 0)
        enc.localFixedMatrix(runtime, fixedMatrix, sharedRow).get();
    else
        enc.remoteFixedMatrix(runtime, sharedRow).get();

    f64Matrix<D8> prodMtx(sizeX, sizeZ);
    for(u64 i = 0; i < sizeZ; ++i) {
        f64Matrix<D8> repeatedCol(sizeX, sizeY);
        for(u64 j = 0; j < sizeX; ++j)
            for(u64 k = 0; k < sizeY; ++k)
                repeatedCol(j, k) = fixedMatrix2(i, k);
        
        sf64Matrix<D8> sharedCol(sizeX, sizeY);
        if(role == 0)
            enc.localFixedMatrix(runtime, repeatedCol, sharedCol).get();
        else
            enc.remoteFixedMatrix(runtime, sharedCol).get();
        
        sf64Matrix<D8> sharedProd(sizeX, sizeY);
        std::cerr << "HERE" << std::endl;
        eval.asyncMul(runtime, sharedRow, sharedCol, sharedProd).get();
        for(u64 j = 0; j < sizeX; ++j) {
            sf64<D8> sum = sharedProd(j, 0);
            for(u64 k = 1; k < sizeY; ++k)
                sum = sum + sharedProd(j, k);
            enc.revealAll(runtime, sum, prodMtx(j, i)).get();
        }
    }

    // for(u64 i = 0; i < sizeX; ++i) {
    //     for(u64 j = 0; j < sizeZ; ++j) {
    //         std::cerr << prodMtx(i, j) << " ";
    //     }
    //     std::cerr << std::endl;
    // }

    return 0;
}

// int splitted_fixed_matrix_multiplication_test(oc::CLP& cmd) {

//     SPLITTED_TEST_INIT

//     u64 sizeX = cmd.getMany<int>("sizeX")[0],
//         sizeY = cmd.getMany<int>("sizeY")[0],
//         sizeZ = cmd.getMany<int>("sizeZ")[0];
    
//     f64Matrix<D8> fixedMatrix(sizeX, sizeY);
//     for (u64 i = 0; i < sizeX; ++i)
//         for (u64 j = 0; j < sizeY; ++j)
//             fixedMatrix(i, j) = double(i + 1) / (i + j + 1);
    
//     f64Matrix<D8> fixedMatrix2(sizeY, sizeZ);
//     for (u64 i = 0; i < sizeY; ++i)
//         for (u64 j = 0; j < sizeZ; ++j)
//             fixedMatrix2(i, j) = double(i + 1) / (i + j + 1);
    
//     std::vector<sf64Matrix<D8>> sharedRow(sizeX),
//                                 sharedCol(sizeZ);

//     for(u64 i = 0; i < sizeX; ++i) {
//         f64Matrix<D8> row(sizeY, 1);
//         for(u64 j = 0; j < sizeY; ++j)
//             row(j, 0) = fixedMatrix(i, j);
//         sharedRow[i] = sf64Matrix<D8>(sizeY, 1);
//         if(role == 0)
//             enc.localFixedMatrix(runtime, row, sharedRow[i]).get();
//         else
//             enc.remoteFixedMatrix(runtime, sharedRow[i]).get();
//     }

//     for(u64 i = 0; i < sizeZ; ++i) {
//         f64Matrix<D8> col(sizeY, 1);
//         for(u64 j = 0; j < sizeY; ++j)
//             col(j, 0) = fixedMatrix2(j, i);
//         sharedCol[i] = sf64Matrix<D8>(sizeY, 1);
//         if(role == 0)
//             enc.localFixedMatrix(runtime, col, sharedCol[i]).get();
//         else
//             enc.remoteFixedMatrix(runtime, sharedCol[i]).get();
//     }

//     f64Matrix<D8> prodMtx(sizeX, sizeZ);
//     for(u64 i = 0; i < sizeX; ++i) {
//         for(u64 j = 0; j < sizeZ; ++j) {
//             sf64Matrix<D8> sharedProd;
//             eval.asyncMul(runtime, sharedRow[i], sharedCol[j], sharedProd).get();
//             sf64<D8> sum = sharedProd(0, 0);
//             for(u64 k = 1; k < sizeY; ++k)
//                 sum = sum + sharedProd(k, 0);
//             enc.revealAll(runtime, sum, prodMtx(i, j)).get();
//         }
//         if(role == 0 && rank == 0 && i % 10 == 0) {
//             std::cerr << "Row " << i << " done" << std::endl;
//         }
//     }
    
//     return 0;
// }

int splitted_fixed_matrix_multiplication_test(oc::CLP& cmd) {

    SPLITTED_TEST_INIT

    u64 sizeX = cmd.getMany<int>("sizeX")[0],
        sizeY = cmd.getMany<int>("sizeY")[0],
        sizeZ = cmd.getMany<int>("sizeZ")[0];
    
    f64Matrix<D8> fixedMatrix(sizeX, sizeY);
    for (u64 i = 0; i < sizeX; ++i)
        for (u64 j = 0; j < sizeY; ++j)
            fixedMatrix(i, j) = double(i + 1) / (i + j + 1);
    
    f64Matrix<D8> fixedMatrix2(sizeZ, sizeY);
    for (u64 i = 0; i < sizeZ; ++i)
        for (u64 j = 0; j < sizeY; ++j)
            fixedMatrix2(i, j) = double(j + 1) / (i + j + 1);
    
    sf64Matrix<D8> sharedRow(sizeX, sizeY);
    if(role == 0)
        enc.localFixedMatrix(runtime, fixedMatrix, sharedRow).get();
    else
        enc.remoteFixedMatrix(runtime, sharedRow).get();

    f64Matrix<D8> prodMtx(sizeX, sizeZ);
    sf64Matrix<D8> sharedProdMtx(sizeX, sizeZ);
    f64Matrix<D8> repeatedCol(sizeX, sizeY);
    sf64Matrix<D8> sharedCol(sizeX, sizeY);
    for(u64 i = 0; i < sizeZ; ++i) {
        for(u64 j = 0; j < sizeX; ++j)
            for(u64 k = 0; k < sizeY; ++k)
                repeatedCol(j, k) = fixedMatrix2(i, k);
        
        if(role == 0)
            enc.localFixedMatrix(runtime, repeatedCol, sharedCol).get();
        else
            enc.remoteFixedMatrix(runtime, sharedCol).get();
        
        sf64Matrix<D8> sharedProd(sizeX, sizeY);
        std::cerr<<"HERE "<<role<<' '<<rank<<std::endl;
        eval.asyncMul(runtime, sharedRow, sharedCol, sharedProd).get();
        std::cerr<<"END "<<role<<' '<<rank<<std::endl;
        for(u64 j = 0; j < sizeX; ++j) {
            auto sum = sharedProd(j, 0);
            for(u64 k = 1; k < sizeY; ++k) {
                sum[0] = sum[0] + sharedProd(j, k)[0];
                sum[1] = sum[1] + sharedProd(j, k)[1];
            }
            sharedProdMtx(j, i)[0] = sum[0];
            sharedProdMtx(j, i)[1] = sum[1];
        }
    }
    enc.revealAll(runtime, sharedProdMtx, prodMtx).get();

    // for(u64 i = 0; i < sizeX; ++i) {
    //     for(u64 j = 0; j < sizeZ; ++j) {
    //         std::cerr << prodMtx(i, j) << " ";
    //     }
    //     std::cerr << std::endl;
    // }

    return 0;
}
