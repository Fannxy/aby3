#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-Basic/Matrix.h"
#include "../aby3-RTR/BuildingBlocks.h"


using namespace oc;
using namespace aby3;

int fixed_matrix_mult_test(oc::CLP& cmd) {

    BASIC_TEST_INIT

    u64 sizeX = 16,
        sizeY = cmd.getMany<int>("dataSize")[0],
        sizeZ = 16;
    
    f64Matrix<D8> A(sizeX, sizeY);
    for (u64 i = 0; i < sizeX; ++i)
        for (u64 j = 0; j < sizeY; ++j)
            A(i, j) = double(i + 1) / (i + j + 1);
    
    f64Matrix<D8> B(sizeY, sizeZ);
    for (u64 i = 0; i < sizeY; ++i)
        for (u64 j = 0; j < sizeZ; ++j)
            B(i, j) = double(i + 1) / (i + j + 1);
    
    sf64Matrix<D8> sharedA(sizeX, sizeY),
                sharedB(sizeY, sizeZ),
                sharedC(sizeX, sizeZ);
    if(role == 0) {
        enc.localFixedMatrix(runtime, A, sharedA).get();
        enc.localFixedMatrix(runtime, B, sharedB).get();
    } else {
        enc.remoteFixedMatrix(runtime, sharedA).get();
        enc.remoteFixedMatrix(runtime, sharedB).get();
    }

    fixed_matrix_mult(sharedA, sharedB, sharedC, 0, enc, eval, runtime);

    f64Matrix<D8> C(sizeX, sizeZ);
    enc.revealAll(runtime, sharedC, C).get();

    if(role == 0) {
        for(u64 i = 0; i < sizeX; ++i) {
            for(u64 j = 0; j < sizeZ; ++j) {
                std::cerr << C(i, j) << " ";
            }
            std::cerr << std::endl;
        }
    }
    
    return 0;
}

int splitted_fixed_matrix_mult_test(oc::CLP& cmd) {

    SPLITTED_TEST_INIT

    u64 sizeX = 16,
        sizeY = cmd.getMany<int>("dataSize")[0],
        sizeZ = 16;
    
    f64Matrix<D8> A(sizeX, sizeY);
    for (u64 i = 0; i < sizeX; ++i)
        for (u64 j = 0; j < sizeY; ++j)
            A(i, j) = double(i + 1) / (i + j + 1);
    
    f64Matrix<D8> B(sizeY, sizeZ);
    for (u64 i = 0; i < sizeY; ++i)
        for (u64 j = 0; j < sizeZ; ++j)
            B(i, j) = double(i + 1) / (i + j + 1);
    
    sf64Matrix<D8> sharedA(sizeX, sizeY),
                sharedB(sizeY, sizeZ),
                sharedC(sizeX, sizeZ);
    if(role == 0) {
        enc.localFixedMatrix(runtime, A, sharedA).get();
        enc.localFixedMatrix(runtime, B, sharedB).get();
    } else {
        enc.remoteFixedMatrix(runtime, sharedA).get();
        enc.remoteFixedMatrix(runtime, sharedB).get();
    }

    fixed_matrix_mult(sharedA, sharedB, sharedC, 0, enc, eval, runtime);

    f64Matrix<D8> C(sizeX, sizeZ);
    enc.revealAll(runtime, sharedC, C).get();
    
    return 0;
}