#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-Basic/Matrix.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-Basic/Shuffle.h"
#include "../aby3-Basic/timer.h"

using namespace oc;
using namespace aby3;

int fake_test(int role, i64Matrix& dataX, Sh3Runtime& runtime){

    u64 sizeX = dataX.size();
    i64Matrix buffer(sizeX, 1);
    if(role == 0){
        runtime.mComm.mNext.asyncSendCopy(dataX.data(), dataX.size());
        runtime.mComm.mNext.recv(buffer.data(), dataX.size());

    }
    if(role == 1){
        runtime.mComm.mPrev.asyncSendCopy(dataX.data(), dataX.size());
        runtime.mComm.mPrev.recv(buffer.data(), dataX.size());
    }

    return 0;
}

int fake_test2(int role, i64Matrix& dataX, Sh3Runtime& runtime){

    u64 sizeX = dataX.size();
    i64Matrix buffer(sizeX, 1);
    if(role == 0){
        runtime.mComm.mNext.asyncSendCopy(dataX.data(), dataX.size());
        // runtime.mComm.mNext.recv(buffer.data(), dataX.size());

    }
    if(role == 1){
        // runtime.mComm.mPrev.asyncSendCopy(dataX.data(), dataX.size());
        runtime.mComm.mPrev.recv(buffer.data(), dataX.size());
    }

    return 0;
}