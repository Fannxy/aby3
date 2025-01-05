#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-Basic/Matrix.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-Basic/Shuffle.h"
#include "../aby3-Basic/timer.h"
#include "../aby3-RTR/GASTest.h"



using namespace oc;
using namespace aby3;

int get_si64_shares(u64 data_size, aby3::si64Matrix& data){
    data.resize(data_size, 1);
    for(i64 i=0; i<data_size; i++){
        data.mShares[0](i) = 0;
        data.mShares[1](i) = 0;
    }
    return 0;
}

int get_sb_shares(u64 data_size, aby3::sbMatrix& data){
    data.resize(data_size, 64);
    for(i64 i=0; i<data_size; i++){
        data.mShares[0](i) = 0;
        data.mShares[1](i) = 0;
    }
    return 0;
}

int splitted_micro_benchmarks(oc::CLP& cmd){
    SPLITTED_TEST_INIT

    u64 sizeX = cmd.getMany<int>("dataSize")[0];
    std::string task = "";

    if(cmd.isSet("ii-mul")){
        task = "ii-mul";
        aby3::si64Matrix dataX, dataY, dataZ;
        get_si64_shares(sizeX, dataX);
        get_si64_shares(sizeX, dataY);
        dataZ.resize(sizeX, 1);
        cipher_mul(role, dataX, dataY, dataZ, eval, enc, runtime);
    }

    if(cmd.isSet("ff-mul")){
        task = "ff-mul";
        aby3::sf64Matrix<D8> dataX, dataY, dataZ;
        get_sf64_shares<D8>(sizeX, dataX);
        get_sf64_shares<D8>(sizeX, dataY);
        dataZ.resize(sizeX, 1);
        cipher_mul<D8>(role, dataX, dataY, dataZ, eval, enc, runtime);
    }

    if(cmd.isSet("ib-mul")){
        task = "ib-mul";
        aby3::si64Matrix dataX, dataZ;
        aby3::sbMatrix dataY;
        get_si64_shares(sizeX, dataX);  
        get_sb_shares(sizeX, dataY);
        dataY.resize(sizeX, 1);
        dataZ.resize(sizeX, 1);
        cipher_mul(role, dataX, dataY, dataZ, eval, enc, runtime);
    }

    if(cmd.isSet("a-gt")){
        task = "a-gt";
        aby3::si64Matrix dataX, dataY;
        aby3::sbMatrix dataZ;
        get_si64_shares(sizeX, dataX);
        get_si64_shares(sizeX, dataY);
        dataZ.resize(sizeX, 1);
        cipher_gt(role, dataX, dataY, dataZ, eval, runtime);
    }

    if(cmd.isSet("a2b")){
        task = "a2b";
        aby3::si64Matrix dataX;
        aby3::sbMatrix dataY(sizeX, 64);    
        get_si64_shares(sizeX, dataX);
        arith2bool(role, dataX, dataY, enc, eval, runtime);
    }

    if(cmd.isSet("b2a")){
        task = "b2a";
        aby3::sbMatrix dataX;
        aby3::si64Matrix dataY(sizeX, 1);
        get_sb_shares(sizeX, dataX);
        bool2arith(role, dataX, dataY, enc, eval, runtime);
    }

    if(cmd.isSet("b2a-single")){
        task = "b2a-single";
        aby3::sbMatrix dataX;
        aby3::si64Matrix dataY(sizeX, 1);
        get_sb_shares(sizeX, dataX);
        dataX.resize(sizeX, 1);
        bool2arith(role, dataX, dataY, enc, eval, runtime);
    }

    if(cmd.isSet("shuffle")){
        task = "shuffle";
        aby3::sbMatrix dataX, dataRes;
        get_sb_shares(sizeX, dataX);
        dataRes.resize(sizeX, 64);
        efficient_shuffle(dataX, role, dataRes, enc, eval, runtime);
    }

    if(cmd.isSet("fake_test")){
        task = "fake";
        aby3::i64Matrix dataX(sizeX, 1);
        for(u64 i=0; i<sizeX; i++) dataX(i, 0) = i;
        fake_test(role, dataX, runtime);
    }

    if(cmd.isSet("fake_test2")){
        task = "fake2";
        aby3::i64Matrix dataX(sizeX, 1);
        for(u64 i=0; i<sizeX; i++) dataX(i, 0) = i;
        fake_test2(role, dataX, runtime);
    }

    Timer& timer = Timer::getInstance();
    std::string stamp_file = "/tmp/aby3-stamp.txt";
    get_value("stampFile", cmd, stamp_file);
    std::ofstream stamp(stamp_file, std::ios::app);
    timer.print_time_stamps(stamp);
    stamp.close();

    return 0;
}