#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Feddb-Core/feddb.h"
#include "../aby3-Feddb-Core/genperm.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/debug.h"

using namespace oc;
using namespace aby3;
using namespace std;


int oblivious_idx_select_test(CLP &cmd) {
    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN OblIdx TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    size_t data_size = 5;
    size_t idx_size = 4;
    i64Matrix data_plain(data_size, 1);
    i64Matrix idx_plain(idx_size, 1);
    i64Matrix res_plain(idx_size, 1);

    data_plain(0,0)=1;
    data_plain(1,0)=4;
    data_plain(2,0)=9;
    data_plain(3,0)=16;
    data_plain(4,0)=25;

    idx_plain(0,0)=2;
    idx_plain(1,0)=1;
    idx_plain(2,0)=3;
    idx_plain(3,0)=2;

    res_plain(0,0)=9;
    res_plain(1,0)=4;
    res_plain(2,0)=16;
    res_plain(3,0)=9;

    // encrypt the inputs.
    si64Matrix data_shared(data_size, 1);
    si64Matrix idx_shared(idx_size, 1);
    si64Matrix res_shared(idx_size, 1);

    if (role == 0) {
        enc.localIntMatrix(runtime, data_plain, data_shared).get();
        enc.localIntMatrix(runtime, idx_plain, idx_shared).get();
    } else {
        enc.remoteIntMatrix(runtime, data_shared).get();
        enc.remoteIntMatrix(runtime, idx_shared).get();
    }

    oblivious_idx_select(role, data_shared, idx_shared, res_shared, enc, eval, runtime);

    i64Matrix res_test(idx_size, 1);
    enc.revealAll(runtime, res_shared, res_test).get();

    if(role == 0){
        check_result("OblIdx Test", res_test, res_plain);
    }
    
    return 0;
}

int genperm_test(CLP &cmd){
    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN GenPerm TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    size_t data_size = 9;
    size_t perm_size = 9;

    i64Matrix data_plain(data_size, 1);
    i64Matrix perm_plain(perm_size, 1);

    
    for(size_t i=0;i<5;i++){
        data_plain(i,0)=i;
    }

    data_plain(5,0)=2;
    data_plain(6,0)=1;
    data_plain(7,0)=3;
    data_plain(8,0)=2;

    perm_plain(0,0)=0;
    perm_plain(1,0)=1;
    perm_plain(2,0)=3;
    perm_plain(3,0)=6;
    perm_plain(4,0)=8;
    perm_plain(5,0)=4;
    perm_plain(6,0)=2;
    perm_plain(7,0)=7;
    perm_plain(8,0)=5;

    si64Matrix data_shared(data_size, 1);

    if (role == 0) {
        enc.localIntMatrix(runtime, data_plain, data_shared).get();
    } else {
        enc.remoteIntMatrix(runtime, data_shared).get();
    }

    si64Matrix perm_shared(perm_size, 1);
    genPerm(role, data_shared, perm_shared, enc, eval, runtime);

    i64Matrix perm_test(perm_size, 1);
    enc.revealAll(runtime, perm_shared, perm_test).get();

    //DEBUG
    // if(role == 0){
    //     std::cout << "perm_test: " << std::endl;
    //     for(size_t i=0;i<perm_size;i++){
    //         std::cout << perm_test(i,0) << " ";
    //     }
    //     std::cout << std::endl ;
    // }

    if(role == 0){
        check_result("GenPerm Test", perm_test, perm_plain);
    }

    return 0;

}