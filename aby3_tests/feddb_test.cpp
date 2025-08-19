#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Feddb-Core/feddb.h"
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
    debug_info("initial finished");

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
    debug_info("oblivious_idx_select finished");

    i64Matrix res_test(idx_size, 1);
    enc.revealAll(runtime, res_shared, res_test).get();

    if(role == 0){
        check_result("OblIdx Test", res_test, res_plain);
    }
    
    return 0;
}