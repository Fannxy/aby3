#include "../aby3-Benchmark/benchmark.h"
#include "../aby3-Basic/Search.h"

using namespace oc;
using namespace aby3;

#define SET_OR_DEFAULT_SIZET(cmd, key, default_value) \
    size_t key = default_value; \
    if(cmd.isSet(#key)){ \
        auto keys = cmd.getMany<size_t>(#key); \
        key = keys[0]; \
    } \
    else{ \
        debug_info(#key " is not set, using default value: " + std::to_string(key)); \
    }

#define SET_OR_DEFAULT_INT(cmd, key, default_value) \
    int key = default_value; \
    if(cmd.isSet(#key)){ \
        auto keys = cmd.getMany<int>(#key); \
        key = keys[0]; \
    } \
    else{ \
        debug_info(#key " is not set, using default value: " + std::to_string(key)); \
    }

#define SET_OR_DEFAULT_STRING(cmd, key, default_value) \
    std::string key = default_value; \
    if(cmd.isSet(#key)){ \
        auto keys = cmd.getMany<std::string>(#key); \
        key = keys[0]; \
    } \
    else{ \
        debug_info(#key " is not set, using default value: " + key); \
    }

void generate_data(std::vector<aby3::si64Matrix> &data, size_t n, size_t m, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    i64Matrix random_data(n*m, 1);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            random_data(i*m + j, 0) = rand() % 100;
        }
    }
    aby3::si64Matrix full_data(n*m, 1);
    large_data_encryption(pIdx, random_data, full_data, enc, runtime);

    for(int i = 0; i < n; ++i){
        data[i].resize(m, 1);
        for(int j = 0; j < m; ++j){
            data[i].mShares[0](j, 0) = full_data.mShares[0](i*m + j, 0);
            data[i].mShares[1](j, 0) = full_data.mShares[1](i*m + j, 0);
        }
    }
    return;
}


void generate_data(std::vector<aby3::sbMatrix> &data, size_t n, size_t m, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    i64Matrix random_data(n*m, 1);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            random_data(i*m + j, 0) = rand() % 100;
        }
    }
    aby3::sbMatrix full_data(n*m, bitsize);
    large_data_encryption(pIdx, random_data, full_data, enc, runtime);

    for(int i = 0; i < n; ++i){
        data[i].resize(m, bitsize);
        for(int j = 0; j < m; ++j){
            data[i].mShares[0](j, 0) = full_data.mShares[0](i*m + j, 0);
            data[i].mShares[1](j, 0) = full_data.mShares[1](i*m + j, 0);
        }
    }
    return;
}


void binary_search_benchmark(oc::CLP& cmd){

    // get the configs.
    #ifdef MPI_APP
    SETUP_PROCESS
    #endif
    #ifndef MPI_APP
    SETUP_SINGLE_PROCESS
    #endif

    // get the corresponding parameters.
    SET_OR_DEFAULT_INT(cmd, N, 1<<10);
    SET_OR_DEFAULT_STRING(cmd, Func, "mcomp");
    SET_OR_DEFAULT_STRING(cmd, record_file, "/root/aby3/record.txt");

    // get the data.
    std::vector<aby3::si64Matrix> data(N);
    std::vector<aby3::sbMatrix> key(N);
    generate_data(data, N, 1, role, enc, eval, runtime);
    generate_data(key, N, 1, 32, role, enc, eval, runtime);

    boolIndex _target_key(N/2, role);
    aby3::sbMatrix target_key = _target_key.to_matrix(32);

    aby3::sbMatrix result;
    aby3::si64Matrix result_data;

    // benchmark the search function.
    Timer& timer = Timer::getInstance();
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    if(Func == "mcomp"){
        timer.start("mcomp");
        mcompBS(key, target_key, result, role, enc, eval, runtime);
        timer.end("mcomp");
    }
    if(Func == "comp"){
        timer.start("comp");
        compBS(data, key, target_key, result_data, role, enc, eval, runtime);
        timer.end("comp");

    }
    if(Func == "mtag"){
        timer.start("mtag");
        mtagBS(key, target_key, result, role, enc, eval, runtime);
        timer.end("mtag");
    }
    if(Func == "tag"){
        timer.start("tag");
        tagBS(data, key, target_key, result_data, role, enc, eval, runtime);
        timer.end("tag");
    }
    if(Func == "subH"){
        SET_OR_DEFAULT_INT(cmd, alpha, sqrtToPowerOfTwo(N));
        SET_OR_DEFAULT_INT(cmd, threshold, 1<<8);
        timer.start("subH");
        subHBS(data, key, target_key, result_data, role, enc, eval, runtime, alpha, threshold);
        timer.end("subH");

        if(role == 0){
            std::ofstream stream(record_file, std::ios::app);
            stream << "alpha: " << alpha << ", threshold: " << threshold << std::endl;  
        }
    }

    // print the results.
    #ifndef MPI_APP
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        timer.print_total("milliseconds", stream);
        cmeter.print_total_per_party("MB", stream); 
    }
    #endif

    return;
}