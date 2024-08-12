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

void create_custom_mpi_type(MPI_Datatype* mpi_type, size_t size, MPI_Datatype base_type){
    int block_lengths[1] = {size};
    MPI_Aint displacements[1] = {0};
    MPI_Datatype types[1] = {base_type};
    MPI_Type_create_struct(1, block_lengths, displacements, types, mpi_type);
    MPI_Type_commit(mpi_type);
    return;
}


void generate_data(aby3::si64Matrix &data, size_t n, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    i64Matrix random_data(n, 1);
    for(int i = 0; i < n; ++i){
        random_data(i, 0) = 1;
    }
    data.resize(n, 1);
    large_data_encryption(pIdx, random_data, data, enc, runtime);
    return;
}

void generate_data(std::vector<aby3::si64Matrix> &data, size_t n, size_t m, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    aby3::si64Matrix full_data;
    generate_data(full_data, n*m, pIdx, enc, eval, runtime);

    data.resize(n);
    for(int i = 0; i < n; ++i){
        data[i].resize(m, 1);
        for(int j = 0; j < m; ++j){
            data[i].mShares[0](j, 0) = full_data.mShares[0](i*m + j, 0);
            data[i].mShares[1](j, 0) = full_data.mShares[1](i*m + j, 0);
        }
    }
    return;
}

void generate_data(aby3::sbMatrix &data, size_t n, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    i64Matrix random_data(n, 1);
    for(int i = 0; i < n; ++i){
        random_data(i, 0) = 1;
    }
    data.resize(n, bitsize);
    large_data_encryption(pIdx, random_data, data, enc, runtime);

    return;
}

void generate_data(std::vector<aby3::sbMatrix> &data, size_t n, size_t m, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    aby3::sbMatrix full_data;
    generate_data(full_data, n*m, bitsize, pIdx, enc, eval, runtime);

    data.resize(n);
    for(int i = 0; i < n; ++i){
        data[i].resize(m, bitsize);
        for(int j = 0; j < m; ++j){
            data[i].mShares[0](j, 0) = full_data.mShares[0](i*m + j, 0);
            data[i].mShares[1](j, 0) = full_data.mShares[1](i*m + j, 0);
        }
    }
    return;
}

void generate_data_parallel(std::vector<aby3::si64Matrix> &data, size_t n, size_t m, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    // all processors generate the data and combine to the rank0.
    // the data size is split by n.
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int task_size = n / total_tasks;
    if(rank == (total_tasks - 1)){
        task_size = n - task_size * (total_tasks - 1);
    }

    // each processor generates the data.
    aby3::si64Matrix local_data;
    generate_data(local_data, task_size*m, pIdx, enc, eval, runtime);

    int MPI_BLOCKSIZE = 8;
    MPI_Datatype mpi_type;
    create_custom_mpi_type(&mpi_type, MPI_BLOCKSIZE, MPI_INT64_T);

    // gather the data to rank0.
    // 1. rank0 compute the sizes and offsets.
    if(rank == 0){
        std::vector<int> sizes(total_tasks);
        std::vector<int> offsets(total_tasks);
        for(int i = 0; i < total_tasks-1; ++i){
            sizes[i] = task_size * m / MPI_BLOCKSIZE;
            offsets[i] = i * task_size * m / MPI_BLOCKSIZE;
        }
        sizes[total_tasks-1] = (n - task_size * (total_tasks - 1)) / MPI_BLOCKSIZE;
        offsets[total_tasks-1] = (total_tasks - 1) * task_size / MPI_BLOCKSIZE;
        
        aby3::si64Matrix full_data(n*m, 1);

        std::copy(local_data.mShares[0].data(), local_data.mShares[0].data() + local_data.mShares[0].size(), full_data.mShares[0].data());
        std::copy(local_data.mShares[1].data(), local_data.mShares[1].data() + local_data.mShares[1].size(), full_data.mShares[1].data());

        MPI_Gatherv(MPI_IN_PLACE, 0, mpi_type, full_data.mShares[0].data(), sizes.data(), offsets.data(), mpi_type, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, mpi_type, full_data.mShares[1].data(), sizes.data(), offsets.data(), mpi_type, 0, MPI_COMM_WORLD);

        debug_mpi(rank, pIdx, "finish gather data! data size: " + std::to_string(full_data.rows()));

        data.resize(n);
        for(int i = 0; i < n; ++i){
            data[i].resize(m, 1);
            for(int j = 0; j < m; ++j){
                data[i].mShares[0](j, 0) = full_data.mShares[0](i*m + j, 0);
                data[i].mShares[1](j, 0) = full_data.mShares[1](i*m + j, 0);
            }
        }
    }
    else{
        MPI_Gatherv(local_data.mShares[0].data(), task_size / MPI_BLOCKSIZE, mpi_type, nullptr, nullptr, nullptr, mpi_type, 0, MPI_COMM_WORLD);
        MPI_Gatherv(local_data.mShares[1].data(), task_size / MPI_BLOCKSIZE, mpi_type, nullptr, nullptr, nullptr, mpi_type, 0, MPI_COMM_WORLD);
    }

    MPI_Type_free(&mpi_type);

    return;
}

void generate_data_parallel(std::vector<aby3::sbMatrix> &data, size_t n, size_t m, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    // all processors generate the data and combine to the rank0.
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int task_size = n / total_tasks;
    if(rank == (total_tasks - 1)){
        task_size = n - task_size * (total_tasks - 1);
    }

    // each processor generates the data.
    aby3::sbMatrix local_data;
    generate_data(local_data, task_size*m, bitsize, pIdx, enc, eval, runtime);

    // gather the data to rank0.
    // 1. rank0 compute the sizes and offsets.
    if(rank == 0){
        std::vector<int> sizes(total_tasks);
        std::vector<int> offsets(total_tasks);
        for(int i = 0; i < total_tasks-1; ++i){
            sizes[i] = task_size * m;
            offsets[i] = i * task_size;
        }
        sizes[total_tasks-1] = (n - task_size * (total_tasks - 1));
        offsets[total_tasks-1] = (total_tasks - 1) * task_size;

        // test the data sizes.
        aby3::sbMatrix full_data(n*m, bitsize);

        std::copy(local_data.mShares[0].data(), local_data.mShares[0].data() + local_data.mShares[0].size(), full_data.mShares[0].data());
        std::copy(local_data.mShares[1].data(), local_data.mShares[1].data() + local_data.mShares[1].size(), full_data.mShares[1].data());

        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, full_data.mShares[0].data(), sizes.data(), offsets.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, full_data.mShares[1].data(), sizes.data(), offsets.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);

        data.resize(n);
        for(int i = 0; i < n; ++i){
            data[i].resize(m, bitsize);
            for(int j = 0; j < m; ++j){
                data[i].mShares[0](j, 0) = full_data.mShares[0](i*m + j, 0);
                data[i].mShares[1](j, 0) = full_data.mShares[1](i*m + j, 0);
            }
        }
    }
    else{
        MPI_Gatherv(local_data.mShares[0].data(), task_size, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(local_data.mShares[1].data(), task_size, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
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

void pta_binary_search_benchmark(oc::CLP& cmd){

    SETUP_PROCESS

    // get the corresponding parameters.
    SET_OR_DEFAULT_INT(cmd, N, 1<<10);
    SET_OR_DEFAULT_INT(cmd, optB, 1<<10);
    SET_OR_DEFAULT_STRING(cmd, Func, "mcomp");
    SET_OR_DEFAULT_STRING(cmd, record_file, "/root/aby3/record.txt");

    // get the data.
    std::vector<aby3::si64Matrix> data;
    std::vector<aby3::sbMatrix> key;

    // debug_mpi(rank, role, "begin data generation: data size: " + std::to_string(N));

    generate_data_parallel(data, N, 1, role, enc, eval, runtime);
    generate_data_parallel(key, N, 1, 32, role, enc, eval, runtime);

    // debug_mpi(rank, role, "finish data generation! data size: " + std::to_string(data.size()));

    boolIndex _target_key(N/2, role);
    aby3::sbMatrix target_key = _target_key.to_matrix(32);
    aby3::sbMatrix result;
    aby3::si64Matrix result_data;

    auto ptaBSTask = new ABY3MPITask<sb64, sb64, si64, si64, PtABS>(size, optB, role, enc, runtime, eval);
    auto ptaMBSTask = new ABY3MPIPairOnlyTask<sb64, sb64, sb64, sb64, PtAMBS>(size, optB, role, enc, runtime, eval);
    auto ptaMBSABTask = new ABY3MPIPairOnlyTask<sb64, sb64, std::pair<si64, sb64>, std::pair<si64, sb64>, PtAMBSAB>(size, optB, role, enc, runtime, eval);

    // benchmark the search function.
    Timer& timer = Timer::getInstance();
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    if(Func == "mcomp"){
        MPI_Barrier(MPI_COMM_WORLD);
        timer.start("mcomp");
        ptaMBS(key, target_key, result, ptaMBSTask);
        MPI_Barrier(MPI_COMM_WORLD);
        timer.end("mcomp");
    }
    if(Func == "comp"){
        MPI_Barrier(MPI_COMM_WORLD);
        timer.start("comp");
        ptaBS(data, key, target_key, result_data, ptaBSTask);
        MPI_Barrier(MPI_COMM_WORLD);
        timer.end("comp");
    }
    if(Func == "subH"){
        SET_OR_DEFAULT_INT(cmd, alpha, sqrtToPowerOfTwo(N));
        SET_OR_DEFAULT_INT(cmd, threshold, 1<<8);
        MPI_Barrier(MPI_COMM_WORLD);
        timer.start("subH");
        subHBS_with_PtA(data, key, target_key, result_data, role, enc, eval, runtime, alpha, threshold, ptaBSTask, ptaMBSABTask);
        MPI_Barrier(MPI_COMM_WORLD);
        timer.end("subH");

        if(role == 0 && rank == 0){
            std::ofstream stream(record_file, std::ios::app);
            stream << "alpha: " << alpha << ", threshold: " << threshold << std::endl;  
        }
    }

    // print the results.
    if(role == 0 && rank == 0){
        std::ofstream stream(record_file, std::ios::app);
        timer.print_total("milliseconds", stream);
        cmeter.print_total_per_party("MB", stream); 
    }

}