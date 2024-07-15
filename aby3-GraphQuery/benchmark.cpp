#include "benchmark.h"

using namespace oc;
using namespace aby3;

#define SET_OR_DEFAULT(cmd, key, default_value) \
    size_t key = default_value; \
    if(cmd.isSet(#key)){ \
        auto keys = cmd.getMany<size_t>(#key); \
        key = keys[0]; \
    } \
    else{ \
        debug_info(#key " is not set, using default value: " + std::to_string(key)); \
    }

#define CONFIG_INIT \
    int role = -1; \
    if (cmd.isSet("role")) { \
        auto keys = cmd.getMany<int>("role"); \
        role = keys[0]; \
    } \
    if (role == -1) { \
        throw std::runtime_error(LOCATION); \
    } \
    IOService ios; \
    Sh3Encryptor enc; \
    Sh3Evaluator eval; \
    Sh3Runtime runtime; \
    basic_setup((u64)role, ios, enc, eval, runtime); \
    aby3Info party_info(role, enc, eval, runtime); \
    Timer& timer = Timer::getInstance(); \
    timer.clear_records();


#define MPI_INIT \ 
    int role = -1; \
    if (cmd.isSet("role")) { \
        auto keys = cmd.getMany<int>("role"); \
        role = keys[0]; \
    } \
    if (role == -1) { \
        throw std::runtime_error(LOCATION); \
    } \
    IOService ios; \
    Sh3Encryptor enc; \
    Sh3Evaluator eval; \
    Sh3Runtime runtime; \
    int rank; \
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
    int total_tasks; \
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks); \
    multi_processor_setup((u64)role, rank, ios, enc, eval, runtime); \
    aby3Info party_info(role, enc, eval, runtime); \
    Timer& timer = Timer::getInstance(); \
    timer.clear_records();


size_t get_sending_bytes(aby3Info &party_info){
    size_t send_next = party_info.runtime->mComm.mNext.getTotalDataSent();
    size_t send_prev = party_info.runtime->mComm.mPrev.getTotalDataSent();
    return send_next + send_prev;
}

size_t get_receiving_bytes(aby3Info &party_info){
    size_t recv_next = party_info.runtime->mComm.mNext.getTotalDataRecv();
    size_t recv_prev = party_info.runtime->mComm.mPrev.getTotalDataRecv();
    return recv_next + recv_prev;
}

void clear_sending_reveiving_bytes(aby3Info &party_info){
    party_info.runtime->mComm.mNext.resetStats();
    party_info.runtime->mComm.mPrev.resetStats();
    return;
}

void communication_synchronize(aby3Info &party_info){

    CommunicationMeter& cmeter = CommunicationMeter::getInstance();
    // synchronize the communications to the 0-party.
    std::vector<uint64_t> communications_value;
    std::vector<std::string> key_set;
    for(const auto& pair : cmeter.communications){
        std::string key = pair.first;
        uint64_t comm = pair.second;
        communications_value.push_back(comm);
        key_set.push_back(key);
    }

    std::vector<uint64_t> comm_next(communications_value.size());
    std::vector<uint64_t> comm_prev(communications_value.size());

    if(party_info.pIdx == 0){
        auto tmp_next = party_info.runtime->mComm.mNext.asyncRecv(comm_next.data(), comm_next.size());
        auto tmp_prev = party_info.runtime->mComm.mPrev.asyncRecv(comm_prev.data(), comm_prev.size());
        tmp_next.get();
        tmp_prev.get();
    }
    if(party_info.pIdx == 1){
        auto tmp = party_info.runtime->mComm.mPrev.asyncSendFuture<uint64_t>(communications_value.data(), communications_value.size());
        tmp.get();
    }
    if(party_info.pIdx == 2){
        auto tmp = party_info.runtime->mComm.mNext.asyncSendFuture<uint64_t>(communications_value.data(), communications_value.size());
        tmp.get();
    }

    if(party_info.pIdx == 0){
        for(size_t i=0; i<key_set.size(); i++){
            std::string key = key_set[i];
            cmeter.totalCommunications[key] = {communications_value[i], comm_next[i], comm_prev[i]};
        }
    }

    return;
}

int privGraph_performance_profiling(oc::CLP& cmd){

    #ifdef MPI_APP
    MPI_INIT
    #else
    CONFIG_INIT
    if(role == 0) debug_info("in this function!!!!");
    #endif

    // get the graph file parameters.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/baseline/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/privGraph/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }

    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_2dpartition.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    if (role == 0) debug_info("Profiling the privGraph Query Engine...");

    // get the oram configurations, todo - define the parameters.
    size_t noram_stash_size = 1 << 4;
    size_t noram_pack_size = 1 << 2;
    size_t eoram_stash_size = 1 << 8;
    size_t eoram_pack_size = 1 << 4;

    if(cmd.isSet("noram_stash_size")){
        auto keys = cmd.getMany<size_t>("noram_stash_size");
        noram_stash_size = keys[0];
    }
    else{
        debug_info("noram_stash_size is not set, using default value: " + std::to_string(noram_stash_size));
    }
    if(cmd.isSet("noram_pack_size")){
        auto keys = cmd.getMany<size_t>("noram_pack_size");
        noram_pack_size = keys[0];
    }
    else{
        debug_info("noram_pack_size is not set, using default value: " + std::to_string(noram_pack_size));
    }
    if(cmd.isSet("eoram_stash_size")){
        auto keys = cmd.getMany<size_t>("eoram_stash_size");
        eoram_stash_size = keys[0];
    }
    else{
        debug_info("eoram_stash_size is not set, using default value: " + std::to_string(eoram_stash_size));
    }
    if(cmd.isSet("eoram_pack_size")){
        auto keys = cmd.getMany<size_t>("eoram_pack_size");
        eoram_pack_size = keys[0];
    }
    else{
        debug_info("eoram_pack_size is not set, using default value: " + std::to_string(eoram_pack_size));
    }

    if(role == 0) debug_info("Getting oram configs success");

    if(role == 0) debug_info("Environment setup success");

    // load the graph.
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));

    if(role == 0) debug_info("sending_bytes = " + std::to_string(get_sending_bytes(party_info)) + " | reveiving_bytes = " + std::to_string(get_receiving_bytes(party_info)));

    timer.start("GraphLoad");

    GraphQueryEngine secGraphEngine(party_info, meta_file, graph_data_file);

    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);


    if(role == 0) debug_info("Graph loaded successfully!");

    // oram construction.
    if(role == 0) debug_info("Eoram construction...");

    cmeter.start("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("EdgeOramInit_recv", get_receiving_bytes(party_info));
    timer.start("EdgeOramInit");

    secGraphEngine.edge_block_oram_initialization(eoram_stash_size, eoram_pack_size);

    timer.end("EdgeOramInit");
    cmeter.end("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("EdgeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Eoram init success\nNoram construction...");

    size_t b = secGraphEngine.graph->b;
    size_t b2 = secGraphEngine.graph->edge_list_size;
    size_t v = secGraphEngine.graph->v;

    eoram_stash_size = secGraphEngine.edge_block_oram->S;

    // basic graph query process.
    size_t ps = 0, pe = v-1;
    boolIndex snode = boolIndex(secGraphEngine.get_block_index(ps), role);
    boolIndex enode = boolIndex(secGraphEngine.get_block_index(pe), role);
    boolIndex edge_log_idx = boolIndex(secGraphEngine.get_edge_block_index(ps, pe), role);
    boolIndex snode_log_idx = boolIndex(secGraphEngine.get_block_index(ps), role);

    // 1) edge existence query.
    cmeter.start("EdgeExistQuery_send", get_sending_bytes(party_info));
    cmeter.start("EdgeExistQuery_recv", get_receiving_bytes(party_info));

    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    timer.start("EdgeExistQuery");

    for(int i=0; i<1; i++){
        boolShare flag = edge_existance(snode, enode, edge_log_idx, secGraphEngine);
    }
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("EdgeExistQuery");
    cmeter.end("EdgeExistQuery_send", get_sending_bytes(party_info));
    cmeter.end("EdgeExistQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Edge existence query success");

    delete secGraphEngine.edge_block_oram;

    cmeter.start("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("NodeOramInit_recv", get_receiving_bytes(party_info));
    timer.start("NodeOramInit");

    secGraphEngine.node_edges_oram_initialization(noram_stash_size, noram_pack_size); 
    noram_stash_size = secGraphEngine.node_edges_oram->S;

    timer.end("NodeOramInit");
    cmeter.end("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("NodeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Noram init success");

    // 2) outting edges count query.
    cmeter.start("OuttingEdgesCountQuery_send", get_sending_bytes(party_info));
    cmeter.start("OuttingEdgesCountQuery_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("OuttingEdgesCountQuery");
    for(int i=0; i<1; i++){
        aby3::si64Matrix out_edges = outting_edge_count(snode, snode_log_idx, secGraphEngine);
    }
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("OuttingEdgesCountQuery");
    cmeter.end("OuttingEdgesCountQuery_send", get_sending_bytes(party_info));
    cmeter.end("OuttingEdgesCountQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Outting edges count query success");

    if(role == 0) debug_info("Neighbors get query begin!");

    cmeter.start("NeighborsGetQuery_send", get_sending_bytes(party_info));
    cmeter.start("NeighborsGetQuery_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("NeighborsGetQuery");
    for(size_t i=0; i<1; i++){
        aby3::sbMatrix neighbors = outting_neighbors_sorted(snode, snode_log_idx, secGraphEngine);
    }
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("NeighborsGetQuery");
    cmeter.end("NeighborsGetQuery_send", get_sending_bytes(party_info));
    cmeter.end("NeighborsGetQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Neighbors get query success");

    communication_synchronize(party_info);

    #ifdef MPI_APP
    if((role == 0) && (rank == 0)){
    #else
    if(role == 0){
    #endif
        std::ofstream stream(record_file, std::ios::app);
        secGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int adj_performance_profiling(oc::CLP& cmd){

    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    // get the graph file parameters.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/baseline/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/adjmat/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }

    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_edge_list.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    if (role == 0) debug_info("Profiling the Mat Query Engine...");

    // get the oram configurations, todo - define the parameters.
    size_t noram_stash_size = 1 << 4;
    size_t noram_pack_size = 1 << 2;
    size_t eoram_stash_size = 1 << 8;
    size_t eoram_pack_size = 1 << 4;

    if(cmd.isSet("noram_stash_size")){
        auto keys = cmd.getMany<size_t>("noram_stash_size");
        noram_stash_size = keys[0];
    }
    else{
        debug_info("noram_stash_size is not set, using default value: " + std::to_string(noram_stash_size));
    }
    if(cmd.isSet("noram_pack_size")){
        auto keys = cmd.getMany<size_t>("noram_pack_size");
        noram_pack_size = keys[0];
    }
    else{
        debug_info("noram_pack_size is not set, using default value: " + std::to_string(noram_pack_size));
    }
    if(cmd.isSet("eoram_stash_size")){
        auto keys = cmd.getMany<size_t>("eoram_stash_size");
        eoram_stash_size = keys[0];
    }
    else{
        debug_info("eoram_stash_size is not set, using default value: " + std::to_string(eoram_stash_size));
    }
    if(cmd.isSet("eoram_pack_size")){
        auto keys = cmd.getMany<size_t>("eoram_pack_size");
        eoram_pack_size = keys[0];
    }
    else{
        debug_info("eoram_pack_size is not set, using default value: " + std::to_string(eoram_pack_size));
    }

    if(role == 0) debug_info("Getting oram configs success");

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    aby3Info party_info(role, enc, eval, runtime);

    if(role == 0) debug_info("Environment setup success");

    // get the timer.
    Timer& timer = Timer::getInstance();
    timer.clear_records();
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // graph loading.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");
    AdjGraphQueryEngine adjGraphEngine(party_info, meta_file, graph_data_file);
    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");

    // oram construction.
    if(role == 0) debug_info("Eoram construction...");
    cmeter.start("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("EdgeOramInit_recv", get_receiving_bytes(party_info));
    timer.start("EdgeOramInit");
    adjGraphEngine.edge_oram_initialization(eoram_stash_size, eoram_pack_size);
    timer.end("EdgeOramInit");
    cmeter.end("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("EdgeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Eoram init success\nNoram construction...");

    cmeter.start("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("NodeOramInit_recv", get_receiving_bytes(party_info));
    timer.start("NodeOramInit");
    adjGraphEngine.node_oram_initialization(noram_stash_size, noram_pack_size);
    timer.end("NodeOramInit");
    cmeter.end("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("NodeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    eoram_stash_size = adjGraphEngine.edge_oram->S;
    noram_stash_size = adjGraphEngine.node_oram->S;

    if(role == 0) debug_info("Noram init success");

    // 1) edge existence query.
    size_t v = adjGraphEngine.graph->v;
    size_t ps = 0, pe = v-1;
    boolIndex snode = boolIndex(ps, role);
    boolIndex enode = boolIndex(pe, role);

    cmeter.start("EdgeExistQuery_send", get_sending_bytes(party_info));
    cmeter.start("EdgeExistQuery_recv", get_receiving_bytes(party_info));
    timer.start("EdgeExistQuery");
    for(int i=0; i<eoram_stash_size; i++){
        boolShare flag = edge_existance(snode, enode, adjGraphEngine);
    }
    timer.end("EdgeExistQuery");
    cmeter.end("EdgeExistQuery_send", get_sending_bytes(party_info));
    cmeter.end("EdgeExistQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Edge existence query success");

    // 2) outting edges count query.
    cmeter.start("OuttingEdgesCountQuery_send", get_sending_bytes(party_info));
    cmeter.start("OuttingEdgesCountQuery_recv", get_receiving_bytes(party_info));
    timer.start("OuttingEdgesCountQuery");
    for(int i=0; i<noram_stash_size; i++){
        aby3::sbMatrix out_edges = outting_edge_count(snode, adjGraphEngine);
    }
    timer.end("OuttingEdgesCountQuery");
    cmeter.end("OuttingEdgesCountQuery_send", get_sending_bytes(party_info));
    cmeter.end("OuttingEdgesCountQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    // 3) neighbors get query.
    adjGraphEngine.node_oram_initialization(noram_stash_size, noram_pack_size);
    cmeter.start("NeighborsGetQuery_send", get_sending_bytes(party_info));
    cmeter.start("NeighborsGetQuery_recv", get_receiving_bytes(party_info));
    timer.start("NeighborsGetQuery");
    for(size_t i=0; i<noram_stash_size; i++){
        aby3::sbMatrix neighbors = outting_neighbors(snode, adjGraphEngine);
    }
    timer.end("NeighborsGetQuery");
    cmeter.end("NeighborsGetQuery_send", get_sending_bytes(party_info));
    cmeter.end("NeighborsGetQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    // print the timer records.
    communication_synchronize(party_info);
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        adjGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int list_performance_profiling(oc::CLP& cmd){

    // get the configs.
    // int role = -1;
    // if (cmd.isSet("role")) {
    //     auto keys = cmd.getMany<int>("role");
    //     role = keys[0];
    // }
    // if (role == -1) {
    //     throw std::runtime_error(LOCATION);
    // }

    #ifdef MPI_APP
    MPI_INIT
    #else
    CONFIG_INIT
    if(role == 0) debug_info("in this function!!!!");
    #endif

    // get the graph file parameters.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/baseline/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/edgelist/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }
    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_edge_list.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    if (role == 0) debug_info("Profiling the Edgelist Query Engine...");

    // // setup communications.
    // IOService ios;
    // Sh3Encryptor enc;
    // Sh3Evaluator eval;
    // Sh3Runtime runtime;
    // basic_setup((u64)role, ios, enc, eval, runtime);
    // aby3Info party_info(role, enc, eval, runtime);

    if(role == 0) debug_info("Environment setup success");

    // // get the timer.
    // Timer& timer = Timer::getInstance();
    // timer.clear_records();
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // graph loading.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");
    plainGraphList plainGraph(meta_file, graph_data_file);
    ListGraphQueryEngine listGraphEngine(party_info, plainGraph);
    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");

    // 1) edge existence query.
    size_t v = plainGraph.v;
    size_t ps = 0, pe = v-1;
    boolIndex snode = boolIndex(ps, role);
    boolIndex enode = boolIndex(pe, role);

    cmeter.start("EdgeExistQuery_send", get_sending_bytes(party_info));
    cmeter.start("EdgeExistQuery_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("EdgeExistQuery");
    boolShare flag = edge_existance(snode, enode, listGraphEngine);
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("EdgeExistQuery");
    cmeter.end("EdgeExistQuery_send", get_sending_bytes(party_info));
    cmeter.end("EdgeExistQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Edge existence query success");

    // 2) outting edges count query.
    cmeter.start("OuttingEdgesCountQuery_send", get_sending_bytes(party_info));
    cmeter.start("OuttingEdgesCountQuery_recv", get_receiving_bytes(party_info)); 
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif  
    timer.start("OuttingEdgesCountQuery");
    // aby3::sbMatrix out_edges = outting_edge_count(snode, listGraphEngine);
    aby3::si64Matrix out_edges = outting_edge_count_arith(snode, listGraphEngine);
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("OuttingEdgesCountQuery");
    cmeter.end("OuttingEdgesCountQuery_send", get_sending_bytes(party_info));
    cmeter.end("OuttingEdgesCountQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Outting edges count query success");

    // // 3) neighbors get query.
    // timer.start("NeighborsGetQuery");
    // aby3::sbMatrix neighbors = outting_neighbors(snode, listGraphEngine);
    // timer.end("NeighborsGetQuery");

    if(role == 0) debug_info("Neighbors get query success");

    // 4) sorted neighbors get query.
    plainGraph.list_sort();
    listGraphEngine.rebuild(party_info, plainGraph);

    cmeter.start("NeighborsGetQuery_send", get_sending_bytes(party_info));
    cmeter.start("NeighborsGetQuery_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("NeighborsGetQuery");
    aby3::sbMatrix neighbors_sorted = outting_neighbors_sorted(snode, listGraphEngine);
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("NeighborsGetQuery");
    cmeter.end("NeighborsGetQuery_send", get_sending_bytes(party_info));
    cmeter.end("NeighborsGetQuery_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Sorted Neighbors get query success");

    communication_synchronize(party_info);
    // print the timer records.
    #ifdef MPI_APP
    if((role == 0) && (rank == 0)){
    #else
    if(role == 0){
    #endif
        std::ofstream stream(record_file, std::ios::app);
        listGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }


    return 0;
}

int privGraph_integration_profiling(oc::CLP& cmd){
    
    // get the configs.
    CONFIG_INIT
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // get the graph file parameters.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/multiparty/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record_offline/privGraph/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    std::string meta_file = graph_data_folder + file_prefix + "_meta_multiparty.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    file_prefix = graph_data_folder + file_prefix;

    SET_OR_DEFAULT(cmd, noram_stash_size, 1 << 4);
    SET_OR_DEFAULT(cmd, noram_pack_size, 1 << 2);
    SET_OR_DEFAULT(cmd, eoram_stash_size, 1 << 8);
    SET_OR_DEFAULT(cmd, eoram_pack_size, 1 << 4);

    // std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    if (role == 0) debug_info("Profiling the privGraph Query Engine...");

    // graph integration.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");
    GraphQueryEngine secGraphEngine(party_info, meta_file, file_prefix, true);
    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");

    // oram construction.
    if(role == 0) debug_info("Eoram construction...");
    cmeter.start("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("EdgeOramInit_recv", get_receiving_bytes(party_info));
    timer.start("EdgeOramInit");
    secGraphEngine.edge_block_oram_initialization(eoram_stash_size, eoram_pack_size);
    timer.end("EdgeOramInit");
    cmeter.end("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("EdgeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Eoram init success\nNoram construction...");

    cmeter.start("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("NodeOramInit_recv", get_receiving_bytes(party_info));
    timer.start("NodeOramInit");
    secGraphEngine.node_edges_oram_initialization(noram_stash_size, noram_pack_size);
    timer.end("NodeOramInit");
    cmeter.end("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("NodeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Noram init success");

    communication_synchronize(party_info);
    // print the timer records.
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        debug_info("Printing configs... " + record_file);
        secGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
        stream.close();
    }

    return 0;
}

int adj_integration_profiling(oc::CLP& cmd){

    // get the configs.
    CONFIG_INIT
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // get the graph file parameters.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/multiparty/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record_offline/adjmat/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta_multiparty.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";
    file_prefix = graph_data_folder + file_prefix;

    SET_OR_DEFAULT(cmd, noram_stash_size, 1 << 4);
    SET_OR_DEFAULT(cmd, noram_pack_size, 1 << 2);
    SET_OR_DEFAULT(cmd, eoram_stash_size, 1 << 8);
    SET_OR_DEFAULT(cmd, eoram_pack_size, 1 << 4);

    if (role == 0) debug_info("Profiling the Mat Query Engine...");

    // graph integration.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");
    AdjGraphQueryEngine adjGraphEngine(party_info, meta_file, file_prefix, true);
    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");

    // oram construction.
    if(role == 0) debug_info("Eoram construction...");
    cmeter.start("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("EdgeOramInit_recv", get_receiving_bytes(party_info));
    timer.start("EdgeOramInit");
    adjGraphEngine.edge_oram_initialization(eoram_stash_size, eoram_pack_size);
    timer.end("EdgeOramInit");
    cmeter.end("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("EdgeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Eoram init success\nNoram construction...");

    cmeter.start("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("NodeOramInit_recv", get_receiving_bytes(party_info));
    timer.start("NodeOramInit");
    adjGraphEngine.node_oram_initialization(noram_stash_size, noram_pack_size);
    timer.end("NodeOramInit");
    cmeter.end("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("NodeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Noram init success");

    communication_synchronize(party_info);

    // print the timer records.
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        adjGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int list_integration_profiling(oc::CLP& cmd){

    // get the configs.
    CONFIG_INIT
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // get the graph file parameters.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/multiparty/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record_offline/edgelist/";

    int record_counter = -1;
    
    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta_multiparty.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";
    file_prefix = graph_data_folder + file_prefix;

    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");
    ListGraphQueryEngine listGraphEngine(party_info, meta_file, file_prefix, true);
    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");
    communication_synchronize(party_info);

    // print the timer records.
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        listGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int cycle_detection_profiling(oc::CLP& cmd){

    #ifdef MPI_APP
    MPI_INIT
    #else
    CONFIG_INIT
    if(role == 0) debug_info("in this function!!!!");
    #endif

    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    if(role == 0) debug_info("in this function!!!");

    // init the EORAM.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/adv_application/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/cycle_detection/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }
    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }

    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_2dpartition.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    SET_OR_DEFAULT(cmd, eoram_stash_size, 1 << 8);
    SET_OR_DEFAULT(cmd, eoram_pack_size, 1 << 4);

    if (role == 0) debug_info("Profiling the cycle detection...");

    if(role == 0){
        debug_info("meta_file = " + meta_file);
        debug_info("graph_data_file = " + graph_data_file);
        debug_info("record_file = " + record_file);
    }

    // load the graph and construct the eoram.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));

    if(role == 0) debug_info("sending_bytes = " + std::to_string(get_sending_bytes(party_info)) + " | reveiving_bytes = " + std::to_string(get_receiving_bytes(party_info)));

    timer.start("GraphLoad");

    GraphQueryEngine secGraphEngine(party_info, meta_file, graph_data_file);

    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");

    // oram construction.
    cmeter.start("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("EdgeOramInit_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("EdgeOramInit");

    secGraphEngine.edge_block_oram_initialization(eoram_stash_size, eoram_pack_size);
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("EdgeOramInit");
    cmeter.end("EdgeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("EdgeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Eoram init success");

    // cycle detection.
    // set the six edges.
    size_t node1 = 1, node2 = 2, node3 = 3;
    boolIndex node1_ind = boolIndex(secGraphEngine.get_block_index(node1), role);
    boolIndex node2_ind = boolIndex(secGraphEngine.get_block_index(node2), role);
    boolIndex node3_ind = boolIndex(secGraphEngine.get_block_index(node3), role);

    std::vector<std::tuple<boolIndex, boolIndex, boolIndex>> target_edges(6);
    target_edges[0] = {boolIndex(secGraphEngine.get_edge_block_index(node1, node2), role), node1_ind, node2_ind};
    target_edges[1] = {boolIndex(secGraphEngine.get_edge_block_index(node2, node3), role), node2_ind, node3_ind};
    target_edges[2] = {boolIndex(secGraphEngine.get_edge_block_index(node3, node1), role), node3_ind, node1_ind};
    target_edges[3] = {boolIndex(secGraphEngine.get_edge_block_index(node2, node1), role), node2_ind, node1_ind};
    target_edges[4] = {boolIndex(secGraphEngine.get_edge_block_index(node3, node2), role), node3_ind, node2_ind};
    target_edges[5] = {boolIndex(secGraphEngine.get_edge_block_index(node1, node3), role), node1_ind, node3_ind};

    std::vector<boolShare> edge_flags;

    if(role == 0) debug_info("Cycle detection...");

    cmeter.start("CycleDetection_send", get_sending_bytes(party_info));
    cmeter.start("CycleDetection_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("CycleDetection");

    for(auto& edge_info : target_edges){
        auto& edge_index = std::get<0>(edge_info);
        auto& node1_index = std::get<1>(edge_info);
        auto& node2_index = std::get<2>(edge_info);

        boolShare flag = edge_existance(node1_index, node2_index, edge_index, secGraphEngine);
        edge_flags.push_back(flag);
    }
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("CycleDetection");
    cmeter.end("CycleDetection_send", get_sending_bytes(party_info));
    cmeter.end("CycleDetection_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Cycle detection success");

    communication_synchronize(party_info);

    // print the timer records.
    #ifdef MPI_APP
    if((role == 0) && (rank == 0)){
    #else
    if(role == 0){
    #endif
        std::ofstream stream(record_file, std::ios::app);
        secGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int twohop_neighbor_profiling(oc::CLP& cmd){

    #ifdef MPI_APP
    MPI_INIT
    #else
    CONFIG_INIT
    if(role == 0) debug_info("in this function!!!!");
    #endif
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // init the NORAM.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/adv_application/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/two_hop_neighbor/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }
    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }

    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_2dpartition.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    SET_OR_DEFAULT(cmd, noram_stash_size, 1 << 8);
    SET_OR_DEFAULT(cmd, noram_pack_size, 1 << 4);

    if(role == 0){
        debug_info("meta_file = " + meta_file);
        debug_info("graph_data_file = " + graph_data_file);
        debug_info("record_file = " + record_file);
    }

    if (role == 0) debug_info("Profiling the two-hop neighbor query...");

    // load the graph and construct the noram.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");

    GraphQueryEngine secGraphEngine(party_info, meta_file, graph_data_file);

    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");

    // oram construction.
    cmeter.start("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("NodeOramInit_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("NodeOramInit");

    secGraphEngine.node_edges_oram_initialization(noram_stash_size, noram_pack_size);
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("NodeOramInit");
    cmeter.end("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("NodeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Noram init success");

    // // 2-htop neighbor detection.
    // size_t node1 = 1;
    // boolIndex snode = boolIndex(node1, role);
    // boolIndex node1_ind = boolIndex(secGraphEngine.get_block_index(node1), role);

    // cmeter.start("twohop_neighbor_send", get_sending_bytes(party_info));
    // cmeter.start("twohop_neighbor_recv", get_receiving_bytes(party_info));
    // #ifdef MPI_APP
    // MPI_Barrier(MPI_COMM_WORLD);
    // #endif
    // timer.start("twohop_neighbor");    

    // aby3::sbMatrix direct_neighbors = outting_neighbors_sorted(snode, node1_ind, secGraphEngine);

    // aby3::i64Matrix direct_neighbors_plain(direct_neighbors.rows(), 1);
    // enc.revealAll(runtime, direct_neighbors, direct_neighbors_plain).get();

    // size_t neighbor_upper_bound = 10;
    // if(file_prefix.find("slashdot") == 0){
    //     neighbor_upper_bound = 7;
    // }
    // if(file_prefix.find("dblp") == 0){
    //     neighbor_upper_bound = 4;
    // }
    // if(file_prefix.find("twitter") == 0){
    //     neighbor_upper_bound = 5;
    // }

    // // size_t direct_neighbor_number = 0;
    // std::vector<size_t> direct_neighbors_list;
    // // for(size_t i=0; i<direct_neighbors_plain.rows(); i++){
    // //     if(direct_neighbors_plain(i, 0) != 0){
    // //         direct_neighbor_number++;
    // //         direct_neighbors_list.push_back(i);
    // //     }
    // // }
    // while(direct_neighbors_list.size() < neighbor_upper_bound){
    //     direct_neighbors_list.push_back(10);
    // }

    // // 2-hop neighbor detection.
    // std::vector<size_t> two_hop_neighbors_list;
    // for(auto& neighbor : direct_neighbors_list){
    //     boolIndex neighbor_id_ = boolIndex(neighbor, role);
    //     boolIndex neighbor_ind_ = boolIndex(secGraphEngine.get_edge_block_index(node1, neighbor), role);
    //     aby3::sbMatrix two_hop_neighbors = outting_neighbors_sorted(neighbor_id_, neighbor_ind_, secGraphEngine);
    // }
    // #ifdef MPI_APP
    // MPI_Barrier(MPI_COMM_WORLD);
    // #endif
    // timer.end("twohop_neighbor");
    // cmeter.end("twohop_neighbor_send", get_sending_bytes(party_info));
    // cmeter.end("twohop_neighbor_recv", get_receiving_bytes(party_info));
    // clear_sending_reveiving_bytes(party_info);

    // // print the timer records.

    communication_synchronize(party_info);  

    #ifdef MPI_APP
    if((role == 0) && (rank == 0)){
    #else
    if(role == 0){
    #endif
        std::ofstream stream(record_file, std::ios::app);
        // debug_info("Direct neighbors of node " + std::to_string(direct_neighbors_list.size()), stream);
        secGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int neighbor_statistics_profiling(oc::CLP& cmd){

    #ifdef MPI_APP
    MPI_INIT
    #else
    CONFIG_INIT
    if(role == 0) debug_info("in this function!!!!");
    #endif
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // set the graph.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/adv_application/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/statistics_neighbor/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }
    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }

    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_2dpartition.txt";
    std::string graph_property_file = graph_data_folder + file_prefix + "_2dproperty.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    SET_OR_DEFAULT(cmd, noram_stash_size, 1 << 8);
    SET_OR_DEFAULT(cmd, noram_pack_size, 1 << 4);

    if(role == 0){
        debug_info("meta_file = " + meta_file);
        debug_info("graph_data_file = " + graph_data_file);
        debug_info("property_graph_file = " + graph_property_file);
        debug_info("record_file = " + record_file);
    }

    if (role == 0) debug_info("Profiling the neighbor statis query...");

    // get the property graph.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");

    PropertyGraph2d property_graph(meta_file, graph_data_file, party_info);
    property_graph.add_property(graph_property_file, party_info);
    property_graph.get_node_edges_property();


    // construct the query_engine.
    GraphQueryEngine secGraphEngine;
    secGraphEngine.party_info = &party_info;
    secGraphEngine.set_graph(property_graph);

    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");


    // oram construction.
    cmeter.start("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.start("NodeOramInit_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("NodeOramInit");

    if(role == 0) debug_info("Noram construction...");
    secGraphEngine.node_edges_property_oram_initialization(noram_stash_size, noram_pack_size);
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("NodeOramInit");
    cmeter.end("NodeOramInit_send", get_sending_bytes(party_info));
    cmeter.end("NodeOramInit_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Noram init success");

    // run the statistic_func.
    size_t node1 = 0;
    boolIndex snode = boolIndex(node1, role);
    boolIndex node1_ind = boolIndex(secGraphEngine.get_block_index(node1), role);
    size_t upper_bound = 10;
    boolIndex sec_upper_bound = boolIndex(upper_bound, role);

    cmeter.start("statistic_send", get_sending_bytes(party_info));
    cmeter.start("statistic_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("statistic");

    aby3::si64Matrix statistics = outting_edge_range_statistics(snode, node1_ind, sec_upper_bound, secGraphEngine); 
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("statistic");
    cmeter.end("statistic_send", get_sending_bytes(party_info));
    cmeter.end("statistic_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Statistics success");

    // print the timer records.
    communication_synchronize(party_info);
    #ifdef MPI_APP
    if((role == 0) && (rank == 0)){
    #else
    if(role == 0){
    #endif
        std::ofstream stream(record_file, std::ios::app);
        secGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int cycle_detection_profiling_edgelist(oc::CLP& cmd){

    #ifdef MPI_APP
    MPI_INIT
    #else
    CONFIG_INIT
    if(role == 0) debug_info("in this function!!!!");
    #endif
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // load graph.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/adv_application/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/cycle_detection/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }
    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_edge_list.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    // graph loading.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");
    plainGraphList plainGraph(meta_file, graph_data_file);
    ListGraphQueryEngine listGraphEngine(party_info, plainGraph);
    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    // cycle detection.
    size_t node1 = 1, node2 = 2, node3 = 3;
    boolIndex node1_ind = boolIndex(node1, role);   
    boolIndex node2_ind = boolIndex(node2, role);
    boolIndex node3_ind = boolIndex(node3, role);

    std::vector<std::pair<boolIndex, boolIndex>> target_edges(6);
    target_edges[0] = {node1_ind, node2_ind};
    target_edges[1] = {node2_ind, node3_ind};
    target_edges[2] = {node3_ind, node1_ind};
    target_edges[3] = {node2_ind, node1_ind};
    target_edges[4] = {node3_ind, node2_ind};
    target_edges[5] = {node1_ind, node3_ind};

    std::vector<boolShare> edge_flags;

    if(role == 0) debug_info("Cycle detection...");

    cmeter.start("CycleDetection_send", get_sending_bytes(party_info));
    cmeter.start("CycleDetection_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("CycleDetection");

    for(auto& edge_info : target_edges){
        auto& node1_index = edge_info.first;
        auto& node2_index = edge_info.second;

        boolShare flag = edge_existance(node1_index, node2_index, listGraphEngine);
        edge_flags.push_back(flag);
    }
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("CycleDetection");
    cmeter.end("CycleDetection_send", get_sending_bytes(party_info));
    cmeter.end("CycleDetection_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Cycle detection success");

    communication_synchronize(party_info);

    // print the timer records.
    #ifdef MPI_APP
    if((role == 0) && (rank == 0)){
    #else
    if(role == 0){
    #endif
        std::ofstream stream(record_file, std::ios::app);
        stream << "===== Edge list ====" <<std::endl;
        listGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int twohop_neighbor_profiling_edgelist(oc::CLP& cmd){

    #ifdef MPI_APP
    MPI_INIT
    #else
    CONFIG_INIT
    if(role == 0) debug_info("in this function!!!!");
    #endif
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // load graph.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/adv_application/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/cycle_detection/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }
    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_edge_list.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    // graph loading.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");
    plainGraphList plainGraph(meta_file, graph_data_file);
    ListGraphQueryEngine listGraphEngine(party_info, plainGraph);
    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Graph loaded successfully!");

    // two-hop neighbors.
    size_t node1 = 0;
    boolIndex snode = boolIndex(node1, role);
    
    cmeter.start("twohop_neighbor_send", get_sending_bytes(party_info));
    cmeter.start("twohop_neighbor_recv", get_receiving_bytes(party_info));
    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.start("twohop_neighbor");    

    aby3::sbMatrix direct_neighbors = outting_neighbors_sorted(snode, listGraphEngine);

    if(role == 0) debug_info("get the first direct neighbors! size = " + std::to_string(direct_neighbors.rows()) + " x " + std::to_string(direct_neighbors.i64Size()));

    aby3::i64Matrix direct_neighbors_plain(direct_neighbors.i64Size(), 1);
    // enc.revealAll(runtime, direct_neighbors, direct_neighbors_plain).get();
    large_data_decryption(role, direct_neighbors, direct_neighbors_plain, enc, runtime);


    if(role == 0) debug_info("revealed the first direct neighbors!");

    size_t neighbor_upper_bound = 10;
    if(file_prefix.find("slashdot") == 0){
        neighbor_upper_bound = 7;
    }
    if(file_prefix.find("dblp") == 0){
        neighbor_upper_bound = 4;
    }
    if(file_prefix.find("twitter") == 0){
        neighbor_upper_bound = 5;
    }


    if(role == 0) debug_info("in the second stage!");
    size_t direct_neighbor_number = 0;
    std::vector<size_t> direct_neighbors_list;
    for(size_t i=0; i<direct_neighbors_plain.rows(); i++){
        if(direct_neighbor_number >= neighbor_upper_bound){
            break;
        }
        if(direct_neighbors_plain(i, 0) != 0){
            direct_neighbor_number++;
            direct_neighbors_list.push_back(i);
        }
    }
    while(direct_neighbors_list.size() < neighbor_upper_bound){
        direct_neighbors_list.push_back(0);
    }

    // 2-hop neighbor detection.
    size_t neighbor_count = direct_neighbors_list.size();
    std::vector<size_t> two_hop_neighbors_list;
    for(auto& neighbor : direct_neighbors_list){
        boolIndex neighbor_id_ = boolIndex(neighbor, role);
        aby3::sbMatrix two_hop_neighbors = outting_neighbors_sorted(neighbor_id_, listGraphEngine); 
        if(role == 0) debug_info("getting neighbors of " + std::to_string(neighbor_count) + " / " + std::to_string(direct_neighbor_number));
    }

    #ifdef MPI_APP
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    timer.end("twohop_neighbor");
    cmeter.end("twohop_neighbor_send", get_sending_bytes(party_info));
    cmeter.end("twohop_neighbor_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    // print the timer records.
    communication_synchronize(party_info);
    #ifdef MPI_APP
    if((role == 0) && (rank == 0)){
    #else
    if(role == 0){
    #endif
        std::ofstream stream(record_file, std::ios::app);
        stream << "===== Edge list ====" <<std::endl;
        debug_info("Direct neighbors of node " + std::to_string(direct_neighbor_number), stream);
        listGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}

int neighbor_statistics_profiling_edgelist(oc::CLP& cmd){

    #ifdef MPI_APP
    MPI_INIT
    #else
    CONFIG_INIT
    if(role == 0) debug_info("in this function!!!!");
    #endif
    CommunicationMeter& cmeter = CommunicationMeter::getInstance();

    // load graph.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/adv_application/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/cycle_detection/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }

    if(cmd.isSet("data_folder")){
        auto keys = cmd.getMany<std::string>("data_folder");
        graph_data_folder = keys[0];
        if(role == 0) debug_info("data_folder: " + graph_data_folder);
    }
    if(cmd.isSet("record_folder")){
        auto keys = cmd.getMany<std::string>("record_folder");
        record_folder = keys[0];
        if(role == 0) debug_info("record_folder: " + record_folder);
    }

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_edge_list.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";
    std::string graph_property_file = graph_data_folder + file_prefix + "_edge_list_property.txt";

    // graph loading.
    cmeter.start("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.start("GraphLoad_recv", get_receiving_bytes(party_info));
    timer.start("GraphLoad");
    plainGraphList plainGraph(meta_file, graph_data_file);
    ListGraphQueryEngine listGraphEngine(party_info, plainGraph);
    listGraphEngine.add_property(party_info, graph_property_file);  
    timer.end("GraphLoad");
    cmeter.end("GraphLoad_send", get_sending_bytes(party_info));
    cmeter.end("GraphLoad_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    // statistic analysis.
    // run the statistic_func.
    size_t node1 = 0;
    boolIndex snode = boolIndex(node1, role);
    size_t upper_bound = 10;
    boolIndex sec_upper_bound = boolIndex(upper_bound, role);

    cmeter.start("statistic_send", get_sending_bytes(party_info));
    cmeter.start("statistic_recv", get_receiving_bytes(party_info));
    timer.start("statistic");

    aby3::si64Matrix statistics = outting_edge_range_statistics(snode, sec_upper_bound, listGraphEngine);

    timer.end("statistic");
    cmeter.end("statistic_send", get_sending_bytes(party_info));
    cmeter.end("statistic_recv", get_receiving_bytes(party_info));
    clear_sending_reveiving_bytes(party_info);

    if(role == 0) debug_info("Statistics success");

    // print the timer records.
    communication_synchronize(party_info);
    #ifdef MPI_APP
    if((role == 0) && (rank == 0)){
    #else
    if(role == 0){
    #endif
        std::ofstream stream(record_file, std::ios::app);
        stream << "===== Edge list ====" <<std::endl;
        listGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // cmeter.print_total("MB", stream);
        cmeter.print_total_per_party("MB", stream);
    }

    return 0;
}