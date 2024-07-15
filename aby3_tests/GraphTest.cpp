#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-Basic/timer.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-GraphQuery/Graph.h"

using namespace oc;
using namespace aby3;
using namespace std;

int graph_loading_test(oc::CLP &cmd){

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
        debug_info("RUN Graph Loading TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    aby3Info party_info(role, enc, eval, runtime);

    // filenames.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/micro_benchmark/";
    std::string meta_file = "tmp_graph_meta.txt";
    std::string graph_data_file = "tmp_graph_2dpartition.txt";

    // load the graph.
    Graph2d secGraph(graph_data_folder + meta_file, graph_data_folder + graph_data_file, party_info);

    // check the graph.
    secGraph.check_graph(graph_data_folder + meta_file, graph_data_folder + graph_data_file, party_info);

    return 0;
}

int graph_sort_test(oc::CLP &cmd){
    
    TEST_INIT

    if(role == 0){
        debug_info("RUN Multiple Graph Sort TEST");
    }

    // filenames for 2d-partition data.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/multiparty/";
    std::string file_prefix = graph_data_folder + "random_n-16_k-2";
    std::string meta_file = file_prefix + "_meta_multiparty.txt";

    // load the graph.
    Graph2d secGraph(meta_file, file_prefix, party_info, true);
    // check the graph.
    secGraph.check_graph(meta_file, file_prefix, party_info, true);

    // check the edgelist graph.
    file_prefix = graph_data_folder + "random_n-16";
    meta_file = file_prefix + "_edge_list_meta_multiparty.txt";
    ListGraphQueryEngine listEngine(party_info, meta_file, file_prefix, true);

    if(party_info.pIdx == 0) debug_info("list generation ok!");
    listEngine.check_graph(meta_file, file_prefix, party_info, true);

    return 0;
}

int graph_block_fetch_test(oc::CLP& cmd){
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
        debug_info("RUN Graph Block Fetch TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    aby3Info party_info(role, enc, eval, runtime);

    // filenames.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/micro_benchmark/";
    std::string meta_file = "tmp_graph_meta.txt";
    std::string graph_data_file = "tmp_graph_2dpartition.txt";

    // get the plain Graph.
    plainGraph2d pGraph(graph_data_folder + meta_file, graph_data_folder + graph_data_file);

    size_t stash_size = 8;
    size_t pack_size = 4;

    // construct the graph query engine.
    GraphQueryEngine secGraphEngine(party_info, graph_data_folder + meta_file, graph_data_folder + graph_data_file, stash_size, pack_size, stash_size, pack_size);

    // query the graph.
    int starting_node = 10, ending_node = 11;
    int starting_chunk = starting_node / pGraph.k, ending_chunk = ending_node / pGraph.k;
    int target_edge_block = starting_chunk * pGraph.b + ending_chunk;
    boolIndex cipher_edge_block = boolIndex(target_edge_block, role);
    boolIndex cipher_node_chunk = boolIndex(starting_chunk, role);

    // fetch the target edge block.
    aby3::sbMatrix target_edge_block_enc = secGraphEngine.get_edge_block(cipher_edge_block);
    aby3::sbMatrix target_node_chunk_enc = secGraphEngine.get_node_edges(cipher_node_chunk);

    // check the result.
    aby3::i64Matrix test_edge_block;
    aby3::i64Matrix test_node_chunk;
    enc.revealAll(runtime, target_edge_block_enc, test_edge_block).get();
    enc.revealAll(runtime, target_node_chunk_enc, test_node_chunk).get();


    std::vector<int> true_edge_block_vec = pGraph.get_edge_block(starting_chunk, ending_chunk);
    std::vector<int> true_node_chunk_vec = pGraph.get_node_chunk(starting_chunk);

    aby3::i64Matrix true_edge_block(2*pGraph.l, 1);
    aby3::i64Matrix true_node_chunk(2*pGraph.l*pGraph.b, 1);
    for(size_t i=0; i<2*pGraph.l; i++){
        true_edge_block(i, 0) = true_edge_block_vec[i];
    }
    for(size_t i=0; i<2*pGraph.l*pGraph.b; i++){
        true_node_chunk(i, 0) = true_node_chunk_vec[i];
    }

    // check the result.
    if(role == 0){
        check_result("Edge block fetch test", test_edge_block, true_edge_block);
        check_result("Node chunk fetch test", test_node_chunk, true_node_chunk);
    }

    return 0;
}

int basic_graph_query_test(oc::CLP& cmd){
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
        debug_info("RUN Basic Graph Query TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    aby3Info party_info(role, enc, eval, runtime);

    // graph filename, using star graph for testing, 0-other edge exists while other edges do not exist.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/micro_benchmark/";
    std::string meta_file = "star_meta.txt";
    std::string graph_data_file = "star_2dpartition.txt";

    size_t stash_size = 8;
    size_t pack_size = 4;

    // construct the graph query engine.
    GraphQueryEngine GQEngine(party_info, graph_data_folder + meta_file, graph_data_folder + graph_data_file, stash_size, pack_size, stash_size, pack_size);

    // query the graph.
    int starting_node = 0, ending_node = 33;
    int logical_index = GQEngine.get_edge_block_index(starting_node, ending_node);

    boolIndex priv_logical_index = boolIndex(logical_index, role);
    boolIndex priv_starting_node = boolIndex(starting_node, role);
    boolIndex priv_ending_node = boolIndex(ending_node, role);

    // query the edge existence.
    boolShare res1 = edge_existance(priv_starting_node, priv_ending_node, priv_logical_index, GQEngine);

    starting_node = 1; 
    logical_index = GQEngine.get_edge_block_index(starting_node, ending_node);
    priv_logical_index = boolIndex(logical_index, role);
    priv_starting_node = boolIndex(starting_node, role);
    priv_ending_node = boolIndex(ending_node, role);

    boolShare res2 = edge_existance(priv_starting_node, priv_ending_node, priv_logical_index, GQEngine);

    // check the result.
    bool test_res1 = back2plain(role, res1, enc, eval, runtime);
    bool test_res2 = back2plain(role, res2, enc, eval, runtime);

    if(role == 0){
        check_result("Basic graph query exist edge test", test_res1, true);
        check_result("Basic graph query fake edge test", test_res2, false);
    }


    // outting edges count.
    int target_node = 0;
    int target_chunk = GQEngine.get_block_index(target_node);
    boolIndex priv_target_chunk = boolIndex(target_chunk, role);
    boolIndex priv_target_node = boolIndex(target_node, role);
    aby3::si64Matrix outting_edges_count_node0 = outting_edge_count(priv_target_node, priv_target_chunk, GQEngine);


    target_node = 10;
    target_chunk = GQEngine.get_block_index(target_node);
    priv_target_chunk = boolIndex(target_chunk, role);  
    priv_target_node = boolIndex(target_node, role);
    aby3::si64Matrix outting_edges_count_node10 = outting_edge_count(priv_target_node, priv_target_chunk, GQEngine);

    if(role == 0) debug_info("finished two outting edges count");

    // check the result.
    aby3::i64Matrix test_outting_edges_count_node0(1, 1);
    aby3::i64Matrix test_outting_edges_count_node10(1, 1);
    enc.revealAll(runtime, outting_edges_count_node0, test_outting_edges_count_node0).get();
    enc.revealAll(runtime, outting_edges_count_node10, test_outting_edges_count_node10).get();

    if(role == 0){
        check_result("Outting edges count node 0 test", test_outting_edges_count_node0(0, 0), GQEngine.graph->e);
        check_result("Outting edges count node 10 test", test_outting_edges_count_node10(0, 0), 0);
    }

    return 0;
}

int neighbors_find_test(oc::CLP& cmd){

    // set up the environment. 
    TEST_INIT
    if(role == 0){
        debug_info("RUN Basic Graph Neighbor Find Query TEST");
    }

    // graph filename, using star graph for testing, 0-other edge exists while other edges do not exist.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/micro_benchmark/";
    std::string meta_file = "tmp_graph_meta.txt";
    std::string graph_data_file = "tmp_graph_2dpartition.txt";

    size_t stash_size = 8;
    size_t pack_size = 4;

    // construct the graph query engine.
    GraphQueryEngine GQEngine(party_info, graph_data_folder + meta_file, graph_data_folder + graph_data_file, stash_size, pack_size, stash_size, pack_size);

    // query the graph.
    int starting_node = 0;
    int logical_index = GQEngine.get_block_index(starting_node);
    boolIndex priv_logical_index = boolIndex(logical_index, role);
    boolIndex priv_node_index = boolIndex(starting_node, role);

    // query the neighbors.
    aby3::sbMatrix masked_neighbors = outting_neighbors(priv_node_index, priv_logical_index, GQEngine);

    // check the result.
    aby3::i64Matrix test_neighbors(masked_neighbors.rows(), 1);
    enc.revealAll(runtime, masked_neighbors, test_neighbors).get();

    std::vector<int> queried_neighbors;
    for(size_t i=0; i<test_neighbors.rows(); i++){
        if(test_neighbors(i, 0) != 0){
            queried_neighbors.push_back(test_neighbors(i, 0));
        }
    }

    plainGraph2d pGraph(graph_data_folder + meta_file, graph_data_folder + graph_data_file);
    std::vector<int> true_neighbors = pGraph.get_outting_neighbors(starting_node);

    std::sort(queried_neighbors.begin(), queried_neighbors.end(), std::greater<int>());
    std::sort(true_neighbors.begin(), true_neighbors.end(), std::greater<int>());
    if(role == 0){
        check_result("Basic graph query neighbors test", queried_neighbors, true_neighbors);
    }

    // load the data and sort.
    pGraph.per_block_sort();
    GraphQueryEngine GQEngineSort(party_info, pGraph, stash_size, pack_size, stash_size, pack_size);

    // query for the result.
    aby3::sbMatrix masked_neighbors_sort = outting_neighbors_sorted(priv_node_index,    priv_logical_index, GQEngineSort);

    std::ofstream party_fs(PARTY_FILE + std::to_string(role) + ".txt", std::ios::app);
    debug_info("finished neighbors find query.", party_fs);

    aby3::i64Matrix test_neighbors_sort(masked_neighbors_sort.rows(), 1);
    enc.revealAll(runtime, masked_neighbors_sort, test_neighbors_sort).get();
    std::vector<int> queried_neighbors_sort;
    for(size_t i=0; i<test_neighbors_sort.rows(); i++){
        if(test_neighbors_sort(i, 0) != 0){
            queried_neighbors_sort.push_back(test_neighbors_sort(i, 0));
        }
    }

    std::sort(queried_neighbors_sort.begin(), queried_neighbors_sort.end(), std::greater<int>());
    if(role == 0){
        debug_output_vector(queried_neighbors_sort);
        check_result("Basic graph query neighbors sort test", queried_neighbors_sort, true_neighbors);
    }
    
    return 0;
}

// adj graph query tests.
int adj_graph_loading_test(oc::CLP& cmd){

    TEST_INIT
    
    if(role == 0){
        debug_info("RUN AdjGraph loading TEST");
    }

    // filenames.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/micro_benchmark/";
    std::string meta_file = "adj_tmp_edge_list_meta.txt";
    std::string data_file = "adj_tmp_edge_list.txt";

    // load the graph.
    GraphAdj adjGraph(graph_data_folder + meta_file, graph_data_folder + data_file, party_info);

    // check the graph.
    adjGraph.check_graph(graph_data_folder + meta_file, graph_data_folder + data_file, party_info);

    graph_data_folder = "/root/aby3/aby3-GraphQuery/data/multiparty/";
    std::string file_prefix = graph_data_folder + "random_n-16";
    meta_file = file_prefix + "_edge_list_meta_multiparty.txt";

    GraphAdj adjGraph2(meta_file, file_prefix, party_info, true);
    adjGraph2.check_graph(meta_file, file_prefix, party_info, true);

    return 0;
}

int adj_basic_graph_query_test(oc::CLP& cmd){

    TEST_INIT

    if(role == 0){
        debug_info("RUN AdjGraph basic query TEST");
    }

    // filenames.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/micro_benchmark/";
    std::string meta_file = "adj_tmp_edge_list_meta.txt";
    std::string data_file = "adj_tmp_edge_list.txt";

    // load the graph.
    plainGraphAdj plainGraph(graph_data_folder + meta_file, graph_data_folder + data_file);
    plainGraph.generate_adj_list();

    GraphAdj adjGraph(graph_data_folder + meta_file, graph_data_folder + data_file, party_info);

    // build the query engine.
    AdjGraphQueryEngine adjGQEngine(party_info, graph_data_folder + meta_file, graph_data_folder + data_file);

    size_t stash_size = 8;
    size_t pack_size = 4;

    adjGQEngine.edge_oram_initialization(stash_size, pack_size);
    adjGQEngine.node_oram_initialization(stash_size, pack_size);

    // edge existence query.
    int starting_node = 0, ending_node = 1;
    boolIndex priv_starting_node = boolIndex(starting_node, role);
    boolIndex priv_ending_node = boolIndex(ending_node, role);
    
    boolShare res1 = edge_existance(priv_starting_node, priv_ending_node, adjGQEngine);
    bool ref_res1 = plainGraph.edge_existence(starting_node, ending_node);

    // outting edges count query
    aby3::sbMatrix neighbor_count = outting_edge_count(priv_starting_node, adjGQEngine);
    int ref_neighbor_count = plainGraph.outting_neighbors_count(starting_node);
    aby3::si64Matrix arith_neighbor_cound = outting_edge_count_arith(priv_starting_node, adjGQEngine);

    // neighbor find query.
    aby3::sbMatrix neighbors = outting_neighbors(priv_starting_node, adjGQEngine);
    std::vector<int> ref_neighbors_ = plainGraph.get_outting_neighbors(starting_node);
    aby3::i64Matrix ref_neighbors(ref_neighbors_.size(), 1);
    for(size_t i=0; i<ref_neighbors_.size(); i++){
        ref_neighbors(i, 0) = ref_neighbors_[i];
    }

    // check the results.
    bool test_res1 = back2plain(role, res1, enc, eval, runtime);
    aby3::i64Matrix test_neighbor_count = back2plain(role, neighbor_count, enc, eval, runtime);
    aby3::i64Matrix test_neighbor_count_arith;
    enc.revealAll(runtime, arith_neighbor_cound, test_neighbor_count_arith).get();
    aby3::i64Matrix test_neighbors = back2plain(role, neighbors, enc, eval, runtime);

    if(role == 0){
        check_result("AdjGraph basic query edge existence test", test_res1, ref_res1);
        check_result("AdjGraph basic query neighbor count test", test_neighbor_count(0, 0), ref_neighbor_count);
        check_result("AdjGraph basic query neighbor count arith test", test_neighbor_count_arith(0, 0), ref_neighbor_count);
        check_result("AdjGraph basic query neighbor find test", test_neighbors, ref_neighbors);
    }

    return 0;
}

// edge-list query test.
int node_edge_list_basic_graph_query_test(oc::CLP& cmd){

    TEST_INIT

    if(role == 0){
        debug_info("RUN Node-Edge List Basic Query TEST");
    }

    // filenames.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/micro_benchmark/";
    std::string meta_file = "adj_tmp_edge_list_meta.txt";
    std::string data_file = "adj_tmp_edge_list.txt";

    // load the graph.
    plainGraphList plainGraph(graph_data_folder + meta_file, graph_data_folder + data_file);
    ListGraphQueryEngine nodeEdgeListGraph(party_info, graph_data_folder + meta_file, graph_data_folder + data_file);

    // basic tests.
    // 1) edge-existence query.
    int starting_node = 0, ending_node = 1;
    boolIndex priv_starting_node = boolIndex(starting_node, role);
    boolIndex priv_ending_node = boolIndex(ending_node, role);
    
    boolShare res1 = edge_existance(priv_starting_node, priv_ending_node, nodeEdgeListGraph);
    bool ref_res1 = plainGraph.edge_existence(starting_node, ending_node);

    // 2) neighbors count query.
    aby3::sbMatrix neighbor_count = outting_edge_count(priv_starting_node, nodeEdgeListGraph);
    int ref_neighbor_count = plainGraph.outting_neighbors_count(starting_node);

    // 3) neighbors find query.
    aby3::sbMatrix neighbors = outting_neighbors(priv_starting_node, nodeEdgeListGraph);
    std::vector<int> ref_neighbors_ = plainGraph.get_outting_neighbors(starting_node);

    bool test1 = back2plain(role, res1, enc, eval, runtime);
    aby3::i64Matrix test2 = back2plain(role, neighbor_count, enc, eval, runtime);
    aby3::i64Matrix test_neighbors = back2plain(role, neighbors, enc, eval, runtime);
    std::vector<int> test_vec_neighbors(ref_neighbors_.size());
    for(size_t i=0; i<ref_neighbors_.size(); i++){
        test_vec_neighbors[i] = test_neighbors(i, 0);
    }
    std::sort(test_vec_neighbors.begin(), test_vec_neighbors.end());

    // 4) sorted neighbors find query.
    plainGraph.list_sort();
    ListGraphQueryEngine nodeEdgeListGraphSort(party_info, plainGraph);
    aby3::sbMatrix neighbors_sorted = outting_neighbors_sorted(priv_starting_node, nodeEdgeListGraphSort);
    aby3::i64Matrix test_neighbors_sorted = back2plain(role, neighbors_sorted, enc, eval, runtime);
    std::vector<int> test_vec_neighbors_sorted(ref_neighbors_.size());
    for(size_t i=0; i<ref_neighbors_.size(); i++){
        test_vec_neighbors_sorted[i] = test_neighbors_sorted(i, 0);
    }
    std::sort(test_vec_neighbors_sorted.begin(), test_vec_neighbors_sorted.end());

    if(role == 0){
        check_result("Node-Edge List basic query edge existence test", test1, ref_res1);
        check_result("Node-Edge List basic query neighbor count test", test2(0, 0), ref_neighbor_count);
        check_result("Node-Edge List basic query neighbor find test", test_vec_neighbors, ref_neighbors_);
        check_result("Node-Edge List basic query neighbor find sort test", test_vec_neighbors_sorted, ref_neighbors_);
    }

    return 0;
}