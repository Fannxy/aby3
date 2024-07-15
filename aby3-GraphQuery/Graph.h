#pragma once
#include <mpi.h>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include "pGraph.h"
#include "aby3-Basic/SqrtOram.h"
#include "aby3-Basic/Sort.h"

static size_t MAX_COMM_SIZE = 1 << 25;

#define PROPERTY_COMPRESS
// #define MPI

struct aby3Info{
    int pIdx;
    aby3::Sh3Encryptor *enc;
    aby3::Sh3Evaluator *eval;
    aby3::Sh3Runtime *runtime;

    aby3Info(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime):
        pIdx(pIdx), enc(&enc), eval(&eval), runtime(&runtime){}
};

struct Graph2d {
    // each 2d-partation graph includes one node chuck list (length b) and one edge block list (length b^2).
    size_t v, e;
    size_t b, k, l;
    size_t edge_list_size;
    // the node_chuck_list is a list of sbMatrix, each sbMatrix is a list of k nodes.
    std::vector<aby3::sbMatrix> node_chuck_list;

    // the edge_block_list is a list of sbMatrix, each sbMatrix is a list of 2*l node tags (the first l is starting nodes and the last l is ending nodes).
    std::vector<aby3::sbMatrix> edge_block_list;

    // the node_edges_list is a list of sbMatrix, each sbMatrix is a list of 2*b*l node tags (the first b*l is starting nodes and the last b*l is ending nodes).
    std::vector<aby3::sbMatrix> node_edges_list;

    void graph_encryption(const plainGraph2d& plain_graph, aby3Info &party_info){
        v = plain_graph.v;
        e = plain_graph.e;
        b = plain_graph.b;
        k = plain_graph.k;
        l = plain_graph.l;
        edge_list_size = plain_graph.edge_list_size;

        // convert plain graph to secure graph.
        edge_block_list.resize(edge_list_size);

        // encrypt the edge block list in a single round of communication.
        aby3::i64Matrix full_edge_block_list(2*l*edge_list_size, 1);
        aby3::sbMatrix full_edge_block_list_sec(2*l*edge_list_size, BITSIZE);
        for(size_t i=0; i<edge_list_size; i++){
            for(size_t j=0; j<l; j++){
                full_edge_block_list(i*2*l+j, 0) = plain_graph.edge_block_list[i][j][0]; // starting_nodes.
                full_edge_block_list(i*2*l+j+l, 0) = plain_graph.edge_block_list[i][j][1]; // ending_nodes.
            }
        }

        // encrypt the edge block list.
        size_t round = (size_t)ceil(edge_list_size * 2 * l / (double)MAX_COMM_SIZE);
        size_t last_len = edge_list_size * 2 * l - (round - 1) * MAX_COMM_SIZE;
        for(size_t i=0; i<round; i++){
            size_t len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix edge_block_list = full_edge_block_list.block(i*MAX_COMM_SIZE, 0, len, 1);
            aby3::sbMatrix edge_block_list_sec(len, BITSIZE);
            if(party_info.pIdx == 0){
                party_info.enc->localBinMatrix(*(party_info.runtime), edge_block_list, edge_block_list_sec).get();
            }
            else{
                party_info.enc->remoteBinMatrix(*(party_info.runtime), edge_block_list_sec).get();
            }
            std::memcpy(full_edge_block_list_sec.mShares[0].data() + i*MAX_COMM_SIZE, edge_block_list_sec.mShares[0].data(), len * sizeof(full_edge_block_list_sec.mShares[0](0, 0)));
            std::memcpy(full_edge_block_list_sec.mShares[1].data() + i*MAX_COMM_SIZE, edge_block_list_sec.mShares[1].data(), len * sizeof(full_edge_block_list_sec.mShares[0](0, 0)));
        }

        // split the edge block list to the standard format.
        for(size_t i=0; i<edge_list_size; i++){
            edge_block_list[i].resize(2*l, BITSIZE);
            std::memcpy(edge_block_list[i].mShares[0].data(), full_edge_block_list_sec.mShares[0].data() + i*2*l, 2*l * sizeof(full_edge_block_list_sec.mShares[0](0, 0)));
            std::memcpy(edge_block_list[i].mShares[1].data(), full_edge_block_list_sec.mShares[1].data() + i*2*l, 2*l * sizeof(full_edge_block_list_sec.mShares[1](0, 0)));
        }

        return;
    }

    void graph_encrption(const std::vector<plainGraph2d>& plain_graph_vec, aby3Info &party_info){
        // get the configurations.
        int N = plain_graph_vec.size(); // how many graphs -> data_providers.
        e = 0, l = 0;
        for(size_t i=0; i<N; i++){
            // e += plain_graph_vec[i].e;
            l += plain_graph_vec[i].l;
        }
        e = plain_graph_vec[0].e; // e is the same (for debugging).
        v = plain_graph_vec[0].v; // v, b and k are the same.
        b = plain_graph_vec[0].b; 
        k = plain_graph_vec[0].k;
        edge_list_size = b * b;

        // encrypt the origional graphs.
        aby3::i64Matrix full_edge_block_list(l*edge_list_size, 1);
        aby3::sbMatrix full_edge_block_list_sec(l*edge_list_size, BITSIZE);
        size_t offset;
        for(size_t i=0; i<edge_list_size; i++){
            offset = 0;
            aby3::i64Matrix edge_block_list(l, 1);
            // generate l-length edge block.
            for(size_t k=0; k<N; k++){
                size_t lk = plain_graph_vec[k].l;
                for(size_t j=0; j<lk; j++){ // start + end node.
                    edge_block_list(offset+j, 0) = ((int64_t)plain_graph_vec[k].edge_block_list[i][j][0] << 30) | ( (int64_t) plain_graph_vec[k].edge_block_list[i][j][1]);
                }
                offset += lk;
            }
            // append to the edge_block_list.
            full_edge_block_list.block(i*l, 0, l, 1) = edge_block_list;
        }

        // encrypt the edge block list.
        large_data_encryption(party_info.pIdx, full_edge_block_list, full_edge_block_list_sec, *(party_info.enc), *(party_info.runtime));

        // generate the k (source-dest)-nodes list in encrypted format.
        offset = 0;
        std::vector<std::vector<aby3::sbMatrix>> target_edge_list(edge_list_size);
        for(size_t i=0; i<edge_list_size; i++){
            target_edge_list[i].resize(N);
            for(size_t k=0; k<N; k++){
                size_t lk = plain_graph_vec[k].l;
                target_edge_list[i][k].resize(lk, BITSIZE);
                std::memcpy(target_edge_list[i][k].mShares[0].data(), full_edge_block_list_sec.mShares[0].data() + offset, lk * sizeof(full_edge_block_list_sec.mShares[0](0, 0)));
                std::memcpy(target_edge_list[i][k].mShares[1].data(), full_edge_block_list_sec.mShares[1].data() + offset, lk * sizeof(full_edge_block_list_sec.mShares[0](0, 0)));
                offset += lk;
            }
        }

        // merge sort.
        std::vector<aby3::sbMatrix> combined_sorted_edge_list(edge_list_size);
        high_dimensional_odd_even_multi_merge(target_edge_list, combined_sorted_edge_list, party_info.pIdx, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

        // split it to soure, dests and fit with the Graph2d format.
        edge_block_list.resize(edge_list_size);
        for(size_t i=0; i<edge_list_size; i++){
            aby3::sbMatrix edge_block(2*l, BITSIZE);
            for(size_t j=0; j<l; j++){
                edge_block.mShares[0](j, 0) = combined_sorted_edge_list[i].mShares[0](j, 0) >> 30;
                edge_block.mShares[1](j, 0) = combined_sorted_edge_list[i].mShares[1](j, 0) >> 30;
                edge_block.mShares[0](j+l, 0) = combined_sorted_edge_list[i].mShares[0](j, 0) & ((1<<30) - 1);
                edge_block.mShares[1](j+l, 0) = combined_sorted_edge_list[i].mShares[1](j, 0) & ((1<<30) - 1);
            }
            edge_block_list[i] = edge_block;
        }

        // fit with the Graph2d format.
        get_node_edges_list();
        return;
    };

    Graph2d(){}; // default constructor

    Graph2d(const plainGraph2d& plain_graph, aby3Info &party_info){
        graph_encryption(plain_graph, party_info);
        return;
    };

    // Graph2d(const std::string& meta_data_file, const std::string& edge_block_file, aby3Info& party_info) : Graph2d(plainGraph2d(meta_data_file, edge_block_file), party_info){
    //     return;
    // }
    Graph2d(const std::string& meta_data_file, const std::string& edge_block_file, aby3Info& party_info){
        plainGraph2d *plain_graph = new plainGraph2d(meta_data_file, edge_block_file);
        graph_encryption(*plain_graph, party_info);
        delete plain_graph;
        return;
    }

    Graph2d(const std::string& meta_data_file, const std::string& edge_block_file, const std::string& node_chunk_file, aby3Info &party_info) : Graph2d(meta_data_file, edge_block_file, party_info) {
        // plainGraph2d plain_graph(meta_data_file, edge_block_file, node_chunk_file);
        plainGraph2d *plain_graph = new plainGraph2d(meta_data_file, edge_block_file, node_chunk_file);
        graph_encryption(*plain_graph, party_info);
        delete plain_graph;
        return;
    }

    Graph2d(const std::string& meta_data_file, const std::string& file_prefix, aby3Info& party_info, bool multiparty){
        if(!multiparty){
            plainGraph2d plain_graph(meta_data_file, file_prefix);
            graph_encryption(plain_graph, party_info);
            return;
        }

        std::ifstream meta(meta_data_file);
        size_t N_providers;
        meta >> v >> N_providers >> b >> k;
        std::vector<plainGraph2d> plain_graph_vec(N_providers);
        for(size_t i=0; i<N_providers; i++){
            std::string party_edge_block_file = file_prefix + "_2dpartition_party-" + std::to_string(i) + ".txt";
            std::string party_meta_file = file_prefix + "_meta_party-" + std::to_string(i) + ".txt";
            plain_graph_vec[i] = plainGraph2d(party_meta_file, party_edge_block_file);
        }

        // encrypt the multi-party graphs.
        graph_encrption(plain_graph_vec, party_info);

        return;
    }

    Graph2d(const std::vector<plainGraph2d>& plain_graph_vec, aby3Info &party_info){
        graph_encrption(plain_graph_vec, party_info);
        return; 
    }

    void get_node_edges_list(){
        node_edges_list.resize(b);
        for(size_t i=0; i<b; i++){
            node_edges_list[i].resize(b * 2 * l, BITSIZE);
            for(size_t j=0; j<b; j++){
                for(size_t k=0; k<l; k++){
                    // starting node list.
                    node_edges_list[i].mShares[0](j*(l)+k, 0) = edge_block_list[i*b+j].mShares[0](k, 0);
                    node_edges_list[i].mShares[1](j*(l)+k, 0) = edge_block_list[i*b+j].mShares[1](k, 0);
                    // ending node list.
                    node_edges_list[i].mShares[0](j*(l)+k+ (b*l), 0) = edge_block_list[i*b+j].mShares[0](k+l, 0);
                    node_edges_list[i].mShares[1](j*(l)+k+ (b*l), 0) = edge_block_list[i*b+j].mShares[1](k+l, 0);
                }
            }   
        }
        return;
    }

    void check_graph(const plainGraph2d& plain_graph, aby3Info &party_info){
        // check the graph.
        // reveal the secure edge list back to plaintext.
        aby3::i64Matrix start_nodes(edge_list_size * l, 1);
        aby3::i64Matrix end_nodes(edge_list_size * l, 1);

        aby3::sbMatrix start_nodes_sec(edge_list_size * l, 1);
        aby3::sbMatrix end_nodes_sec(edge_list_size * l, 1);

        for(size_t i=0; i<edge_list_size; i++){
            for(size_t j=0; j<l; j++){
                start_nodes_sec.mShares[0](i*l+j, 0) = edge_block_list[i].mShares[0](j, 0);
                start_nodes_sec.mShares[1](i*l+j, 0) = edge_block_list[i].mShares[1](j, 0);
                end_nodes_sec.mShares[0](i*l+j, 0) = edge_block_list[i].mShares[0](j + l, 0);
                end_nodes_sec.mShares[1](i*l+j, 0) = edge_block_list[i].mShares[1](j + l, 0);
            }
        }

        party_info.enc->revealAll(*(party_info.runtime), start_nodes_sec, start_nodes).get();
        party_info.enc->revealAll(*(party_info.runtime), end_nodes_sec, end_nodes).get();
        
        // check whether the secure graph is the same as the plaintext graph.
        bool check_flag = true;
        bool break_outter_loop = false;
        for(size_t i=0; i<edge_list_size; i++){
            for(size_t j=0; j<l; j++){
                if(start_nodes(i*l+j, 0) != plain_graph.edge_block_list[i][j][0] || end_nodes(i*l+j, 0) != plain_graph.edge_block_list[i][j][1]){
                    check_flag = false;
                    break_outter_loop = true;
                    break;
                }
            }
            if(break_outter_loop) break;
        }

        if(party_info.pIdx == 0){
            if(check_flag){
                debug_info("\033[32m The secure graph is the same as the plaintext graph. \033[0m\n");
            }
            else{
                debug_info("\033[31m Error: the secure graph is not the same as the plaintext graph. \033[0m\n");

                // print the secure graph.
                for(size_t i=0; i<edge_list_size; i++){
                    for(size_t j=0; j<l; j++){
                        debug_info("edge_block_list[" + std::to_string(i) + "][" + std::to_string(j) + "]: " + std::to_string(start_nodes(i*l+j, 0)) + " " + std::to_string(end_nodes(i*l+j, 0)));
                        debug_info("plain_graph[" + std::to_string(i) + "][" + std::to_string(j) + "]: " + std::to_string(plain_graph.edge_block_list[i][j][0]) + " " + std::to_string(plain_graph.edge_block_list[i][j][1]));
                    }
                }
            }
        }
        return;
    }

    virtual std::vector<aby3::sbMatrix> get_node_edges_with_property(){
        THROW_RUNTIME_ERROR("This function is NOT supported in Graph2d!");
        std::vector<aby3::sbMatrix> tmp;
        return tmp;
    }

    void check_graph(const std::string& meta_data_file, const std::string& edge_block_file, aby3Info &party_info){
        plainGraph2d plain_graph(meta_data_file, edge_block_file);
        check_graph(plain_graph, party_info);
    }

    void check_graph(const std::string& meta_data_file, const std::string& file_prefix, aby3Info &party_info, bool multiparty){
        if(!multiparty){
            check_graph(meta_data_file, file_prefix, party_info);
            return;
        }

        // generate multiple plain graphs.
        std::ifstream meta(meta_data_file);
        size_t N_providers;
        meta >> v >> N_providers >> b >> k;
        std::vector<plainGraph2d> plain_graph_vec(N_providers);
        for(size_t i=0; i<N_providers; i++){
            std::string party_edge_block_file = file_prefix + "_2dpartition_party-" + std::to_string(i) + ".txt";
            std::string party_meta_file = file_prefix + "_meta_party-" + std::to_string(i) + ".txt";
            plain_graph_vec[i] = plainGraph2d(party_meta_file, party_edge_block_file);
        }

        // construct an integrated plain graph.
        plainGraph2d whole_plain_graph;
        whole_plain_graph.v = v;
        whole_plain_graph.b = b;
        whole_plain_graph.k = k;
        whole_plain_graph.l = 0;
        whole_plain_graph.edge_list_size = b * b;
        for(size_t i=0; i<N_providers; i++){
            whole_plain_graph.l += plain_graph_vec[i].l;
        }

        whole_plain_graph.edge_block_list.resize(edge_list_size);
        for(size_t i=0; i<whole_plain_graph.edge_list_size; i++){
            size_t offset = 0;
            whole_plain_graph.edge_block_list[i].resize(whole_plain_graph.l);
            for(size_t k=0; k<N_providers; k++){
                size_t lk = plain_graph_vec[k].l;
                for(size_t j=0; j<lk; j++){
                    whole_plain_graph.edge_block_list[i][offset + j][0] = plain_graph_vec[k].edge_block_list[i][j][0];
                    whole_plain_graph.edge_block_list[i][offset + j][1] = plain_graph_vec[k].edge_block_list[i][j][1];
                }
                offset += lk;
            }
        }

        whole_plain_graph.per_block_sort();
        check_graph(whole_plain_graph, party_info);
        
        return;
    }
};

struct GraphAdj {
    size_t v, e;
    size_t adj_size;
    std::vector<aby3::sbMatrix> adj_list;

    GraphAdj(){}; // default constructor

    GraphAdj(const std::string& meta_data_file, const std::string& data_file, aby3Info &party_info){
        plainGraphAdj plain_graph(meta_data_file, data_file);

        // check the adj graph.
        plain_graph.generate_adj_list();
        v = plain_graph.v;
        e = plain_graph.e;
        adj_size = plain_graph.adj_list.size();
        adj_list.resize(adj_size);

        aby3::i64Matrix full_adj_matrix(adj_size, 1);
        for(size_t i=0; i<adj_size; i++){
            full_adj_matrix(i, 0) = plain_graph.adj_list[i];
        }

        aby3::sbMatrix full_adj_matrix_sec(adj_size, BITSIZE);

        // encrypt the adj graph
        size_t round = (size_t)ceil(adj_size / (double)MAX_COMM_SIZE);
        size_t last_len = adj_size - (round - 1) * MAX_COMM_SIZE;
        for(size_t i=0; i<round; i++){
            size_t len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix adj_matrix = full_adj_matrix.block(i*MAX_COMM_SIZE, 0, len, 1);
            aby3::sbMatrix adj_matrix_sec(len, BITSIZE);
            if(party_info.pIdx == 0){
                party_info.enc->localBinMatrix(*(party_info.runtime), adj_matrix, adj_matrix_sec).get();
            }
            else{
                party_info.enc->remoteBinMatrix(*(party_info.runtime), adj_matrix_sec).get();
            }
            std::memcpy(full_adj_matrix_sec.mShares[0].data() + i*MAX_COMM_SIZE, adj_matrix_sec.mShares[0].data(), len * sizeof(full_adj_matrix_sec.mShares[0](0, 0)));
            std::memcpy(full_adj_matrix_sec.mShares[1].data() + i*MAX_COMM_SIZE, adj_matrix_sec.mShares[1].data(), len * sizeof(full_adj_matrix_sec.mShares[0](0, 0)));
        }

        for(size_t i=0; i<adj_size; i++){
            adj_list[i].resize(1, BITSIZE);
            adj_list[i].mShares[0](0, 0) = full_adj_matrix_sec.mShares[0](i, 0);
            adj_list[i].mShares[1](0, 0) = full_adj_matrix_sec.mShares[1](i, 0);
        }
        return;
    }

    GraphAdj(const std::string& meta_data_file, const std::string& file_prefix, aby3Info& party_info, bool multiparty){
        if(!multiparty){
            GraphAdj(meta_data_file, file_prefix, party_info);
            return;
        }

        std::ifstream meta(meta_data_file);
        size_t N_providers;
        meta >> v >> N_providers;
        e = 0;
        std::vector<plainGraphAdj> plain_graph_vec(N_providers);
        adj_size = v * v;
        for(size_t i=0; i<N_providers; i++){
            std::string party_meta_file = file_prefix + "_edge_list_meta_party-" + std::to_string(i) + ".txt";
            std::string party_data_file = file_prefix + "_edge_list_party-" + std::to_string(i) + ".txt";
            plain_graph_vec[i] = plainGraphAdj(party_meta_file, party_data_file);
            plain_graph_vec[i].generate_adj_list();
            e += plain_graph_vec[i].e;
        }

        // encrypt the multi-party graphs.
        aby3::i64Matrix full_adj_matrix(N_providers * adj_size, 1); // each party has a v*v adj matrix.
        for(size_t i=0; i<N_providers; i++){
            // std::memcpy(full_adj_matrix.data() + i*adj_size, plain_graph_vec[i].adj_list.data(), adj_size * sizeof(plain_graph_vec[i].adj_list[0]));
            for(size_t j=0; j<adj_size; j++){
                full_adj_matrix(i*adj_size + j, 0) = plain_graph_vec[i].adj_list[j];
            }
        }
        aby3::sbMatrix full_adj_matrix_sec(N_providers * adj_size, BITSIZE);

        // encrypt the adj graph.
        large_data_encryption(party_info.pIdx, full_adj_matrix, full_adj_matrix_sec, *(party_info.enc), *(party_info.runtime));
        
        size_t mix_len = N_providers;
        if(!checkPowerOfTwo(mix_len)){
            // THROW_RUNTIME_ERROR("The number of data providers must be power of 2.");
            // pad to the power of 2.
            mix_len = roundUpToPowerOfTwo(N_providers);
            size_t left_size = mix_len - N_providers;

            aby3::sbMatrix _adj_mixup(left_size * adj_size, BITSIZE);
            _adj_mixup.mShares[0].setZero();
            _adj_mixup.mShares[1].setZero();
            full_adj_matrix_sec.resize(mix_len * adj_size, BITSIZE);

            std::copy(_adj_mixup.mShares[0].begin(), _adj_mixup.mShares[0].end(), full_adj_matrix_sec.mShares[0].data() + N_providers * adj_size);
            std::copy(_adj_mixup.mShares[1].begin(), _adj_mixup.mShares[1].end(), full_adj_matrix_sec.mShares[1].data() + N_providers * adj_size);
        }

        size_t round = log2(mix_len);
        size_t mid_len = mix_len * adj_size;

        for(size_t i=0; i<round; i++){
            mid_len /= 2;

            aby3::sbMatrix _left_val(mid_len, BITSIZE);
            aby3::sbMatrix _right_val(mid_len, BITSIZE);
            aby3::sbMatrix _res_val(mid_len, BITSIZE);

            std::copy(full_adj_matrix_sec.mShares[0].data(), full_adj_matrix_sec.mShares[0].data() + mid_len, _left_val.mShares[0].data());
            std::copy(full_adj_matrix_sec.mShares[1].data(), full_adj_matrix_sec.mShares[1].data() + mid_len, _left_val.mShares[1].data());
            std::copy(full_adj_matrix_sec.mShares[0].data() + mid_len, full_adj_matrix_sec.mShares[0].end(), _right_val.mShares[0].data());
            std::copy(full_adj_matrix_sec.mShares[1].data() + mid_len, full_adj_matrix_sec.mShares[1].end(), _right_val.mShares[1].data());

            // bool_cipher_add(party_info.pIdx, _left_val, _right_val, _res_val, *(party_info.enc), *(party_info.eval), *(party_info.runtime));
            // rounded addition
            size_t round_inner = (size_t)ceil(mid_len / (double)MAX_COMM_SIZE);

            size_t last_len = mid_len - (round_inner - 1) * MAX_COMM_SIZE;
            for(size_t j=0; j<round_inner; j++){
                size_t len = (j == round_inner - 1) ? last_len : MAX_COMM_SIZE;

                aby3::sbMatrix left_val(len, BITSIZE);
                aby3::sbMatrix right_val(len, BITSIZE);
                std::memcpy(left_val.mShares[0].data(), _left_val.mShares[0].data() + j*MAX_COMM_SIZE, len * sizeof(_left_val.mShares[0](0, 0)));
                std::memcpy(left_val.mShares[1].data(), _left_val.mShares[1].data() + j*MAX_COMM_SIZE, len * sizeof(_left_val.mShares[1](0, 0)));
                std::memcpy(right_val.mShares[0].data(), _right_val.mShares[0].data() + j*MAX_COMM_SIZE, len * sizeof(_right_val.mShares[0](0, 0)));
                std::memcpy(right_val.mShares[1].data(), _right_val.mShares[1].data() + j*MAX_COMM_SIZE, len * sizeof(_right_val.mShares[1](0, 0)));

                aby3::sbMatrix res_val(len, BITSIZE);
                
                bool_cipher_add(party_info.pIdx, left_val, right_val, res_val, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

                std::memcpy(_res_val.mShares[0].data() + j*MAX_COMM_SIZE, res_val.mShares[0].data(), len * sizeof(res_val.mShares[0](0, 0)));
                std::memcpy(_res_val.mShares[1].data() + j*MAX_COMM_SIZE, res_val.mShares[1].data(), len * sizeof(res_val.mShares[1](0, 0)));
            }

            full_adj_matrix_sec = _res_val;
        }
        
        adj_list.resize(adj_size);
        for(size_t i=0; i<adj_size; i++){
            adj_list[i].resize(1, BITSIZE);
            adj_list[i].mShares[0](0, 0) = full_adj_matrix_sec.mShares[0](i, 0);
            adj_list[i].mShares[1](0, 0) = full_adj_matrix_sec.mShares[1](i, 0);
        }

        return;
    }

    void check_graph(plainGraphAdj plain_graph, aby3Info &party_info){
        // check the graph.
        plain_graph.generate_adj_list();
        // reveal the secure edge list back to plaintext.
        aby3::i64Matrix adj_matrix(adj_size, 1);
        aby3::sbMatrix adj_matrix_sec(adj_size, 1);

        for(size_t i=0; i<adj_size; i++){
            adj_matrix_sec.mShares[0](i, 0) = adj_list[i].mShares[0](0, 0);
            adj_matrix_sec.mShares[1](i, 0) = adj_list[i].mShares[1](0, 0);
        }

        large_data_decryption(party_info.pIdx, adj_matrix_sec, adj_matrix, *(party_info.enc), *(party_info.runtime));
        
        // check whether the secure graph is the same as the plaintext graph.
        bool check_flag = true;
        for(size_t i=0; i<adj_size; i++){
            if(adj_matrix(i, 0) != plain_graph.adj_list[i]){
                check_flag = false;
                break;
            }
        }

        if(party_info.pIdx == 0){
            if(check_flag){
                debug_info("\033[32m The secure graph is the same as the plaintext graph. \033[0m\n");
            }
            else{
                debug_info("\033[31m Error: the secure graph is not the same as the plaintext graph. \033[0m\n");

                // print the secure graph.
                debug_info("plain graph: ");
                for(size_t i=0; i<v; i++){
                    std::string info = " ";
                    for(size_t j=0; j<v; j++){
                        info += std::to_string(plain_graph.adj_list[i*v+j]) + " ";
                    }
                    debug_info(info);
                }

                debug_info("secure graph: ");
                for(size_t i=0; i<v; i++){
                    std::string info = " ";
                    for(size_t j=0; j<v; j++){
                        info += std::to_string(adj_matrix(i*v+j, 0)) + " ";
                    }
                    debug_info(info);
                }
            }
        }
        return;
    }

    void check_graph(const std::string& meta_data_file, const std::string& data_file, aby3Info &party_info){
        plainGraphAdj plain_graph(meta_data_file, data_file);
        plain_graph.generate_adj_list();

        // reveal the secure edge list back to plaintext.
        aby3::i64Matrix adj_matrix(plain_graph.adj_list.size(), 1);
        aby3::sbMatrix adj_matrix_sec(plain_graph.adj_list.size(), 1);

        for(size_t i=0; i<plain_graph.adj_list.size(); i++){
            adj_matrix_sec.mShares[0](i, 0) = adj_list[i].mShares[0](0, 0);
            adj_matrix_sec.mShares[1](i, 0) = adj_list[i].mShares[1](0, 0);
        }

        party_info.enc->revealAll(*(party_info.runtime), adj_matrix_sec, adj_matrix).get();
        
        // check whether the secure graph is the same as the plaintext graph.
        bool check_flag = true;
        for(size_t i=0; i<plain_graph.adj_list.size(); i++){
            if(adj_matrix(i, 0) != plain_graph.adj_list[i]){
                check_flag = false;
                break;
            }
        }

        if(party_info.pIdx == 0){
            if(check_flag){
                debug_info("\033[32m The secure graph is the same as the plaintext graph. \033[0m\n");
            }
            else{
                debug_info("\033[31m Error: the secure graph is not the same as the plaintext graph. \033[0m\n");
            }
        }
    }

    void check_graph(const std::string& meta_data_file, const std::string& file_prefix, aby3Info& party_info, bool multiparty){
        
        if(!multiparty){
            check_graph(meta_data_file, file_prefix, party_info);
            return;
        }
        plainGraphAdj whole_plain_graph;
        std::ifstream meta(meta_data_file);
        size_t N_providers;
        meta >> v >> N_providers;
        whole_plain_graph.v = v;

        for(size_t i=0; i<N_providers; i++){
            std::string party_meta_file = file_prefix + "_edge_list_meta_party-" + std::to_string(i) + ".txt";
            std::string party_data_file = file_prefix + "_edge_list_party-" + std::to_string(i) + ".txt";
            plainGraphAdj plain_graph(party_meta_file, party_data_file);
            
            whole_plain_graph.starting_node_list.insert(whole_plain_graph.starting_node_list.end(), plain_graph.starting_node_list.begin(), plain_graph.starting_node_list.end());
            whole_plain_graph.ending_node_list.insert(whole_plain_graph.ending_node_list.end(), plain_graph.ending_node_list.begin(), plain_graph.ending_node_list.end());
        }

        check_graph(whole_plain_graph, party_info);

        return;
    }

};

class PropertyGraph2d : public Graph2d {

    std::vector<std::vector<aby3::sbMatrix>> node_edges_properties;
    std::vector<std::vector<aby3::sbMatrix>> properties;


    public:
        PropertyGraph2d(const std::string& meta_data_file, const std::string& edge_block_file, aby3Info &party_info) : Graph2d(meta_data_file, edge_block_file, party_info){
            return;
        }

        void add_property(const std::string& property_file, aby3Info &party_info){
            std::ifstream property(property_file);
            if(!property.is_open()){
                THROW_RUNTIME_ERROR("The property file " + property_file + " does not exist.");
            }

            std::vector<int> plain_property(edge_list_size * l);
            
            for(size_t i=0; i<edge_list_size; i++){
                for(size_t j=0; j<l; j++){
                    property >> plain_property[i*l+j];
                }
            }

            aby3::i64Matrix property_matrix(edge_list_size * l, 1);
            for(size_t i=0; i<edge_list_size; i++){
                for(size_t j=0; j<l; j++){
                    property_matrix(i*l+j, 0) = plain_property[i*l+j];
                }
            }

            plain_property.resize(0);

            aby3::sbMatrix property_sec(edge_list_size, BITSIZE);
            large_data_encryption(party_info.pIdx, property_matrix, property_sec, *(party_info.enc), *(party_info.runtime));

            std::vector<aby3::sbMatrix> property_list(edge_list_size);
            for(size_t i=0; i<edge_list_size; i++){
                property_list[i].resize(l, BITSIZE);
                for(size_t j=0; j<l; j++){
                    property_list[i].mShares[0](0, 0) = property_sec.mShares[0](i, 0);
                    property_list[i].mShares[1](0, 0) = property_sec.mShares[1](i, 0);
                }
            }

            property_sec.resize(0, 0);

            properties.push_back(property_list);
            return;
        }

        void get_node_edges_property(){

            for(auto& property_list : properties){
                std::vector<aby3::sbMatrix> node_edges_property(b);
                for(size_t i=0; i<b; i++){
                    node_edges_property[i].resize(b * l, BITSIZE);
                    for(size_t j=0; j<b; j++){
                        for(size_t k=0; k<l; k++){
                            node_edges_property[i].mShares[0](j*l+k, 0) = property_list[i*b+j].mShares[0](k, 0);
                            node_edges_property[i].mShares[1](j*l+k, 0) = property_list[i*b+j].mShares[1](k, 0);
                        }
                    }
                }
                node_edges_properties.push_back(node_edges_property);
            }
            properties.resize(0);
            return;
        }

        std::vector<aby3::sbMatrix> get_node_edges_with_property(){
            std::vector<aby3::sbMatrix> node_edges_with_property(b);

            #ifndef PROPERTY_COMPRESS
            for(size_t i=0; i<b; i++){
                node_edges_with_property[i].resize(b * 3 * l, BITSIZE);
                for(size_t j=0; j<b; j++){
                    for(size_t k=0; k<l; k++){
                        // starting node list.
                        node_edges_with_property[i].mShares[0](j*(l)+k, 0) = edge_block_list[i*b+j].mShares[0](k, 0);
                        node_edges_with_property[i].mShares[1](j*(l)+k, 0) = edge_block_list[i*b+j].mShares[1](k, 0);
                        // ending node list.
                        node_edges_with_property[i].mShares[0](j*(l)+k+ (b*l), 0) = edge_block_list[i*b+j].mShares[0](k+l, 0);
                        node_edges_with_property[i].mShares[1](j*(l)+k+ (b*l), 0) = edge_block_list[i*b+j].mShares[1](k+l, 0);
                        // property_list.
                        node_edges_with_property[i].mShares[0](j*(l)+k+ (2*b*l), 0) = node_edges_properties[0][i].mShares[0](j*l+k, 0);
                        node_edges_with_property[i].mShares[1](j*(l)+k+ (2*b*l), 0) = node_edges_properties[0][i].mShares[1](j*l+k, 0);
                    }
                }
            }
            #endif

            #ifdef PROPERTY_COMPRESS

            debug_info("Compress the properties.");

            for(size_t i=0; i<b; i++){
                node_edges_with_property[i].resize(b*2*l, BITSIZE);
                for(size_t j=0; j<b; j++){
                    for(size_t k=0; k<l; k++){
                        node_edges_with_property[i].mShares[0](j*l+k, 0) = (edge_block_list[i*b+j].mShares[0](k, 0) << 30) | (edge_block_list[i*b+j].mShares[0](k+l, 0));
                        node_edges_with_property[i].mShares[1](j*l+k, 0) = (edge_block_list[i*b+j].mShares[1](k, 0) << 30) | (edge_block_list[i*b+j].mShares[1](k+l, 0));

                        // properties, in the second l.
                        node_edges_with_property[i].mShares[0](j*l+k + b*l, 0) = node_edges_properties[0][i].mShares[0](j*l+k, 0);
                        node_edges_with_property[i].mShares[1](j*l+k + b*l, 0) = node_edges_properties[0][i].mShares[1](j*l+k, 0);
                    }
                }
            }

            #endif

            node_edges_properties.resize(0);
            return node_edges_with_property;
        }

};

class GraphQueryEngine{
    public:
        Graph2d *graph;
        ABY3SqrtOram *edge_block_oram;
        ABY3SqrtOram *node_edges_oram;
        ABY3SqrtOram *node_edges_property_oram;
        aby3Info *party_info;

        size_t logb, logk;
        size_t se = 0, pe = 0, sn = 0, pn = 0;

        GraphQueryEngine(){}

        GraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& edge_block_file){
            graph = new Graph2d(meta_data_file, edge_block_file, party_info);
            this->party_info = &party_info;

            if(!checkPowerOfTwo(graph->b) || !checkPowerOfTwo(graph->k)){
                if(party_info.pIdx == 0){
                    debug_info("b = " + std::to_string(graph->b) + " k = " + std::to_string(graph->k) + "\n");
                }
                THROW_RUNTIME_ERROR("The block size and the chunk size must be power of 2.");
            }

            logb = log2(graph->b);
            logk = log2(graph->k);

            return;
        }

        GraphQueryEngine(plainGraph2d& plain_graph, aby3Info& party_info){
            graph = new Graph2d(plain_graph, party_info);
            this->party_info = &party_info;

            if(!checkPowerOfTwo(graph->b) || !checkPowerOfTwo(graph->k)){
                THROW_RUNTIME_ERROR("The block size and the chunk size must be power of 2.");
            }

            logb = log2(graph->b);
            logk = log2(graph->k);

            return;
        }

        GraphQueryEngine(aby3Info& party_info, const std::string& meta_data_file, const std::string& file_prefix, bool multiparty){
            if(!multiparty){
                GraphQueryEngine(party_info, meta_data_file, file_prefix);
                return;
            }
            graph = new Graph2d(meta_data_file, file_prefix, party_info, multiparty);
            this->party_info = &party_info;

            if(!checkPowerOfTwo(graph->b) || !checkPowerOfTwo(graph->k)){
                THROW_RUNTIME_ERROR("The block size and the chunk size must be power of 2.");
            }

            logb = log2(graph->b);
            logk = log2(graph->k);

            return;
        }

        void set_graph(PropertyGraph2d &property_graph){
            this->graph = &property_graph;

            if(!checkPowerOfTwo(graph->b) || !checkPowerOfTwo(graph->k)){
                if(party_info->eval == 0){
                    debug_info("b = " + std::to_string(graph->b) + " k = " + std::to_string(graph->k) + "\n");
                }
                THROW_RUNTIME_ERROR("The block size and the chunk size must be power of 2.");
            }


            logb = log2(graph->b);
            logk = log2(graph->k);

            if(party_info->pIdx == 0){
                debug_info("logb = " + std::to_string(logb) + " logk = " + std::to_string(logk) + "\n");
            }

            return;
        }

        void edge_block_oram_initialization(const int stash_size, const int pack_size){
            edge_block_oram = new ABY3SqrtOram(graph->edge_list_size, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));
            edge_block_oram->initiate(graph->edge_block_list);
            this->se = edge_block_oram->S;
            this->pe = edge_block_oram->pack;
            return;
        }

        void node_edges_oram_initialization(const int stash_size, const int pack_size){
            // iniitialize the node edges oram.
            node_edges_oram = new ABY3SqrtOram(graph->b, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));

            // construct the node edges data.
            graph->get_node_edges_list();
            node_edges_oram->initiate(graph->node_edges_list);
            this->sn = node_edges_oram->S;
            this->pn = node_edges_oram->pack;
            return;
        }

        void node_edges_property_oram_initialization(const int stash_size, const int pack_size){

            if(party_info->pIdx == 0) {
                debug_info("before node_edges_property_oram initialization!");
            }

            node_edges_property_oram = new ABY3SqrtOram(graph->b, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));

            if(party_info->pIdx == 0){
                debug_info("before get property!");
            }

            std::vector<aby3::sbMatrix> property_graph_node_edges = graph->get_node_edges_with_property();        

            node_edges_property_oram->initiate(property_graph_node_edges);

            this->sn = node_edges_property_oram->S;
            this->pn = node_edges_property_oram->pack;

            return;
        }

        GraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& edge_block_file, const size_t edge_oram_stash_size, const size_t edge_oram_pack_size, const size_t node_oram_stash_size, const size_t node_oram_pack_size) : GraphQueryEngine(party_info, meta_data_file, edge_block_file)
        {
            edge_block_oram_initialization(edge_oram_stash_size, edge_oram_pack_size);
            node_edges_oram_initialization(node_oram_stash_size, node_oram_pack_size);
            return;
        }

        GraphQueryEngine(aby3Info& party_info, plainGraph2d& plain_graph, const size_t edge_oram_stash_size, const size_t edge_oram_pack_size, const size_t node_oram_stash_size, const size_t node_oram_pack_size) : GraphQueryEngine(plain_graph, party_info)
        {
            edge_block_oram_initialization(edge_oram_stash_size, edge_oram_pack_size);
            node_edges_oram_initialization(node_oram_stash_size, node_oram_pack_size);
            return;
        }

        ~GraphQueryEngine(){
            if(graph != nullptr) {
                debug_info("delete graph!");
                // delete graph;
                graph = nullptr;
                // debug_info("success delete graph!");
            }
            if(edge_block_oram != nullptr) {
                // delete edge_block_oram;
                // debug_info("delete edge_block_oram!");
                edge_block_oram = nullptr;
                // debug_info("success delete edge_block_oram!");
            }
            if(node_edges_oram != nullptr) {
                // delete node_edges_oram;
                // debug_info("delete node_edges_oram!");
                node_edges_oram = nullptr;
            }
            if(node_edges_property_oram != nullptr) {
                // delete node_edges_property_oram;
                // debug_info("delete node_edges_property_oram!");
                node_edges_property_oram = nullptr;
            }
        }

        void rebuild(aby3Info& party_info, plainGraph2d& plain_graph){
            if(graph != nullptr) {
                delete graph;
                graph = nullptr;
            }
            if(edge_block_oram != nullptr) {
                delete edge_block_oram;
                edge_block_oram = nullptr;
            }
            if(node_edges_oram != nullptr) {
                delete node_edges_oram;
                node_edges_oram = nullptr;
            }

            graph = new Graph2d(plain_graph, party_info);
            this->party_info = &party_info;

            if(!checkPowerOfTwo(graph->b) || !checkPowerOfTwo(graph->k)){
                THROW_RUNTIME_ERROR("The block size and the chunk size must be power of 2.");
            }

            logb = log2(graph->b);
            logk = log2(graph->k);

            return;
        }

        boolIndex get_block_index(boolIndex node_index){
            boolIndex block_index, left_size;
            bool_shift(party_info->pIdx, node_index, logk, block_index, true); // right shift
            return block_index;
        }

        boolIndex get_edge_block_index(boolIndex starting_node, boolIndex ending_node){
            // this function should be called in plaintext phase.
            // this function is only used for test.
            boolIndex starting_block_index = get_block_index(starting_node);
            boolIndex ending_block_index = get_block_index(ending_node);
            
            boolIndex block_index;
            bool_shift(party_info->pIdx, block_index, logb, starting_block_index, false); // left shift
            aby3::sbMatrix edge_block_index_mat;
            aby3::sbMatrix starting_block_index_mat = starting_block_index.to_matrix();
            aby3::sbMatrix ending_block_index_mat = ending_block_index.to_matrix();
            bool_cipher_add(party_info->pIdx, starting_block_index_mat, ending_block_index_mat, edge_block_index_mat, *(party_info->enc), *(party_info->eval), *(party_info->runtime));

            block_index.from_matrix(edge_block_index_mat);

            return block_index;
        }

        int get_block_index(int node_index){
            return node_index >> logk;
        }

        int get_edge_block_index(int starting_node, int ending_node){
            return ((starting_node >> logk) * graph->b) + (ending_node >> logk);
        }

        aby3::sbMatrix get_edge_block(boolIndex edge_block_idx){
            return edge_block_oram->access(edge_block_idx);
        }

        aby3::sbMatrix get_node_edges(boolIndex node_idx){
            return node_edges_oram->access(node_idx);
        }

        aby3::sbMatrix get_node_edges_property(boolIndex node_idx){
            return node_edges_property_oram->access(node_idx);
        }

        void print_configs(std::ostream& stream){
            stream << "v : " << graph->v << std::endl;
            stream << "e : " << graph->e << std::endl;
            stream << "b : " << graph->b << std::endl;
            stream << "l : " << graph->l << std::endl;
            stream << "k : " << graph->k << std::endl;
            stream << "se : " << this->se << std::endl;
            stream << "pe : " << this->pe << std::endl;
            stream << "sn : " << this->sn << std::endl;
            stream << "pn : " << this->pn << std::endl;
        }
};

class AdjGraphQueryEngine{

public:
    GraphAdj *graph;
    ABY3SqrtOram *edge_oram;
    ABY3SqrtOram *node_oram;
    aby3Info *party_info;

    size_t v;

    AdjGraphQueryEngine(){}

    AdjGraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& data_file){
        graph = new GraphAdj(meta_data_file, data_file, party_info);
        this->party_info = &party_info;
        v = graph->v;

        if(!checkPowerOfTwo(v)){ // as we only support 2^m indexed ORAM for efficiency.
            THROW_RUNTIME_ERROR("The number of nodes must be power of 2.");
        }
        return;
    }

    AdjGraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& file_prefix, bool multiparty){
        if(!multiparty){
            THROW_RUNTIME_ERROR("The multiparty mode is not enabled.");
            return;
        }
        graph = new GraphAdj(meta_data_file, file_prefix, party_info, true);
        this->party_info = &party_info;
        v = graph->v;

        if(!checkPowerOfTwo(v)){ // as we only support 2^m indexed ORAM for efficiency.
            THROW_RUNTIME_ERROR("The number of nodes must be power of 2.");
        }
        return;
    }


    void edge_oram_initialization(const int stash_size, const int pack_size){
        edge_oram = new ABY3SqrtOram(graph->adj_size, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));
        edge_oram->initiate(graph->adj_list);
        return;
    }

    void node_oram_initialization(const int stash_size, const int pack_size){

        // organize the node-edges data structure.
        std::vector<aby3::sbMatrix> node_edges_list(v);

        for(size_t i=0; i<v; i++){
            aby3::sbMatrix node_edges(v, BITSIZE);
            for(size_t j=0; j<v; j++){
                node_edges.mShares[0](j, 0) = graph->adj_list[i*v+j].mShares[0](0, 0);
                node_edges.mShares[1](j, 0) = graph->adj_list[i*v+j].mShares[1](0, 0);
            }
            node_edges_list[i] = node_edges;
        }

        node_oram = new ABY3SqrtOram(v, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));
        node_oram->initiate(node_edges_list);

        return;
    }

    aby3::sbMatrix get_target_edge(boolIndex edge_index){
        return edge_oram->access(edge_index);
    }

    aby3::sbMatrix get_target_node(boolIndex node_index){
        return node_oram->access(node_index);
    }

    boolIndex get_logical_edge_index(boolIndex start_node, boolIndex end_node){
        boolIndex logical_edge_index;
        bool_shift(party_info->pIdx, start_node, log2(v), logical_edge_index, false); // left shift
        aby3::sbMatrix edge_index_mat;
        aby3::sbMatrix starting_index_mat = logical_edge_index.to_matrix();
        aby3::sbMatrix ending_index_mat = end_node.to_matrix();
        bool_cipher_add(party_info->pIdx, starting_index_mat, ending_index_mat, edge_index_mat, *(party_info->enc), *(party_info->eval), *(party_info->runtime));

        logical_edge_index.from_matrix(edge_index_mat);

        return logical_edge_index;
    
    }

    void print_configs(std::ostream& stream){
        stream << "v : " << graph->v << std::endl;
        stream << "e : " << graph->e << std::endl;
        stream << "se : " << edge_oram->S << std::endl;
        stream << "pe : " << edge_oram->pack << std::endl;
        stream << "sn : " << node_oram->S << std::endl;
        stream << "pn : " << node_oram->pack << std::endl;
    }

};

class ListGraphQueryEngine{
    public:
        aby3Info *party_info;
        size_t v;
        size_t e;
        aby3::sbMatrix starting_node_list;
        aby3::sbMatrix ending_node_list;
        aby3::sbMatrix edge_property_list;

        void graph_encryption(const plainGraphList& plain_graph, aby3Info &party_info){
            v = plain_graph.v;
            e = plain_graph.e;
            this->party_info = &party_info;

            // convert the plain graph to secure graph.
            aby3::i64Matrix plain_starting_nodes(e, 1);
            aby3::i64Matrix plain_ending_nodes(e, 1);
            for(size_t i=0; i<e; i++){
                plain_starting_nodes(i, 0) = plain_graph.starting_node_list[i];
                plain_ending_nodes(i, 0) = plain_graph.ending_node_list[i];
            }
            
            starting_node_list.resize(e, BITSIZE);
            ending_node_list.resize(e, BITSIZE);

            // encrypt the nodes.
            large_data_encryption(party_info.pIdx, plain_starting_nodes, starting_node_list, *(party_info.enc), *(party_info.runtime));
            large_data_encryption(party_info.pIdx, plain_ending_nodes, ending_node_list, *(party_info.enc), *(party_info.runtime));

            return;
        }

        void graph_encryption(const std::vector<plainGraphList>& plain_graph_vec, aby3Info &party_info){
            size_t N = plain_graph_vec.size();
            this->v = plain_graph_vec[0].v;
            e = 0;
            for(size_t i=0; i<N; i++){
                e += plain_graph_vec[i].e;
            }

            aby3::i64Matrix full_edges(e, 1);
            aby3::sbMatrix full_edges_sec(e, BITSIZE);

            size_t offset = 0;
            for(size_t i=0; i<N; i++){
                for(size_t j=0; j<plain_graph_vec[i].e; j++){
                    full_edges(offset + j, 0) = ((int64_t)plain_graph_vec[i].starting_node_list[j] << 30) | ((int64_t) plain_graph_vec[i].ending_node_list[j]);
                }
                offset += plain_graph_vec[i].e;
            }

            // encrypt the edges.
            large_data_encryption(party_info.pIdx, full_edges, full_edges_sec, *(party_info.enc), *(party_info.runtime));

            // split the secure edges to the secure starting and ending nodes.
            std::vector<aby3::sbMatrix> vec_edge_list(N);
            offset = 0;
            for(size_t i=0; i<N; i++){

                vec_edge_list[i].resize(plain_graph_vec[i].e, BITSIZE);

                std::memcpy(vec_edge_list[i].mShares[0].data(), full_edges_sec.mShares[0].data() + offset, plain_graph_vec[i].e * sizeof(full_edges_sec.mShares[0](0, 0)));
                std::memcpy(vec_edge_list[i].mShares[1].data(), full_edges_sec.mShares[1].data() + offset, plain_graph_vec[i].e * sizeof(full_edges_sec.mShares[0](0, 0)));

                offset += plain_graph_vec[i].e;
            }

            // merge sort.
            aby3::sbMatrix sorted_edge_list(e, BITSIZE);
            odd_even_multi_merge(vec_edge_list, sorted_edge_list, party_info.pIdx, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

            // split it to soure, dests and fit with the Graph2d format.
            this->starting_node_list.resize(e, BITSIZE);
            this->ending_node_list.resize(e, BITSIZE);
            for(size_t i=0; i<e; i++){
                starting_node_list.mShares[0](i, 0) = sorted_edge_list.mShares[0](i, 0) >> 30;
                starting_node_list.mShares[1](i, 0) = sorted_edge_list.mShares[1](i, 0) >> 30;
                ending_node_list.mShares[0](i, 0) = sorted_edge_list.mShares[0](i, 0) & ((1<<30) - 1);
                ending_node_list.mShares[1](i, 0) = sorted_edge_list.mShares[1](i, 0) & ((1<<30) - 1);
            }

            return;
        }

        ListGraphQueryEngine(){}

        ListGraphQueryEngine(aby3Info& party_info, const plainGraphList& plain_graph){
            graph_encryption(plain_graph, party_info);
            return;
        }

        ListGraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& data_file) : ListGraphQueryEngine(party_info, plainGraphList(meta_data_file, data_file)){
            return;
        }

        ListGraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& file_prefix, bool multiparty){
            if(!multiparty){
                graph_encryption(plainGraphList(meta_data_file, file_prefix), party_info);
                return;
            }

            std::vector<plainGraphList> plain_graph_vec;
            size_t N_providers;
            std::ifstream meta(meta_data_file);
            meta >> v >> N_providers;
            for(size_t i=0; i<N_providers; i++){
                std::string party_meta_file = file_prefix + "_edge_list_meta_party-" + std::to_string(i) + ".txt";
                std::string party_data_file = file_prefix + "_edge_list_party-" + std::to_string(i) + ".txt";
                plain_graph_vec.push_back(plainGraphList(party_meta_file, party_data_file));
            }

            graph_encryption(plain_graph_vec, party_info);

            return;
        }

        void add_property(aby3Info &party_info, const std::string& edge_property_file){
            std::vector<int> edge_property(this->e);
            std::ifstream property(edge_property_file);
            if(!property.is_open()){
                THROW_RUNTIME_ERROR("The edge property file " + edge_property_file + " does not exist.");
            }
            for(size_t i=0; i<e; i++){
                property >> edge_property[i];
            }

            aby3::i64Matrix property_matrix(e, 1);
            for(size_t i=0; i<e; i++){
                property_matrix(i, 0) = edge_property[i];
            }
            edge_property_list.resize(e, BITSIZE);
            large_data_encryption(party_info.pIdx, property_matrix, edge_property_list, *(party_info.enc), *(party_info.runtime));
            return;
        }

        void rebuild(aby3Info& party_info, plainGraphList& plain_graph){
            v = plain_graph.v;
            e = plain_graph.e;
            this->party_info = &party_info;

            // convert the plain graph to secure graph.
            aby3::i64Matrix plain_starting_nodes(e, 1);
            aby3::i64Matrix plain_ending_nodes(e, 1);
            for(size_t i=0; i<e; i++){
                plain_starting_nodes(i, 0) = plain_graph.starting_node_list[i];
                plain_ending_nodes(i, 0) = plain_graph.ending_node_list[i];
            }
            
            starting_node_list.resize(e, BITSIZE);
            ending_node_list.resize(e, BITSIZE);

            // encrypt the nodes.
            large_data_encryption(party_info.pIdx, plain_starting_nodes, starting_node_list, *(party_info.enc), *(party_info.runtime));
            large_data_encryption(party_info.pIdx, plain_ending_nodes, ending_node_list, *(party_info.enc), *(party_info.runtime));
            return;
        }

        void print_configs(std::ostream& stream){
            stream << "v : " << v << std::endl;
            stream << "e : " << e << std::endl;
        }

        void check_graph(plainGraphList& plain_graph, aby3Info& party_info){
            
            plain_graph.list_sort();

            // reveal the secure edge list back to plaintext.
            aby3::i64Matrix start_nodes(e, 1);
            aby3::i64Matrix end_nodes(e, 1);
            large_data_decryption(party_info.pIdx, starting_node_list, start_nodes, *(party_info.enc), *(party_info.runtime));
            large_data_decryption(party_info.pIdx, ending_node_list, end_nodes, *(party_info.enc), *(party_info.runtime));

            // check whether the secure graph is the same as the plaintext graph.
            bool check_flag = true;
            bool break_outter_loop = false;

            for(size_t i=0; i<e; i++){
                if(start_nodes(i, 0) != plain_graph.starting_node_list[i] || end_nodes(i, 0) != plain_graph.ending_node_list[i]){
                    check_flag = false;
                    break;
                }
            }
            if(party_info.pIdx == 0){
                if(check_flag){
                    debug_info("\033[32m The secure graph is the same as the plaintext graph. \033[0m\n");
                }
                else{
                    debug_info("\033[31m Error: the secure graph is not the same as the plaintext graph. \033[0m\n");

                    // print the secure graph.
                    for(size_t i=0; i<e; i++){
                        debug_info("edge_list[" + std::to_string(i) + "]: " + std::to_string(start_nodes(i, 0)) + " " + std::to_string(end_nodes(i, 0)));
                        debug_info("plain_graph[" + std::to_string(i) + "]: " + std::to_string(plain_graph.starting_node_list[i]) + " " + std::to_string(plain_graph.ending_node_list[i]));
                    }
                }
            }

            return;
        }

        void check_graph(const std::string& meta_data_file, const std::string& file_prefix, aby3Info &party_info, bool multiparty){
            if(!multiparty){
                THROW_RUNTIME_ERROR("Not implemented!");
                return;
            }
            
            plainGraphList plain_graph_whole;
            plain_graph_whole.e = 0;
            size_t N_providers;
            std::ifstream meta(meta_data_file);
            meta >> v >> N_providers;

            for(size_t i=0; i<N_providers; i++){
                std::string party_meta_file = file_prefix + "_edge_list_meta_party-" + std::to_string(i) + ".txt";
                std::string party_data_file = file_prefix + "_edge_list_party-" + std::to_string(i) + ".txt";
                plainGraphList plain_graph(party_meta_file, party_data_file);
                plain_graph_whole.starting_node_list.insert(plain_graph_whole.starting_node_list.end(), plain_graph.starting_node_list.begin(), plain_graph.starting_node_list.end());
                plain_graph_whole.ending_node_list.insert(plain_graph_whole.ending_node_list.end(), plain_graph.ending_node_list.begin(), plain_graph.ending_node_list.end());
                plain_graph_whole.e += plain_graph.e;
            }

            check_graph(plain_graph_whole, party_info);

            return;
        }

};

size_t get_sending_bytes(aby3Info &party_info);

size_t get_receiving_bytes(aby3Info &party_info);

void clear_sending_reveiving_bytes(aby3Info &party_info);

aby3::sbMatrix get_target_node_mask(boolIndex target_start_node, aby3::sbMatrix& node_block, aby3Info &party_info);

aby3::sbMatrix get_target_edge_mask(boolIndex target_start_node, boolIndex target_end_node, aby3::sbMatrix& edge_block, aby3Info &party_info);

boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,
                         boolIndex logical_edge_block_index,
                         GraphQueryEngine &GQEngine);

aby3::si64Matrix outting_edge_count(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine);

aby3::si64Matrix outting_edge_range_statistics(boolIndex node_index, boolIndex logical_node_block_index, boolIndex upper_bound, GraphQueryEngine &GQEngine);

aby3::si64Matrix outting_edge_range_statistics(boolIndex node_index, boolIndex upper_bound, ListGraphQueryEngine &GQEngine);

aby3::sbMatrix outting_neighbors(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine);

aby3::sbMatrix outting_neighbors_sorted(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine);

boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,AdjGraphQueryEngine &GQEngine);

aby3::sbMatrix outting_edge_count(boolIndex boolIndex, AdjGraphQueryEngine &GQEngine);

aby3::si64Matrix outting_edge_count_arith(boolIndex boolIndex, AdjGraphQueryEngine &GQEngine);

aby3::sbMatrix outting_neighbors(boolIndex node_index,AdjGraphQueryEngine &GQEngine);

boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,ListGraphQueryEngine &GQEngine);

aby3::sbMatrix outting_edge_count(boolIndex boolIndex, ListGraphQueryEngine &GQEngine);

aby3::si64Matrix outting_edge_count_arith(boolIndex boolIndex, ListGraphQueryEngine &GQEngine);

aby3::sbMatrix outting_neighbors(boolIndex node_index, ListGraphQueryEngine &GQEngine);

aby3::sbMatrix outting_neighbors_sorted(boolIndex node_index, ListGraphQueryEngine &GQEngine);