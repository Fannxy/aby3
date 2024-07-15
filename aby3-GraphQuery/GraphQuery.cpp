#include "Graph.h"

static const size_t MAX_UNIT_SIZE = 1 << 25;

aby3::sbMatrix get_target_node_mask(boolIndex target_start_node, aby3::sbMatrix& node_block, aby3Info &party_info){

    // check whether node_block contains even rows.
    if(checkEven(node_block.rows()) == false){
        THROW_RUNTIME_ERROR("The node block should contain even rows.");
    }

    size_t edge_num = node_block.rows() / 2;

    // process the node block for the final result.
    aby3::sbMatrix expanded_node_mat(edge_num, BITSIZE);
    aby3::sbMatrix starting_node_mat(edge_num, BITSIZE);

    for(int i=0; i<2; i++){ // iterate for the two shares.
        std::copy(node_block.mShares[i].begin(), node_block.mShares[i].begin() + edge_num, expanded_node_mat.mShares[i].begin());
        for(size_t j=0; j<edge_num; j++){
            starting_node_mat.mShares[i](j, 0) = target_start_node.indexShares[i];
        }
    }

    // comparing the starting nodes with the target starting node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(party_info.pIdx, expanded_node_mat, starting_node_mat, eq_res, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    return eq_res;
}

aby3::sbMatrix get_unique_ending_nodes_in_edge_list(boolIndex target_start_node, aby3::sbMatrix& starting_nodes, aby3::sbMatrix& ending_nodes, aby3Info& party_info){

    size_t edge_num = starting_nodes.rows();    
    aby3::sbMatrix target_starting_nodes(edge_num, BITSIZE);
    for(size_t i=0; i<edge_num; i++){
        for(int j=0; j<2; j++){
            target_starting_nodes.mShares[j](i, 0) = target_start_node.indexShares[j];
        }
    }

    // get the mask indicating which elements are the target node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(party_info.pIdx, starting_nodes, target_starting_nodes, eq_res, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    // expand the length to match with the node id.
    aby3::sbMatrix aligned_target_node_mat(edge_num, BITSIZE);
    for(int i=0; i<2; i++){
        for(size_t j=0; j<edge_num; j++){
            aligned_target_node_mat.mShares[i](j, 0) = (eq_res.mShares[i](j, 0) == 1) ? -1 : 0;
        }
    }

    // multiply with the real starting_nodes to extract the target nodes.
    bool_cipher_and(party_info.pIdx, aligned_target_node_mat, ending_nodes, aligned_target_node_mat, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    // sort the nodes for grouping.
    std::vector<aby3::sbMatrix> masked_nodes(edge_num);
    for(size_t i=0; i<edge_num; i++){
        masked_nodes[i].resize(1, BITSIZE);
        for(int j=0; j<2; j++){
            masked_nodes[i].mShares[j](0, 0) = aligned_target_node_mat.mShares[j](i, 0);
        }
    }

    quick_sort(masked_nodes, party_info.pIdx, *(party_info.enc), *(party_info.eval), *(party_info.runtime), 32);

    // differente EQ to filter out the fliping nodes.
    aby3::sbMatrix left_shifted_nodes(edge_num-1, BITSIZE);
    aby3::sbMatrix right_shifted_nodes(edge_num-1, BITSIZE);

    for (size_t i = 0; i < edge_num-1; i++)
    {
        left_shifted_nodes.mShares[0](i, 0) = masked_nodes[i].mShares[0](0, 0);
        left_shifted_nodes.mShares[1](i, 0) = masked_nodes[i].mShares[1](0, 0);
        right_shifted_nodes.mShares[0](i, 0) = masked_nodes[i+1].mShares[0](0, 0);
        right_shifted_nodes.mShares[1](i, 0) = masked_nodes[i+1].mShares[1](0, 0);
    }

    bool_cipher_eq(party_info.pIdx, left_shifted_nodes, right_shifted_nodes, eq_res, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    eq_res.resize(eq_res.rows(), BITSIZE);
    for(size_t i=0; i<eq_res.rows(); i++){
        for(int j=0; j<2; j++){
            eq_res.mShares[j](i, 0) = (eq_res.mShares[j](i, 0) == 1) ? -1 : 0;
        }
    }
    eq_res.resize(eq_res.rows()+1, BITSIZE);
    boolIndex true_share(0, party_info.pIdx);
    eq_res.mShares[0](eq_res.rows()-1, 0) = true_share.indexShares[0];
    eq_res.mShares[1](eq_res.rows()-1, 0) = true_share.indexShares[1];

    bool_cipher_not(party_info.pIdx, eq_res, eq_res);

    for(size_t i=0; i<edge_num; i++){
        for(int j=0; j<2; j++){
            aligned_target_node_mat.mShares[j](i, 0) = masked_nodes[i].mShares[j](0, 0);
        }
    }

    aby3::sbMatrix filtered_nodes(edge_num, BITSIZE);
    bool_cipher_and(party_info.pIdx, eq_res, aligned_target_node_mat, filtered_nodes, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    // shuffle the masks and target nodes for privacy.
    efficient_shuffle(filtered_nodes, party_info.pIdx, filtered_nodes, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    return filtered_nodes;
}

aby3::sbMatrix get_unique_ending_nodes_in_sored_edge_list(boolIndex target_start_node, aby3::sbMatrix& starting_nodes, aby3::sbMatrix& ending_nodes, aby3Info& party_info){
    // the starting_nodes and ending_nodes are sorted by (starting, ending) pairs.
    size_t edge_num = starting_nodes.rows();
    aby3::sbMatrix target_starting_nodes(edge_num, BITSIZE);
    for(size_t i=0; i<edge_num; i++){
        for(int j=0; j<2; j++){
            target_starting_nodes.mShares[j](i, 0) = target_start_node.indexShares[j];
        }
    }

    // get the mask indicating which elements are the target node.
    aby3::sbMatrix eq_res(edge_num, BITSIZE);
    size_t round = (size_t)ceil(edge_num / (double)MAX_UNIT_SIZE);
    size_t last_len = edge_num - (round - 1) * MAX_UNIT_SIZE;

    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix start_node_(unit_len, BITSIZE);
        aby3::sbMatrix target_start_node_(unit_len, BITSIZE);

        std::memcpy(start_node_.mShares[0].data(), starting_nodes.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(starting_nodes.mShares[0](0, 0)));
        std::memcpy(start_node_.mShares[1].data(), starting_nodes.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(starting_nodes.mShares[1](0, 0)));
        std::memcpy(target_start_node_.mShares[0].data(), target_starting_nodes.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(target_starting_nodes.mShares[0](0, 0)));
        std::memcpy(target_start_node_.mShares[1].data(), target_starting_nodes.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(target_starting_nodes.mShares[1](0, 0)));

        aby3::sbMatrix eq_res_;

        bool_cipher_eq(party_info.pIdx, start_node_, target_start_node_, eq_res_, *(party_info.enc), *(party_info.eval), *(party_info.runtime));
        std::memcpy(eq_res.mShares[0].data() + i * MAX_UNIT_SIZE, eq_res_.mShares[0].data(), unit_len * sizeof(eq_res_.mShares[0](0, 0)));
        std::memcpy(eq_res.mShares[1].data() + i * MAX_UNIT_SIZE, eq_res_.mShares[1].data(), unit_len * sizeof(eq_res_.mShares[1](0, 0)));
    }
    // bool_cipher_eq(party_info.pIdx, starting_nodes, target_starting_nodes, eq_res, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    // expand the length to match with the node id.
    aby3::sbMatrix aligned_target_node_mat(edge_num, BITSIZE);
    for(int i=0; i<2; i++){
        for(size_t j=0; j<edge_num; j++){
            aligned_target_node_mat.mShares[i](j, 0) = (eq_res.mShares[i](j, 0) == 1) ? -1 : 0;
        }
    }

    // multiply with the real starting_nodes to extract the target nodes.

    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix ending_nodes_(unit_len, BITSIZE);
        aby3::sbMatrix aligned_target_node_mat_(unit_len, BITSIZE);

        std::memcpy(ending_nodes_.mShares[0].data(), ending_nodes.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(ending_nodes.mShares[0](0, 0)));
        std::memcpy(ending_nodes_.mShares[1].data(), ending_nodes.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(ending_nodes.mShares[1](0, 0)));
        std::memcpy(aligned_target_node_mat_.mShares[0].data(), aligned_target_node_mat.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(aligned_target_node_mat.mShares[0](0, 0)));
        std::memcpy(aligned_target_node_mat_.mShares[1].data(), aligned_target_node_mat.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(aligned_target_node_mat.mShares[1](0, 0)));

        bool_cipher_and(party_info.pIdx, aligned_target_node_mat_, ending_nodes_, aligned_target_node_mat_, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

        std::memcpy(aligned_target_node_mat.mShares[0].data() + i * MAX_UNIT_SIZE, aligned_target_node_mat_.mShares[0].data(), unit_len * sizeof(aligned_target_node_mat_.mShares[0](0, 0)));
        std::memcpy(aligned_target_node_mat.mShares[1].data() + i * MAX_UNIT_SIZE, aligned_target_node_mat_.mShares[1].data(), unit_len * sizeof(aligned_target_node_mat_.mShares[1](0, 0)));
    }

    // bool_cipher_and(party_info.pIdx, aligned_target_node_mat, ending_nodes, aligned_target_node_mat, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    // differente EQ to filter out the fliping nodes.
    aby3::sbMatrix left_shifted_nodes(edge_num-1, BITSIZE);
    aby3::sbMatrix right_shifted_nodes(edge_num-1, BITSIZE);

    for (size_t i = 0; i < edge_num-1; i++)
    {
        left_shifted_nodes.mShares[0](i, 0) = aligned_target_node_mat.mShares[0](i, 0);
        left_shifted_nodes.mShares[1](i, 0) = aligned_target_node_mat.mShares[1](i, 0);
        right_shifted_nodes.mShares[0](i, 0) = aligned_target_node_mat.mShares[0](i+1, 0);
        right_shifted_nodes.mShares[1](i, 0) = aligned_target_node_mat.mShares[1](i+1, 0);
    }

    // bool_cipher_eq(party_info.pIdx, left_shifted_nodes, right_shifted_nodes, eq_res, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len - 1 : MAX_UNIT_SIZE;
        aby3::sbMatrix left_shifted_nodes_(unit_len, BITSIZE);
        aby3::sbMatrix right_shifted_nodes_(unit_len, BITSIZE);
        aby3::sbMatrix eq_res_(unit_len, BITSIZE);

        std::memcpy(left_shifted_nodes_.mShares[0].data(), left_shifted_nodes.mShares[0].data() + i * MAX_UNIT_SIZE, (unit_len) * sizeof(left_shifted_nodes.mShares[0](0, 0)));
        std::memcpy(left_shifted_nodes_.mShares[1].data(), left_shifted_nodes.mShares[1].data() + i * MAX_UNIT_SIZE, (unit_len) * sizeof(left_shifted_nodes.mShares[1](0, 0)));
        std::memcpy(right_shifted_nodes_.mShares[0].data(), right_shifted_nodes.mShares[0].data() + i * MAX_UNIT_SIZE, (unit_len) * sizeof(right_shifted_nodes.mShares[0](0, 0)));
        std::memcpy(right_shifted_nodes_.mShares[1].data(), right_shifted_nodes.mShares[1].data() + i * MAX_UNIT_SIZE, (unit_len) * sizeof(right_shifted_nodes.mShares[1](0, 0)));

        bool_cipher_eq(party_info.pIdx, left_shifted_nodes_, right_shifted_nodes_, eq_res_, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

        std::memcpy(eq_res.mShares[0].data() + i * MAX_UNIT_SIZE, eq_res_.mShares[0].data(), (unit_len) * sizeof(eq_res_.mShares[0](0, 0)));
        std::memcpy(eq_res.mShares[1].data() + i * MAX_UNIT_SIZE, eq_res_.mShares[1].data(), (unit_len) * sizeof(eq_res_.mShares[1](0, 0)));
    }

    eq_res.resize(eq_res.rows(), BITSIZE);
    for(size_t i=0; i<eq_res.rows(); i++){
        for(int j=0; j<2; j++){
            eq_res.mShares[j](i, 0) = (eq_res.mShares[j](i, 0) == 1) ? -1 : 0;
        }
    }
    eq_res.resize(eq_res.rows()+1, BITSIZE);
    boolIndex true_share(0, party_info.pIdx);
    eq_res.mShares[0](eq_res.rows()-1, 0) = true_share.indexShares[0];
    eq_res.mShares[1](eq_res.rows()-1, 0) = true_share.indexShares[1];

    bool_cipher_not(party_info.pIdx, eq_res, eq_res);

    for(size_t i=0; i<edge_num; i++){
        for(int j=0; j<2; j++){
            aligned_target_node_mat.mShares[j](i, 0) = aligned_target_node_mat.mShares[j](i, 0);
        }
    }

    aby3::sbMatrix filtered_nodes(edge_num, BITSIZE);
    // bool_cipher_and(party_info.pIdx, eq_res, aligned_target_node_mat, filtered_nodes, *(party_info.enc), *(party_info.eval), *(party_info.runtime));
    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix eq_res_(unit_len, BITSIZE);
        aby3::sbMatrix aligned_target_node_mat_(unit_len, BITSIZE);

        std::memcpy(eq_res_.mShares[0].data(), eq_res.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res.mShares[0](0, 0)));
        std::memcpy(eq_res_.mShares[1].data(), eq_res.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res.mShares[1](0, 0)));
        std::memcpy(aligned_target_node_mat_.mShares[0].data(), aligned_target_node_mat.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(aligned_target_node_mat.mShares[0](0, 0)));
        std::memcpy(aligned_target_node_mat_.mShares[1].data(), aligned_target_node_mat.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(aligned_target_node_mat.mShares[1](0, 0)));

        aby3::sbMatrix filtered_nodes_(unit_len, BITSIZE);
        bool_cipher_and(party_info.pIdx, eq_res_, aligned_target_node_mat_, filtered_nodes_, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

        std::memcpy(filtered_nodes.mShares[0].data() + i * MAX_UNIT_SIZE, filtered_nodes_.mShares[0].data(), unit_len * sizeof(filtered_nodes_.mShares[0](0, 0)));
        std::memcpy(filtered_nodes.mShares[1].data() + i * MAX_UNIT_SIZE, filtered_nodes_.mShares[1].data(), unit_len * sizeof(filtered_nodes_.mShares[1](0, 0)));
    }

    // shuffle the masks and target nodes for privacy.
    #ifndef MPI_APP
    efficient_shuffle(filtered_nodes, party_info.pIdx, filtered_nodes, *(party_info.enc), *(party_info.eval), *(party_info.runtime));
    #endif

    return filtered_nodes;
}

// functions using GraphQueryEngine (based on Graph2D).
boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,
                         boolIndex logical_edge_block_index,
                         GraphQueryEngine &GQEngine) {
    // get the edge block from the GraphQueryEngine, which size is 2*l,
    // containing l beginning node and l ending node.
    aby3::sbMatrix edge_block =
        GQEngine.get_edge_block(logical_edge_block_index);

    // process the edge block for the final result.
    aby3::sbMatrix expanded_node_mat(GQEngine.graph->l * 2, BITSIZE);
    for (size_t i = 0; i < GQEngine.graph->l; i++) {
        expanded_node_mat.mShares[0](i, 0) = starting_node.indexShares[0];
        expanded_node_mat.mShares[1](i, 0) = starting_node.indexShares[1];
        expanded_node_mat.mShares[0](i + GQEngine.graph->l, 0) =
            ending_node.indexShares[0];
        expanded_node_mat.mShares[1](i + GQEngine.graph->l, 0) =
            ending_node.indexShares[1];
    }

    // comparing the starting nodes and ending nodes with the target starting
    // node and rge ending node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(GQEngine.party_info->pIdx, edge_block, expanded_node_mat,
                   eq_res, *(GQEngine.party_info->enc),
                   *(GQEngine.party_info->eval),
                   *(GQEngine.party_info->runtime));

    // check whether there exist pair (starting_node, ending_node) in the edge
    // block.
    aby3::sbMatrix starting_match(GQEngine.graph->l, 1),
        ending_match(GQEngine.graph->l, 1);

    for (int i = 0; i < 2; i++) {
        std::copy(eq_res.mShares[i].begin(),
                  eq_res.mShares[i].begin() + GQEngine.graph->l,
                  starting_match.mShares[i].begin());
        std::copy(eq_res.mShares[i].begin() + GQEngine.graph->l,
                  eq_res.mShares[i].end(), ending_match.mShares[i].begin());
    }

    aby3::sbMatrix matching_res(GQEngine.graph->l, 1);
    bool_cipher_and(GQEngine.party_info->pIdx, starting_match, ending_match,
                    matching_res, *(GQEngine.party_info->enc),
                    *(GQEngine.party_info->eval),
                    *(GQEngine.party_info->runtime));

    // aggregate for the final result, all using log-round OP computations.
    aby3::sbMatrix res(1, 1);
    bool_aggregation(GQEngine.party_info->pIdx, matching_res, res,
                        *(GQEngine.party_info->enc),
                        *(GQEngine.party_info->eval),
                        *(GQEngine.party_info->runtime), "OR");

    boolShare query_res;
    query_res.from_matrix(res.mShares[0](0, 0), res.mShares[1](0, 0));

    #ifdef MPI_APP
    // debug_info("mpi exge exist????");
    int role = GQEngine.party_info->pIdx;
    boolShare flag = query_res;
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int left_tasks = total_tasks;
    int start = (total_tasks + 1) / 2;
    while(rank < start && left_tasks > 1){
        int receive_target = rank + start; 
        if(receive_target < left_tasks){
            boolShare flag_recv;
            std::vector<char> charVec(2);
            MPI_Recv(charVec.data(), charVec.size(), MPI_CHAR, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // if(GQEngine.party_info->pIdx == 0) debug_info("my rank is " + std::to_string(rank) + "receive from " + std::to_string(receive_target));
            flag_recv = boolShare(charVec[0], charVec[1]);
            bool_cipher_or(GQEngine.party_info->pIdx, flag, flag_recv, flag, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
        }
        left_tasks = start;
        start = (start + 1) / 2;
    }
    if(rank >= start){
        int send_target = rank - start;
        std::vector<char> charVec = {flag.bshares[0], flag.bshares[1]};
        MPI_Send(charVec.data(), charVec.size(), MPI_CHAR, send_target, 0, MPI_COMM_WORLD);
        if(GQEngine.party_info->pIdx == 0) debug_info("my rank is " + std::to_string(rank) + "send to " + std::to_string(send_target));
    }
    #endif

    return query_res;
}

aby3::si64Matrix outting_edge_count(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine){

    // get the node block from the GraphQueryEngine, which size is b * 2l.
    aby3::sbMatrix node_block = GQEngine.get_node_edges(logical_node_block_index);

    // process the node block for the final result.
    aby3::sbMatrix expanded_node_mat(GQEngine.graph->l * GQEngine.graph->b, BITSIZE);
    aby3::sbMatrix starting_node_mat(GQEngine.graph->l * GQEngine.graph->b, BITSIZE);

    // TODO: node chunk data structure.
    for(int i=0; i<2; i++){
        std::copy(node_block.mShares[i].begin(), node_block.mShares[i].begin() + GQEngine.graph->l * GQEngine.graph->b, expanded_node_mat.mShares[i].begin());
        for(size_t j=0; j<GQEngine.graph->l * GQEngine.graph->b; j++){
            starting_node_mat.mShares[i](j, 0) = node_index.indexShares[i];
        }
    }

    // comparing the starting nodes with the target starting node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(GQEngine.party_info->pIdx, expanded_node_mat, starting_node_mat, eq_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    // trans2 A shares.
    aby3::si64Matrix eq_a_res(eq_res.rows(), 1);
    bool2arith(GQEngine.party_info->pIdx, eq_res, eq_a_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::si64Matrix res(1, 1);
    arith_aggregation(GQEngine.party_info->pIdx, eq_a_res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    #ifdef MPI_APP
    int role = GQEngine.party_info->pIdx;
    aby3::si64Matrix out_edges = res;
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int left_tasks = total_tasks;
    int start = (total_tasks + 1) / 2;
    while(rank < start && left_tasks > 1){
        int receive_target = rank + start; 
        if(receive_target < left_tasks){
            std::vector<int64_t> intVec(2);
            MPI_Recv(intVec.data(), intVec.size(), MPI_INT64_T, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            out_edges.mShares[0](0, 0) = out_edges.mShares[0](0, 0) + intVec[0];
            out_edges.mShares[1](0, 0) = out_edges.mShares[1](0, 0) + intVec[1];
        }
        left_tasks = start;
        start = (start + 1) / 2;
    }
    if(rank >= start){
        int send_target = rank - start;
        std::vector<int64_t> intVec = {out_edges.mShares[0](0, 0), out_edges.mShares[1](0, 0)};
        MPI_Send(intVec.data(), intVec.size(), MPI_INT64_T, send_target, 0, MPI_COMM_WORLD);
    }
    #endif

    return res;
}


aby3::si64Matrix outting_edge_range_statistics(boolIndex node_index, boolIndex logical_node_block_index, boolIndex upper_bound, GraphQueryEngine &GQEngine){

    // get the node block from the GraphQueryEngine, which size is b * 2l.
    aby3::sbMatrix node_block = GQEngine.get_node_edges_property(logical_node_block_index);

    // process the node block for the final result.
    aby3::sbMatrix expanded_node_mat(GQEngine.graph->l * GQEngine.graph->b, BITSIZE);
    aby3::sbMatrix starting_node_mat(GQEngine.graph->l * GQEngine.graph->b, BITSIZE);
    aby3::sbMatrix expanded_bound_mat(GQEngine.graph->l * GQEngine.graph->b, BITSIZE);
    aby3::sbMatrix property_mat(GQEngine.graph->l * GQEngine.graph->b, BITSIZE);

    // TODO: node chunk data structure.
    for(int i=0; i<2; i++){
        std::copy(node_block.mShares[i].begin(), node_block.mShares[i].begin() + GQEngine.graph->l * GQEngine.graph->b, expanded_node_mat.mShares[i].begin());
        #ifndef PROPERTY_COMPRESS
        std::copy(node_block.mShares[i].begin() + 2*GQEngine.graph->l * GQEngine.graph->b, node_block.mShares[i].end(), property_mat.mShares[i].begin());

        for(size_t j=0; j<GQEngine.graph->l * GQEngine.graph->b; j++){
            starting_node_mat.mShares[i](j, 0) = node_index.indexShares[i];
            expanded_bound_mat.mShares[i](j, 0) = upper_bound.indexShares[i];
        }
        #endif
        #ifdef PROPERTY_COMPRESS
        std::copy(node_block.mShares[i].begin() + GQEngine.graph->l * GQEngine.graph->b, node_block.mShares[i].end(), property_mat.mShares[i].begin());
        for(size_t j=0; j<GQEngine.graph->l * GQEngine.graph->b; j++){
            starting_node_mat.mShares[i](j, 0) = (node_index.indexShares[i] >> 30);
            expanded_bound_mat.mShares[i](j, 0) = upper_bound.indexShares[i];
        }
        #endif
    }

    // comparing the starting nodes with the target starting node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(GQEngine.party_info->pIdx, expanded_node_mat, starting_node_mat, eq_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::sbMatrix eq_res_bound;
    bool_cipher_lt(GQEngine.party_info->pIdx, expanded_bound_mat, property_mat, eq_res_bound, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    bool_cipher_and(GQEngine.party_info->pIdx, eq_res, eq_res_bound, eq_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    // trans2 A shares.
    aby3::si64Matrix eq_a_res(eq_res.rows(), 1);
    bool2arith(GQEngine.party_info->pIdx, eq_res, eq_a_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::si64Matrix res(1, 1);
    arith_aggregation(GQEngine.party_info->pIdx, eq_a_res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    #ifdef MPI_APP
    int role = GQEngine.party_info->pIdx;
    aby3::si64Matrix out_edges = res;
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int left_tasks = total_tasks;
    int start = (total_tasks + 1) / 2;
    while(rank < start && left_tasks > 1){
        int receive_target = rank + start; 
        if(receive_target < left_tasks){
            std::vector<int64_t> intVec(2);
            MPI_Recv(intVec.data(), intVec.size(), MPI_INT64_T, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            out_edges.mShares[0](0, 0) = out_edges.mShares[0](0, 0) + intVec[0];
            out_edges.mShares[1](0, 0) = out_edges.mShares[1](0, 0) + intVec[1];
        }
        left_tasks = start;
        start = (start + 1) / 2;
    }
    if(rank >= start){
        int send_target = rank - start;
        std::vector<int64_t> intVec = {out_edges.mShares[0](0, 0), out_edges.mShares[1](0, 0)};
        MPI_Send(intVec.data(), intVec.size(), MPI_INT64_T, send_target, 0, MPI_COMM_WORLD);
    }
    #endif

    return res;
}

aby3::si64Matrix outting_edge_range_statistics(boolIndex node_index, boolIndex upper_bound, ListGraphQueryEngine &GQEngine){
    aby3::sbMatrix expand_starting_node(GQEngine.e, BITSIZE);
    for(size_t i=0; i<GQEngine.e; i++){
        expand_starting_node.mShares[0](i, 0) = node_index.indexShares[0];
        expand_starting_node.mShares[1](i, 0) = node_index.indexShares[1];
    }

    aby3::sbMatrix expand_upper_bound(GQEngine.e, BITSIZE);
    for(size_t i=0; i<GQEngine.e; i++){
        expand_upper_bound.mShares[0](i, 0) = upper_bound.indexShares[0];
        expand_upper_bound.mShares[1](i, 0) = upper_bound.indexShares[1];
    }

    aby3::sbMatrix eq_res(GQEngine.e, BITSIZE);
    size_t round = (size_t)ceil(GQEngine.e / (double)MAX_UNIT_SIZE);
    size_t last_len = GQEngine.e - (round - 1) * MAX_UNIT_SIZE;

    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix start_(unit_len, BITSIZE);
        aby3::sbMatrix expand_(unit_len, BITSIZE);

        std::memcpy(start_.mShares[0].data(), GQEngine.starting_node_list.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.starting_node_list.mShares[0](0, 0)));
        std::memcpy(start_.mShares[1].data(), GQEngine.starting_node_list.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.starting_node_list.mShares[1](0, 0)));
        std::memcpy(expand_.mShares[0].data(), expand_starting_node.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_starting_node.mShares[0](0, 0)));
        std::memcpy(expand_.mShares[1].data(), expand_starting_node.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_starting_node.mShares[1](0, 0)));

        aby3::sbMatrix eq_res_;
        bool_cipher_eq(GQEngine.party_info->pIdx, start_, expand_, eq_res_, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

        aby3::sbMatrix property_(unit_len, BITSIZE);
        aby3::sbMatrix bound_(unit_len, BITSIZE);

        std::memcpy(property_.mShares[0].data(), GQEngine.edge_property_list.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.starting_node_list.mShares[0](0, 0)));
        std::memcpy(property_.mShares[1].data(), GQEngine.edge_property_list.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.starting_node_list.mShares[1](0, 0)));

        std::memcpy(bound_.mShares[0].data(), expand_upper_bound.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_starting_node.mShares[0](0, 0)));
        std::memcpy(bound_.mShares[1].data(), expand_upper_bound.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_starting_node.mShares[1](0, 0)));

        aby3::sbMatrix lt_res_;
        bool_cipher_lt(GQEngine.party_info->pIdx, property_, bound_, lt_res_, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

        bool_cipher_and(GQEngine.party_info->pIdx, eq_res_, lt_res_, eq_res_, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

        std::memcpy(eq_res.mShares[0].data() + i * MAX_UNIT_SIZE, eq_res_.mShares[0].data(), unit_len * sizeof(eq_res_.mShares[0](0, 0)));
        std::memcpy(eq_res.mShares[1].data() + i * MAX_UNIT_SIZE, eq_res_.mShares[1].data(), unit_len * sizeof(eq_res_.mShares[1](0, 0)));
    }

    aby3::si64Matrix eq_res_arith(eq_res.rows(), 1);

    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix data_(unit_len, BITSIZE);
        std::memcpy(data_.mShares[0].data(), eq_res.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res.mShares[0](0, 0)));
        std::memcpy(data_.mShares[1].data(), eq_res.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res.mShares[1](0, 0)));
        aby3::si64Matrix encrypted_data(unit_len, 1);
        bool2arith(GQEngine.party_info->pIdx, data_, encrypted_data, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

        std::memcpy(eq_res_arith.mShares[0].data() + i * MAX_UNIT_SIZE, encrypted_data.mShares[0].data(), unit_len * sizeof(encrypted_data.mShares[0](0, 0)));
        std::memcpy(eq_res_arith.mShares[1].data() + i * MAX_UNIT_SIZE, encrypted_data.mShares[1].data(), unit_len * sizeof(encrypted_data.mShares[1](0, 0)));
    }

    aby3::si64Matrix res(1, 1);
    aby3::i64Matrix res_one(1, 1);
    res_one(0, 0) = 1;
    if(GQEngine.party_info->pIdx == 0){
        GQEngine.party_info->enc->localIntMatrix(*(GQEngine.party_info->runtime), res_one, res).get();
    }
    else{
        GQEngine.party_info->enc->remoteIntMatrix(*(GQEngine.party_info->runtime), res).get();
    }
    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::si64Matrix res_(1, 1);
        aby3::si64Matrix data_(unit_len, 1);
        std::memcpy(data_.mShares[0].data(), eq_res_arith.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res_arith.mShares[0](0, 0)));
        std::memcpy(data_.mShares[1].data(), eq_res_arith.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res_arith.mShares[1](0, 0)));
        arith_aggregation(GQEngine.party_info->pIdx, data_, res_, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

        res.mShares[0](0, 0) = res.mShares[0](0, 0) + res_.mShares[0](0, 0);
        res.mShares[1](0, 0) = res.mShares[1](0, 0) + res_.mShares[1](0, 0);
    }

    #ifdef MPI_APP
    int role = GQEngine.party_info->pIdx;
    aby3::si64Matrix out_edges = res;
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int left_tasks = total_tasks;
    int start = (total_tasks + 1) / 2;
    while(rank < start && left_tasks > 1){
        int receive_target = rank + start; 
        if(receive_target < left_tasks){
            std::vector<int64_t> intVec(2);
            MPI_Recv(intVec.data(), intVec.size(), MPI_INT64_T, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            out_edges.mShares[0](0, 0) = out_edges.mShares[0](0, 0) + intVec[0];
            out_edges.mShares[1](0, 0) = out_edges.mShares[1](0, 0) + intVec[1];
        }
        left_tasks = start;
        start = (start + 1) / 2;
    }
    if(rank >= start){
        int send_target = rank - start;
        std::vector<int64_t> intVec = {out_edges.mShares[0](0, 0), out_edges.mShares[1](0, 0)};
        MPI_Send(intVec.data(), intVec.size(), MPI_INT64_T, send_target, 0, MPI_COMM_WORLD);
    }
    #endif


    return res;
}

aby3::sbMatrix outting_neighbors(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine){

    // get the node block from the GraphQueryEngine, which size is b * 2l.
    aby3::sbMatrix node_block = GQEngine.get_node_edges(logical_node_block_index);
    size_t edge_num = GQEngine.graph->l * GQEngine.graph->b;

    aby3::sbMatrix starting_nodes(edge_num, BITSIZE);
    aby3::sbMatrix ending_nodes(edge_num, BITSIZE);

    for(int i=0; i<2; i++){
        std::copy(node_block.mShares[i].begin(), node_block.mShares[i].begin() + edge_num, starting_nodes.mShares[i].begin());
        std::copy(node_block.mShares[i].begin() + edge_num, node_block.mShares[i].end(), ending_nodes.mShares[i].begin());
    }

    return get_unique_ending_nodes_in_edge_list(node_index, starting_nodes, ending_nodes, *(GQEngine.party_info));
}

aby3::sbMatrix outting_neighbors_sorted(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine){
    
        // get the node block from the GraphQueryEngine, which size is b * 2l.
        aby3::sbMatrix node_block = GQEngine.get_node_edges(logical_node_block_index);
    
        size_t edge_num = GQEngine.graph->l * GQEngine.graph->b;
    
        aby3::sbMatrix starting_nodes(edge_num, BITSIZE);
        aby3::sbMatrix ending_nodes(edge_num, BITSIZE);
    
        for(int i=0; i<2; i++){
            std::copy(node_block.mShares[i].begin(), node_block.mShares[i].begin() + edge_num, starting_nodes.mShares[i].begin());
            std::copy(node_block.mShares[i].begin() + edge_num, node_block.mShares[i].end(), ending_nodes.mShares[i].begin());
        }
    
        // return get_unique_ending_nodes_in_sored_edge_list(node_index, starting_nodes, ending_nodes, *(GQEngine.party_info));
        aby3::sbMatrix neighbors = get_unique_ending_nodes_in_sored_edge_list(node_index, starting_nodes, ending_nodes, *(GQEngine.party_info));

        // debug_info("in this function????? neighbors???");

        #ifdef MPI_APP
        int role = GQEngine.party_info->pIdx;
        aby3::Sh3Encryptor& enc = *(GQEngine.party_info->enc);
        aby3::Sh3Evaluator& eval = *(GQEngine.party_info->eval);
        aby3::Sh3Runtime& runtime = *(GQEngine.party_info->runtime);
        int total_tasks, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        // if(rank == 0) debug_info("in MPI neighbors get!!!");
        int left_tasks = total_tasks;
        int start = (total_tasks + 1) / 2;
        while(rank < start && left_tasks > 1){
            int receive_target = rank + start; 
            if(receive_target < left_tasks){
                // aby3::si64Matrix out_edges_recv;
                // std::vector<int64_t> intVec(out_edges.mShares[0].size());
                uint32_t recv_size;
                MPI_Recv(&recv_size, 1, MPI_UINT32_T, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::vector<int64_t> neighbors_vec(recv_size * 2);
                MPI_Recv(neighbors_vec.data(), neighbors_vec.size(), MPI_INT64_T, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // whether mask out the last element of the neighbors.
                aby3::sbMatrix neighbors_recv(recv_size, BITSIZE);
                std::memcpy(neighbors_recv.mShares[0].data(), neighbors_vec.data(), neighbors_recv.rows() * sizeof(int64_t));
                std::memcpy(neighbors_recv.mShares[1].data(), neighbors_vec.data() + neighbors.mShares[0].size(), neighbors_recv.rows() * sizeof(int64_t));

                aby3::sbMatrix last_element(recv_size, BITSIZE);
                std::fill_n(last_element.mShares[0].data(), last_element.mShares[0].size(), neighbors.mShares[0](neighbors.mShares[0].rows()-1, 0));
                std::fill_n(last_element.mShares[1].data(), last_element.mShares[1].size(), neighbors.mShares[1](neighbors.mShares[1].rows()-1, 0));

                aby3::sbMatrix mask(recv_size, 1);
                bool_cipher_eq(role, last_element, neighbors_recv, mask, enc, eval, runtime);

                aby3::sbMatrix final_mask(1, 1);
                bool_aggregation(role, mask, final_mask, enc, eval, runtime, "OR");
                boolShare final_flag((bool) final_mask.mShares[0](0, 0), (bool)final_mask.mShares[1](0, 0));
                aby3::sbMatrix zero_share(1, BITSIZE);
                zero_share.mShares[0].setZero();
                zero_share.mShares[1].setZero();
                aby3::sbMatrix last_single(1, BITSIZE);
                last_single.mShares[0](0, 0) = neighbors.mShares[0](neighbors.mShares[0].rows()-1, 0);
                last_single.mShares[1](0, 0) = neighbors.mShares[1](neighbors.mShares[1].rows()-1, 0);

                bool_cipher_selector(role, final_flag, zero_share, last_single, last_single, enc, eval, runtime);
                neighbors.mShares[0](neighbors.mShares[0].rows()-1, 0) = last_single.mShares[0](0, 0);
                neighbors.mShares[1](neighbors.mShares[1].rows()-1, 0) = last_single.mShares[1](0, 0);

                aby3::sbMatrix new_neighbors(neighbors.i64Size() + recv_size, BITSIZE);
                std::memcpy(new_neighbors.mShares[0].data(), neighbors.mShares[0].data(), neighbors.mShares[0].size() * sizeof(int64_t));
                std::memcpy(new_neighbors.mShares[1].data(), neighbors.mShares[1].data(), neighbors.mShares[1].size() * sizeof(int64_t));
                std::memcpy(new_neighbors.mShares[0].data() + neighbors.mShares[0].size(), neighbors_recv.mShares[0].data(), neighbors_recv.mShares[0].size() * sizeof(int64_t));
                std::memcpy(new_neighbors.mShares[1].data() + neighbors.mShares[1].size(), neighbors_recv.mShares[1].data(), neighbors_recv.mShares[1].size() * sizeof(int64_t));
                neighbors = new_neighbors;
            }
            left_tasks = start;
            start = (start + 1) / 2;
        }
        if(rank >= start){
            int send_target = rank - start;
            std::vector<int64_t> intVec(neighbors.mShares[0].size() * 2);
            std::memcpy(intVec.data(), neighbors.mShares[0].data(), neighbors.mShares[0].size() * sizeof(int64_t));
            std::memcpy(intVec.data() + neighbors.mShares[0].size(), neighbors.mShares[1].data(), neighbors.mShares[1].size() * sizeof(int64_t));
            uint32_t send_size = neighbors.mShares[0].size();
            MPI_Send(&send_size, 1, MPI_UINT32_T, send_target, 0, MPI_COMM_WORLD);
            MPI_Send(intVec.data(), intVec.size(), MPI_INT64_T, send_target, 0, MPI_COMM_WORLD);
        }
        if(rank == 0){
            efficient_shuffle(neighbors, role, neighbors, enc, eval, runtime);
        }
        #endif

        return neighbors;
}

// functions using AdjGraphQueryEngine.
boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,AdjGraphQueryEngine &GQEngine){
    boolIndex logical_edge_index = GQEngine.get_logical_edge_index(starting_node, ending_node);

    aby3::sbMatrix edge_info = GQEngine.get_target_edge(logical_edge_index);

    // whdther edge_info > 0 or not.
    aby3::sbMatrix res(1, 1);
    aby3::sbMatrix zero_share(edge_info.rows(), edge_info.bitCount());
    for(int i=0; i<edge_info.rows(); i++){
        zero_share.mShares[0](i, 0) = 0;
        zero_share.mShares[1](i, 0) = 0;
    }

    bool_cipher_lt(GQEngine.party_info->pIdx, zero_share, edge_info, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    boolShare query_res;
    query_res.from_matrix(res.mShares[0](0, 0), res.mShares[1](0, 0));

    return query_res;
}

aby3::sbMatrix outting_edge_count(boolIndex node_index, AdjGraphQueryEngine &GQEngine){
    aby3::sbMatrix res(1, 1);

    aby3::sbMatrix node_info = GQEngine.get_target_node(node_index);

    bool_aggregation(GQEngine.party_info->pIdx, node_info, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    return res;
}

aby3::si64Matrix outting_edge_count_arith(boolIndex node_index, AdjGraphQueryEngine &GQEngine){

    aby3::sbMatrix node_info = GQEngine.get_target_node(node_index);

    aby3::si64Matrix node_info_arith(node_info.i64Size(), 1);
    bool2arith(GQEngine.party_info->pIdx, node_info, node_info_arith, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::si64Matrix res(1, 1);
    arith_aggregation(GQEngine.party_info->pIdx, node_info_arith, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    return res;
}

aby3::sbMatrix outting_neighbors(boolIndex node_index, AdjGraphQueryEngine &GQEngine){
    aby3::sbMatrix res = GQEngine.get_target_node(node_index);

    // filterout the > 0 elements all to 1.
    aby3::sbMatrix zero_share(res.rows(), res.bitCount());
    for(int i=0; i<res.rows(); i++){
        zero_share.mShares[0](i, 0) = 0;
        zero_share.mShares[1](i, 0) = 0;
    }
    bool_cipher_lt(GQEngine.party_info->pIdx, zero_share, res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
    return res;
}

boolShare edge_existance(boolIndex starting_node, boolIndex ending_node, ListGraphQueryEngine &GQEngine){
    aby3::sbMatrix expand_starting_node(GQEngine.e, BITSIZE);
    aby3::sbMatrix expand_ending_node(GQEngine.e, BITSIZE);
    for(size_t i=0; i<GQEngine.e; i++){
        expand_starting_node.mShares[0](i, 0) = starting_node.indexShares[0];
        expand_starting_node.mShares[1](i, 0) = starting_node.indexShares[1];
        expand_ending_node.mShares[0](i, 0) = ending_node.indexShares[0];
        expand_ending_node.mShares[1](i, 0) = ending_node.indexShares[1];
    }

    aby3::sbMatrix full_comp_res(GQEngine.e, BITSIZE);

    size_t round = (size_t)ceil(GQEngine.e / (double)MAX_UNIT_SIZE);
    size_t last_len = GQEngine.e - (round - 1) * MAX_UNIT_SIZE;

    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix start_node_(unit_len, BITSIZE);
        aby3::sbMatrix end_node_(unit_len, BITSIZE);
        aby3::sbMatrix graph_node_start_(unit_len, BITSIZE);
        aby3::sbMatrix graph_node_end_(unit_len, BITSIZE);

        std::memcpy(start_node_.mShares[0].data(), expand_starting_node.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_starting_node.mShares[0](0, 0)));
        std::memcpy(start_node_.mShares[1].data(), expand_starting_node.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_starting_node.mShares[1](0, 0)));
        std::memcpy(end_node_.mShares[0].data(), expand_ending_node.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_ending_node.mShares[0](0, 0)));
        std::memcpy(end_node_.mShares[1].data(), expand_ending_node.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_ending_node.mShares[1](0, 0)));
        std::memcpy(graph_node_start_.mShares[0].data(), GQEngine.starting_node_list.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.starting_node_list.mShares[0](0, 0)));
        std::memcpy(graph_node_start_.mShares[1].data(), GQEngine.starting_node_list.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.starting_node_list.mShares[1](0, 0)));
        std::memcpy(graph_node_end_.mShares[0].data(), GQEngine.ending_node_list.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.ending_node_list.mShares[0](0, 0)));
        std::memcpy(graph_node_end_.mShares[1].data(), GQEngine.ending_node_list.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.ending_node_list.mShares[1](0, 0)));

        if(GQEngine.party_info->pIdx == 0){
            debug_info("before eq, round " + std::to_string(i));
        }

        aby3::sbMatrix eq_res_starts, eq_res_ends;
        bool_cipher_eq(GQEngine.party_info->pIdx, start_node_, graph_node_start_, eq_res_starts, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
        bool_cipher_eq(GQEngine.party_info->pIdx, end_node_, graph_node_end_, eq_res_ends, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

        aby3::sbMatrix full_comp_res_;
        bool_cipher_and(GQEngine.party_info->pIdx, eq_res_starts, eq_res_ends, full_comp_res_, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

        std::memcpy(full_comp_res.mShares[0].data() + i * MAX_UNIT_SIZE, full_comp_res_.mShares[0].data(), unit_len * sizeof(full_comp_res_.mShares[0](0, 0)));
        std::memcpy(full_comp_res.mShares[1].data() + i * MAX_UNIT_SIZE, full_comp_res_.mShares[1].data(), unit_len * sizeof(full_comp_res_.mShares[1](0, 0)));
    }

    if(GQEngine.party_info->pIdx == 0){
        debug_info("before aggregation");
        
    }
    // aby3::sbMatrix eq_res_starts, eq_res_ends;
    // bool_cipher_eq(GQEngine.party_info->pIdx, GQEngine.starting_node_list, expand_starting_node, eq_res_starts, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
    // bool_cipher_eq(GQEngine.party_info->pIdx, GQEngine.ending_node_list, expand_ending_node, eq_res_ends, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    // aby3::sbMatrix full_comp_res;
    // bool_cipher_and(GQEngine.party_info->pIdx, eq_res_starts, eq_res_ends, full_comp_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::sbMatrix res(1, 1);
    boolShare true_share(1, GQEngine.party_info->pIdx);
    res.mShares[0](0, 0) = true_share.bshares[0];
    res.mShares[1](0, 0) = true_share.bshares[1];

    for(size_t i=0; i<round; i++){
        if(GQEngine.party_info->pIdx == 0){
            debug_info("before aggregation, round " + std::to_string(i));
        }
        aby3::sbMatrix res_(1, 1);
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix data_(unit_len, 1);
        std::memcpy(data_.mShares[0].data(), full_comp_res.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(full_comp_res.mShares[0](0, 0)));
        std::memcpy(data_.mShares[1].data(), full_comp_res.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(full_comp_res.mShares[1](0, 0)));

        bool_aggregation(GQEngine.party_info->pIdx, data_, res_, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "OR");

        if(GQEngine.party_info->pIdx == 0){
            debug_info("after aggregation, before final - round " + std::to_string(i));
        }
        bool_cipher_or(GQEngine.party_info->pIdx, res, res_, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
    }

    // bool_aggregation(GQEngine.party_info->pIdx, full_comp_res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "OR");

    boolShare query_res;
    query_res.from_matrix(res.mShares[0](0, 0), res.mShares[1](0, 0));

    #ifdef MPI_APP
    int role = GQEngine.party_info->pIdx;
    boolShare flag = query_res;
    aby3::Sh3Encryptor& enc = *(GQEngine.party_info->enc);
    aby3::Sh3Evaluator& eval = *(GQEngine.party_info->eval);
    aby3::Sh3Runtime& runtime = *(GQEngine.party_info->runtime);
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int left_tasks = total_tasks;
    int start = (total_tasks + 1) / 2;
    while(rank < start && left_tasks > 1){
        int receive_target = rank + start; 
        if(receive_target < left_tasks){
            boolShare flag_recv;
            std::vector<char> charVec(2);
            MPI_Recv(charVec.data(), charVec.size(), MPI_CHAR, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            flag_recv = boolShare(charVec[0], charVec[1]);
            bool_cipher_or(role, flag, flag_recv, flag, enc, eval, runtime);
        }
        left_tasks = start;
        start = (start + 1) / 2;
    }
    if(rank >= start){
        int send_target = rank - start;
        std::vector<char> charVec = {flag.bshares[0], flag.bshares[1]};
        MPI_Send(charVec.data(), charVec.size(), MPI_CHAR, send_target, 0, MPI_COMM_WORLD);
    }
    #endif

    return query_res;
}

aby3::sbMatrix outting_edge_count(boolIndex boolIndex, ListGraphQueryEngine &GQEngine){
    aby3::sbMatrix expand_starting_node(GQEngine.e, BITSIZE);
    for(size_t i=0; i<GQEngine.e; i++){
        expand_starting_node.mShares[0](i, 0) = boolIndex.indexShares[0];
        expand_starting_node.mShares[1](i, 0) = boolIndex.indexShares[1];
    }

    aby3::sbMatrix eq_res;
    bool_cipher_eq(GQEngine.party_info->pIdx, GQEngine.starting_node_list, expand_starting_node, eq_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
    eq_res.resize(eq_res.rows(), BITSIZE);
    for(size_t i=0; i<eq_res.rows(); i++){
        for(int j=0; j<2; j++){
            eq_res.mShares[j](i, 0) = (eq_res.mShares[j](i, 0) == 1) ? 1 : 0;
        }
    }

    aby3::sbMatrix res(1, BITSIZE);
    bool_aggregation(GQEngine.party_info->pIdx, eq_res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    return res;
}

aby3::si64Matrix outting_edge_count_arith(boolIndex boolIndex, ListGraphQueryEngine &GQEngine){
    aby3::sbMatrix expand_starting_node(GQEngine.e, BITSIZE);
    for(size_t i=0; i<GQEngine.e; i++){
        expand_starting_node.mShares[0](i, 0) = boolIndex.indexShares[0];
        expand_starting_node.mShares[1](i, 0) = boolIndex.indexShares[1];
    }

    aby3::sbMatrix eq_res(GQEngine.e, BITSIZE);
    size_t round = (size_t)ceil(GQEngine.e / (double)MAX_UNIT_SIZE);
    size_t last_len = GQEngine.e - (round - 1) * MAX_UNIT_SIZE;
    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix start_(unit_len, BITSIZE);
        aby3::sbMatrix expand_(unit_len, BITSIZE);

        std::memcpy(start_.mShares[0].data(), GQEngine.starting_node_list.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.starting_node_list.mShares[0](0, 0)));
        std::memcpy(start_.mShares[1].data(), GQEngine.starting_node_list.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(GQEngine.starting_node_list.mShares[1](0, 0)));
        std::memcpy(expand_.mShares[0].data(), expand_starting_node.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_starting_node.mShares[0](0, 0)));
        std::memcpy(expand_.mShares[1].data(), expand_starting_node.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(expand_starting_node.mShares[1](0, 0)));

        aby3::sbMatrix eq_res_;
        bool_cipher_eq(GQEngine.party_info->pIdx, start_, expand_, eq_res_, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

        std::memcpy(eq_res.mShares[0].data() + i * MAX_UNIT_SIZE, eq_res_.mShares[0].data(), unit_len * sizeof(eq_res_.mShares[0](0, 0)));
        std::memcpy(eq_res.mShares[1].data() + i * MAX_UNIT_SIZE, eq_res_.mShares[1].data(), unit_len * sizeof(eq_res_.mShares[1](0, 0)));
    }
    // bool_cipher_eq(GQEngine.party_info->pIdx, GQEngine.starting_node_list, expand_starting_node, eq_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::si64Matrix eq_res_arith(eq_res.rows(), 1);
    // size_t round = (size_t)ceil(eq_res.rows() / (double)MAX_UNIT_SIZE);
    // size_t last_len = eq_res.rows() - (round - 1) * MAX_UNIT_SIZE;

    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::sbMatrix data_(unit_len, BITSIZE);
        std::memcpy(data_.mShares[0].data(), eq_res.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res.mShares[0](0, 0)));
        std::memcpy(data_.mShares[1].data(), eq_res.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res.mShares[1](0, 0)));
        aby3::si64Matrix encrypted_data(unit_len, 1);
        bool2arith(GQEngine.party_info->pIdx, data_, encrypted_data, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

        std::memcpy(eq_res_arith.mShares[0].data() + i * MAX_UNIT_SIZE, encrypted_data.mShares[0].data(), unit_len * sizeof(encrypted_data.mShares[0](0, 0)));
        std::memcpy(eq_res_arith.mShares[1].data() + i * MAX_UNIT_SIZE, encrypted_data.mShares[1].data(), unit_len * sizeof(encrypted_data.mShares[1](0, 0)));
    }

    aby3::si64Matrix res(1, 1);
    aby3::i64Matrix res_one(1, 1);
    res_one(0, 0) = 1;
    if(GQEngine.party_info->pIdx == 0){
        GQEngine.party_info->enc->localIntMatrix(*(GQEngine.party_info->runtime), res_one, res).get();
    }
    else{
        GQEngine.party_info->enc->remoteIntMatrix(*(GQEngine.party_info->runtime), res).get();
    }
    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_UNIT_SIZE;
        aby3::si64Matrix res_(1, 1);
        aby3::si64Matrix data_(unit_len, 1);
        std::memcpy(data_.mShares[0].data(), eq_res_arith.mShares[0].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res_arith.mShares[0](0, 0)));
        std::memcpy(data_.mShares[1].data(), eq_res_arith.mShares[1].data() + i * MAX_UNIT_SIZE, unit_len * sizeof(eq_res_arith.mShares[1](0, 0)));
        arith_aggregation(GQEngine.party_info->pIdx, data_, res_, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

        res.mShares[0](0, 0) = res.mShares[0](0, 0) + res_.mShares[0](0, 0);
        res.mShares[1](0, 0) = res.mShares[1](0, 0) + res_.mShares[1](0, 0);
    }
    // arith_aggregation(GQEngine.party_info->pIdx, eq_res_arith, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    #ifdef MPI_APP
    int role = GQEngine.party_info->pIdx;
    aby3::si64Matrix out_edges = res;
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int left_tasks = total_tasks;
    int start = (total_tasks + 1) / 2;
    while(rank < start && left_tasks > 1){
        int receive_target = rank + start; 
        if(receive_target < left_tasks){
            std::vector<int64_t> intVec(2);
            MPI_Recv(intVec.data(), intVec.size(), MPI_INT64_T, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            out_edges.mShares[0](0, 0) = out_edges.mShares[0](0, 0) + intVec[0];
            out_edges.mShares[1](0, 0) = out_edges.mShares[1](0, 0) + intVec[1];
        }
        left_tasks = start;
        start = (start + 1) / 2;
    }
    if(rank >= start){
        int send_target = rank - start;
        std::vector<int64_t> intVec = {out_edges.mShares[0](0, 0), out_edges.mShares[1](0, 0)};
        MPI_Send(intVec.data(), intVec.size(), MPI_INT64_T, send_target, 0, MPI_COMM_WORLD);
    }
    #endif

    return res;
}

aby3::sbMatrix outting_neighbors(boolIndex node_index, ListGraphQueryEngine &GQEngine){

    return get_unique_ending_nodes_in_edge_list(node_index, GQEngine.starting_node_list, GQEngine.ending_node_list, *(GQEngine.party_info));
}

aby3::sbMatrix outting_neighbors_sorted(boolIndex node_index, ListGraphQueryEngine &GQEngine){
    // return get_unique_ending_nodes_in_sored_edge_list(node_index, GQEngine.starting_node_list, GQEngine.ending_node_list, *(GQEngine.party_info));
    aby3::sbMatrix neighbors = get_unique_ending_nodes_in_sored_edge_list(node_index, GQEngine.starting_node_list, GQEngine.ending_node_list, *(GQEngine.party_info));
    #ifdef MPI_APP
    int role = GQEngine.party_info->pIdx;
    aby3::Sh3Encryptor& enc = *(GQEngine.party_info->enc);
    aby3::Sh3Evaluator& eval = *(GQEngine.party_info->eval);
    aby3::Sh3Runtime& runtime = *(GQEngine.party_info->runtime);
    int total_tasks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int left_tasks = total_tasks;
    int start = (total_tasks + 1) / 2;
    while(rank < start && left_tasks > 1){
        int receive_target = rank + start; 
        if(receive_target < left_tasks){
            uint32_t recv_size;
            MPI_Recv(&recv_size, 1, MPI_UINT32_T, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // aby3::si64Matrix out_edges_recv;
            // std::vector<int64_t> intVec(out_edges.mShares[0].size());
            std::vector<int64_t> neighbors_vec(recv_size * 2);
            MPI_Recv(neighbors_vec.data(), neighbors_vec.size(), MPI_INT64_T, receive_target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // whether mask out the last element of the neighbors.
            aby3::sbMatrix neighbors_recv(recv_size, BITSIZE);
            std::memcpy(neighbors_recv.mShares[0].data(), neighbors_vec.data(), neighbors_recv.rows() * sizeof(int64_t));
            std::memcpy(neighbors_recv.mShares[1].data(), neighbors_vec.data() + neighbors.mShares[0].size(), neighbors_recv.rows() * sizeof(int64_t));

            aby3::sbMatrix last_element(recv_size, BITSIZE);
            std::fill_n(last_element.mShares[0].data(), last_element.mShares[0].size(), neighbors.mShares[0](neighbors.mShares[0].rows()-1, 0));
            std::fill_n(last_element.mShares[1].data(), last_element.mShares[1].size(), neighbors.mShares[1](neighbors.mShares[1].rows()-1, 0));

            aby3::sbMatrix mask(recv_size, 1);
            bool_cipher_eq(role, last_element, neighbors_recv, mask, enc, eval, runtime);

            aby3::sbMatrix final_mask(1, 1);
            bool_aggregation(role, mask, final_mask, enc, eval, runtime, "OR");
            boolShare final_flag((bool) final_mask.mShares[0](0, 0), (bool)final_mask.mShares[1](0, 0));
            aby3::sbMatrix zero_share(1, BITSIZE);
            zero_share.mShares[0].setZero();
            zero_share.mShares[1].setZero();
            aby3::sbMatrix last_single(1, BITSIZE);
            last_single.mShares[0](0, 0) = neighbors.mShares[0](neighbors.mShares[0].rows()-1, 0);
            last_single.mShares[1](0, 0) = neighbors.mShares[1](neighbors.mShares[1].rows()-1, 0);

            bool_cipher_selector(role, final_flag, zero_share, last_single, last_single, enc, eval, runtime);
            neighbors.mShares[0](neighbors.mShares[0].rows()-1, 0) = last_single.mShares[0](0, 0);
            neighbors.mShares[1](neighbors.mShares[1].rows()-1, 0) = last_single.mShares[1](0, 0);

            aby3::sbMatrix new_neighbors(neighbors.i64Size() + recv_size, BITSIZE);
            std::memcpy(new_neighbors.mShares[0].data(), neighbors.mShares[0].data(), neighbors.mShares[0].size() * sizeof(int64_t));
            std::memcpy(new_neighbors.mShares[1].data(), neighbors.mShares[1].data(), neighbors.mShares[1].size() * sizeof(int64_t));
            std::memcpy(new_neighbors.mShares[0].data() + neighbors.mShares[0].size(), neighbors_recv.mShares[0].data(), neighbors_recv.mShares[0].size() * sizeof(int64_t));
            std::memcpy(new_neighbors.mShares[1].data() + neighbors.mShares[1].size(), neighbors_recv.mShares[1].data(), neighbors_recv.mShares[1].size() * sizeof(int64_t));
            neighbors = new_neighbors;
        }
        left_tasks = start;
        start = (left_tasks + 1) / 2;
    }
    if(rank >= start){
        int send_target = rank - start;
        std::vector<int64_t> intVec(neighbors.mShares[0].size() * 2);
        std::memcpy(intVec.data(), neighbors.mShares[0].data(), neighbors.mShares[0].size() * sizeof(int64_t));
        std::memcpy(intVec.data() + neighbors.mShares[0].size(), neighbors.mShares[1].data(), neighbors.mShares[1].size() * sizeof(int64_t));
        uint32_t send_size = neighbors.mShares[0].size();
        MPI_Send(&send_size, 1, MPI_UINT32_T, send_target, 0, MPI_COMM_WORLD);
        MPI_Send(intVec.data(), intVec.size(), MPI_INT64_T, send_target, 0, MPI_COMM_WORLD);
    }
    if(rank == 0){
        efficient_shuffle(neighbors, role, neighbors, enc, eval, runtime);
    }
    #endif

    return neighbors;
}