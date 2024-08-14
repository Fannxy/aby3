#include "Search.h"
#include <algorithm>
#include "../aby3-Basic/timer.h"

using namespace oc;
using namespace aby3;

#define FUNCTIONAL

static int PER_PARTY_MAX = 1 << 4;

int mcompBS(std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    // get the length.
    int n = keyset.size();
    if(n == 1) {
        boolShare le_res(1, pIdx);
        res.mShares[0](0, 0) = le_res.bshares[0];
        res.mShares[1](0, 0) = le_res.bshares[1];
        return 0;
    }
    // int bitsize = keyset[0].bitCount();
    int bitsize = key.bitCount();

    // comp key with keyset - can be compute in parallel.
    aby3::sbMatrix expand_key_set(n, bitsize);
    aby3::sbMatrix expand_key(n, bitsize);
    for(int i = 0; i < n; i++){
        expand_key_set.mShares[0](i, 0) = keyset[i].mData[0];
        expand_key_set.mShares[1](i, 0) = keyset[i].mData[1];
        expand_key.mShares[0](i, 0) = key.mShares[0](0, 0);
        expand_key.mShares[1](i, 0) = key.mShares[1](0, 0);
    }

    // compute comparison.
    aby3::sbMatrix le_res(n, 1);
    bool_cipher_lt(pIdx, expand_key_set, expand_key, le_res, enc, eval, runtime);
    bool_cipher_not(pIdx, le_res, le_res);

    // differential subtraction.
    res.resize(n, 1);
    res.mShares[0](0, 0) = le_res.mShares[0](0, 0);
    res.mShares[1](0, 0) = le_res.mShares[1](0, 0);

    aby3::sbMatrix diff_res(n-1, 1), left_le(n-1, 1), right_le(n-1, 1);

    for(int i = 0; i < n-1; i++){
        left_le.mShares[0](i, 0) = le_res.mShares[0](i, 0);
        left_le.mShares[1](i, 0) = le_res.mShares[1](i, 0);
        right_le.mShares[0](i, 0) = le_res.mShares[0](i+1, 0);
        right_le.mShares[1](i, 0) = le_res.mShares[1](i+1, 0);
    }

    bool_cipher_not(pIdx, left_le, left_le);
    bool_cipher_and(pIdx, left_le, right_le, diff_res, enc, eval, runtime);

    for(int i = 0; i < n-1; i++){
        res.mShares[0](i+1, 0) = diff_res.mShares[0](i, 0);
        res.mShares[1](i+1, 0) = diff_res.mShares[1](i, 0);
    }

    return 0;
}

int compBS(std::vector<aby3::si64> &data, std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    int n = data.size();
    aby3::sbMatrix comp_flag(n, 1);
    mcompBS(keyset, key, comp_flag, pIdx, enc, eval, runtime);

    // get the result.
    aby3::si64Matrix data_mat(n, 1);
    for(int i = 0; i < n; i++){
        data_mat.mShares[0](i, 0) = data[i].mData[0];
        data_mat.mShares[1](i, 0) = data[i].mData[1];
    }
    // compute the result.
    res.resize(1, 1);
    aby3::si64Matrix res_tmp(n, 1);
    arith_bool_mul(pIdx, data_mat, comp_flag, res_tmp, enc, eval, runtime);

    res.mShares[0](0, 0) = res_tmp.mShares[0](0, 0);
    res.mShares[1](0, 0) = res_tmp.mShares[1](0, 0);
    for(int i=1; i<n; i++){
        res.mShares[0](0, 0) += res_tmp.mShares[0](i, 0);
        res.mShares[1](0, 0) += res_tmp.mShares[1](i, 0);
    }
    return 0;
}

int distribute_mcompBSTGL(std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    int total_tasks, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);

    // currently only processor 0 have the data.
    int n, bitsize;
    if(rank == 0){
        n = keyset.size();
        // bitsize = keyset[0].bitCount();
        bitsize = key.bitCount();
    }
    // broadcast the sizes.
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // define the local task size.
    int local_n = n / total_tasks;
    if(rank == total_tasks - 1) local_n = n - (total_tasks - 1) * local_n;
    std::vector<int> local_counts(total_tasks);
    std::vector<int> local_displs(total_tasks);
    for(int i=0; i<total_tasks; i++){
        local_counts[i] = (i == total_tasks - 1) ? n - i * local_n : local_n;
        local_displs[i] = i * local_n;
    }

    // rank0 scatterv the keyset and bcast the key to all the tasks for parallel compare.
    aby3::sbMatrix local_keyset(local_n, bitsize);
    aby3::sbMatrix target_key(1, bitsize);
    if(rank == 0){
        target_key.mShares[0](0, 0) = key.mShares[0](0, 0);
        target_key.mShares[1](0, 0) = key.mShares[1](0, 0);
    }
    MPI_Bcast(target_key.mShares[0].data(), 1, MPI_INT64_T, 0, MPI_COMM_WORLD);   
    MPI_Bcast(target_key.mShares[1].data(), 1, MPI_INT64_T, 0, MPI_COMM_WORLD);

    if(rank == 0){
        aby3::sbMatrix full_keyset(n, bitsize);
        for(int i=0; i<n; i++){
            full_keyset.mShares[0](i, 0) = keyset[i].mData[0];
            full_keyset.mShares[1](i, 0) = keyset[i].mData[1];
        }
        MPI_Scatterv(full_keyset.mShares[0].data(), local_counts.data(), local_displs.data(), MPI_INT64_T, local_keyset.mShares[0].data(), local_n, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(full_keyset.mShares[1].data(), local_counts.data(), local_displs.data(), MPI_INT64_T, local_keyset.mShares[1].data(), local_n, MPI_INT64_T, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Scatterv(NULL, local_counts.data(), local_displs.data(), MPI_INT64_T, local_keyset.mShares[0].data(), local_n, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(NULL, local_counts.data(), local_displs.data(), MPI_INT64_T, local_keyset.mShares[1].data(), local_n, MPI_INT64_T, 0, MPI_COMM_WORLD);
    }

    // each task compute the LE in parallel.
    aby3::sbMatrix expand_key(local_n, bitsize);
    for(int i=0; i<local_n; i++){
        expand_key.mShares[0](i, 0) = target_key.mShares[0](0, 0);
        expand_key.mShares[1](i, 0) = target_key.mShares[1](0, 0);
    }
    aby3::sbMatrix le_res(local_n, 1);
    bool_cipher_lt(pIdx, local_keyset, expand_key, le_res, enc, eval, runtime);
    bool_cipher_not(pIdx, le_res, le_res);

    // sync and then compute the differential computation, by using shared memory.
    int64_t *share1;
    int64_t *share2;
    MPI_Win win_s1, win_s2;
    if(rank == 0){
        MPI_Win_allocate_shared(n * sizeof(int64_t), sizeof(int64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &share1, &win_s1);
        MPI_Win_allocate_shared(n * sizeof(int64_t), sizeof(int64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &share2, &win_s2);
    }
    else{
        MPI_Aint size;
        int disp_unit;
        MPI_Win_allocate_shared(0, sizeof(int64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &share1, &win_s1);
        MPI_Win_allocate_shared(0, sizeof(int64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &share2, &win_s2);
        MPI_Win_shared_query(win_s1, 0, &size, &disp_unit, &share1);
        MPI_Win_shared_query(win_s2, 0, &size, &disp_unit, &share2);
    }

    // sync and rank0 compute the result.
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Win_lock_all(MPI_MODE_NOCHECK, win_s1);
    for(int i=0; i<local_n; i++){
        share1[rank * (n / total_tasks) + i] = le_res.mShares[0](i, 0);
    }
    MPI_Win_unlock_all(win_s1);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, win_s2);
    for(int i=0; i<local_n; i++){
        share2[rank * (n / total_tasks) + i] = le_res.mShares[1](i, 0);
    }
    MPI_Win_unlock_all(win_s2);
    MPI_Barrier(MPI_COMM_WORLD);

    // in parrallel compute the differential computation.
    aby3::sbMatrix diff_res, left_le, right_le;
    if(rank == 0){
        diff_res.resize(local_n-1, 1);
        left_le.resize(local_n-1, 1);
        right_le.resize(local_n-1, 1);
        for(int i=1; i<local_n; i++){
            left_le.mShares[0](i-1, 0) = share1[i-1];
            left_le.mShares[1](i-1, 0) = share2[i-1];
            right_le.mShares[0](i-1, 0) = share1[i];
            right_le.mShares[1](i-1, 0) = share2[i];
        }
    }
    else if (rank > 0 && rank < total_tasks - 1)
    {
        diff_res.resize(local_n, 1);
        left_le.resize(local_n, 1);
        right_le.resize(local_n, 1);

        int start = rank * (n / total_tasks);
        for(int i=0; i<local_n; i++){
            left_le.mShares[0](i, 0) = share1[start + i - 1];   
            left_le.mShares[1](i, 0) = share2[start + i - 1];
            right_le.mShares[0](i, 0) = share1[start + i];
            right_le.mShares[1](i, 0) = share2[start + i];
        }
    }
    else{
        diff_res.resize(local_n, 1);
        left_le.resize(local_n, 1);
        right_le.resize(local_n, 1);

        int start = rank * (n / total_tasks);
        for(int i=0; i<local_n; i++){
            left_le.mShares[0](i, 0) = share1[start + i - 1];
            left_le.mShares[1](i, 0) = share2[start + i - 1];
            right_le.mShares[0](i, 0) = share1[start + i];
            right_le.mShares[1](i, 0) = share2[start + i];
        }
    }

    bool_cipher_not(pIdx, left_le, left_le);
    bool_cipher_and(pIdx, left_le, right_le, diff_res, enc, eval, runtime);

    // rank0 and the last rank will append the starting and ending result.
    res.resize(local_n, 1);
    if(rank == 0){
        res.mShares[0](0, 0) = share1[0];
        res.mShares[1](0, 0) = share2[0];
        for(int i=1; i<local_n; i++){
            res.mShares[0](i, 0) = diff_res.mShares[0](i-1, 0);
            res.mShares[1](i, 0) = diff_res.mShares[1](i-1, 0);
        }
    }
    else if (rank > 0 && rank < (total_tasks - 1))
    {
        std::copy(diff_res.mShares[0].data(), diff_res.mShares[0].data() + local_n, res.mShares[0].data());
        std::copy(diff_res.mShares[1].data(), diff_res.mShares[1].data() + local_n, res.mShares[1].data());
    }
    else{
        for(int i=0; i<local_n; i++){
            res.mShares[0](i, 0) = diff_res.mShares[0](i, 0);
            res.mShares[1](i, 0) = diff_res.mShares[1](i, 0);
        }
    }

    return 0;

}

int mcompBSTGLAB(std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, aby3::si64Matrix& res_a, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    int total_tasks, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);

    // currently only processor 0 have the data.
    int n, bitsize;
    if(rank == 0){
        n = keyset.size();
        // bitsize = keyset[0].bitCount();
        bitsize = key.bitCount();
    }
    // broadcast the sizes.
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // // define the local task size.
    int local_n = n / total_tasks;
    if(rank == total_tasks - 1) local_n = n - (total_tasks - 1) * local_n;
    std::vector<int> local_counts(total_tasks);
    std::vector<int> local_displs(total_tasks);
    for(int i=0; i<total_tasks; i++){
        local_counts[i] = (i == total_tasks - 1) ? n - i * local_n : local_n;
        local_displs[i] = i * local_n;
    }

    // call the distribute function.
    aby3::sbMatrix final_diff_res(local_n, 1);
    distribute_mcompBSTGL(keyset, key, final_diff_res, pIdx, enc, eval, runtime);

    // converse to the arith result.
    aby3::si64Matrix final_diff_res_a(local_n, 1);
    aby3::si64Matrix ones(local_n, 1);
    init_ones(pIdx, enc, runtime, ones, local_n);
    arith_bool_mul(pIdx,  ones, final_diff_res, final_diff_res_a, enc, eval, runtime);

    // sync and then gather the result -> all the ranks synchronize the result to rank0.
    if(rank == 0){
        res.resize(n, 1);
        res_a.resize(n, 1);
        std::copy(final_diff_res.mShares[0].data(), final_diff_res.mShares[0].data() + local_n, res.mShares[0].data());
        std::copy(final_diff_res.mShares[1].data(), final_diff_res.mShares[1].data() + local_n, res.mShares[1].data());
        std::copy(final_diff_res_a.mShares[0].data(), final_diff_res_a.mShares[0].data() + local_n, res_a.mShares[0].data());
        std::copy(final_diff_res_a.mShares[1].data(), final_diff_res_a.mShares[1].data() + local_n, res_a.mShares[1].data());

        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, res.mShares[0].data(), local_counts.data(), local_displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, res.mShares[1].data(), local_counts.data(), local_displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, res_a.mShares[0].data(), local_counts.data(), local_displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, res_a.mShares[1].data(), local_counts.data(), local_displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Gatherv(final_diff_res.mShares[0].data(), local_n, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD); 
        MPI_Gatherv(final_diff_res.mShares[1].data(), local_n, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD); 
        MPI_Gatherv(final_diff_res_a.mShares[0].data(), local_n, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(final_diff_res_a.mShares[1].data(), local_n, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
    }

    return 0;
}

int mcompBSTGL(std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    int total_tasks, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);

    // currently only processor 0 have the data.
    int n, bitsize;
    if(rank == 0){
        n = keyset.size();
        // bitsize = keyset[0].bitCount();
        bitsize = key.bitCount();
    }
    // broadcast the sizes.
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // // define the local task size.
    int local_n = n / total_tasks;
    if(rank == total_tasks - 1) local_n = n - (total_tasks - 1) * local_n;
    std::vector<int> local_counts(total_tasks);
    std::vector<int> local_displs(total_tasks);
    for(int i=0; i<total_tasks; i++){
        local_counts[i] = (i == total_tasks - 1) ? n - i * local_n : local_n;
        local_displs[i] = i * local_n;
    }

    // call the distribute function.
    aby3::sbMatrix final_diff_res(local_n, 1);
    distribute_mcompBSTGL(keyset, key, final_diff_res, pIdx, enc, eval, runtime);

    // sync and then gather the result -> all the ranks synchronize the result to rank0.
    if(rank == 0){
        res.resize(n, 1);
        std::copy(final_diff_res.mShares[0].data(), final_diff_res.mShares[0].data() + local_n, res.mShares[0].data());
        std::copy(final_diff_res.mShares[1].data(), final_diff_res.mShares[1].data() + local_n, res.mShares[1].data());

        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, res.mShares[0].data(), local_counts.data(), local_displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, res.mShares[1].data(), local_counts.data(), local_displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Gatherv(final_diff_res.mShares[0].data(), local_n, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD); 
        MPI_Gatherv(final_diff_res.mShares[1].data(), local_n, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD); 
    }

    return 0;
}

int compBSTGL(std::vector<aby3::si64> &data, std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    // get the configs.
    int total_tasks, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);

    // first call the mcompBSTGLAB and only proc0 got the result.
    aby3::sbMatrix m_res;
    aby3::si64Matrix m_res_a;
    mcompBSTGLAB(keyset, key, m_res, m_res_a, pIdx, enc, eval, runtime);
    
    // then task0 compute the result.
    // 1) task0 organize the data into the format.
    if(rank == 0){
        std::vector<aby3::si64Matrix> data_mat(1);
        std::vector<aby3::si64Matrix> comp_mat(1);
        data_mat[0].resize(data.size(), 1);
        comp_mat[0].resize(data.size(), 1);
        for(int i=0; i<data.size(); i++){
            data_mat[0].mShares[0](i, 0) = data[i].mData[0];
            data_mat[0].mShares[1](i, 0) = data[i].mData[1];
            comp_mat[0].mShares[0](i, 0) = m_res_a.mShares[0](i, 0);
            comp_mat[0].mShares[1](i, 0) = m_res_a.mShares[1](i, 0);
        }

        // 2) task0 compute the result.
        res.resize(1, 1);
        constant_sint_dot(pIdx, data_mat, comp_mat, res, enc, eval, runtime);
    }

    return 0;
}

// the pure NDSS implementation.
int subHBS_with_TGL(std::vector<aby3::si64> &data, std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, int alpha, int threshold){
    
    // get the configs.
    int total_tasks, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks);
    int n, bitsize;
    if(rank == 0){
        n = data.size();
        // bitsize = keyset[0].bitCount();
        bitsize = key.bitCount();
    }
    // broadcast the sizes.
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bitsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(n < threshold){
        // call the compBSTGL function.
        if(n < total_tasks * PER_PARTY_MAX){
            if(rank == 0){
                return compBS(data, keyset, key, res, pIdx, enc, eval, runtime);
            }
            else{
                return 0;
            }
        }
        // otherwise call the ptaBS.
        return compBSTGL(data, keyset, key, res, pIdx, enc, eval, runtime);
    }

    int beta = n / alpha;
    if(alpha * beta != n) THROW_RUNTIME_ERROR("The size of data should be alpha * beta, while now alpha = " + std::to_string(alpha) + " - beta = " + std::to_string(beta) + " - n = " + std::to_string(n) + " However alpha * beta = " + std::to_string(alpha * beta));

    // recursive call the subHBS.
    std::vector<sb64> keyset_upper;
    aby3::sbMatrix res_upper(alpha, 1);
    aby3::si64Matrix res_upper_a(alpha, 1);
    if(rank == 0){
        // construct the upper-level tree.
        for(size_t i=0; i<alpha; i++){
            keyset_upper.push_back(keyset[i*beta + (beta-1)]);
        }
    }
    mcompBSTGLAB(keyset_upper, key, res_upper, res_upper_a, pIdx, enc, eval, runtime);

    // the p0 task perform the following computation till recursive.
    // parallelize the computation by beta. 
    // bcast res_upper to all parties.
    MPI_Bcast(res_upper.mShares[0].data(), (int)res_upper.i64Size(), MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(res_upper.mShares[1].data(), (int)res_upper.i64Size(), MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(res_upper_a.mShares[0].data(), alpha, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(res_upper_a.mShares[1].data(), alpha, MPI_INT64_T, 0, MPI_COMM_WORLD);

    int local_beta = beta / total_tasks;
    if(rank == total_tasks - 1){
        local_beta = beta - local_beta * (total_tasks - 1);
    }
    std::vector<int> sendcounts(total_tasks);
    std::vector<int> displs(total_tasks);
    std::vector<int> sendcounts_alpha(total_tasks);
    std::vector<int> displs_alpha(total_tasks);
    if(rank == 0){
        for(int i=0; i<total_tasks - 1; i++){
            sendcounts[i] = local_beta;
            displs[i] = i * local_beta;
            sendcounts_alpha[i] = local_beta * alpha;
            displs_alpha[i] = i * local_beta * alpha;
        }
        sendcounts[total_tasks - 1] = (beta - (total_tasks - 1) * local_beta);
        sendcounts_alpha[total_tasks - 1] = (beta - (total_tasks - 1) * local_beta) * alpha;
        displs[total_tasks - 1] = (total_tasks - 1) * local_beta;
        displs_alpha[total_tasks - 1] = (total_tasks - 1) * local_beta * alpha;
    }

    si64Matrix data_to_be_search(local_beta*alpha, 1);
    sbMatrix keyset_to_be_search(local_beta*alpha, bitsize);
    if(rank == 0){
        si64Matrix all_data(beta*alpha, 1);
        for(size_t i=0; i<beta; i++){
            for(size_t j=0; j<alpha; j++){
                all_data.mShares[0](i*alpha + j, 0) = data[j*beta + i].mData[0];
                all_data.mShares[1](i*alpha + j, 0) = data[j*beta + i].mData[1];
            }
        }

        // scatter the data.
        MPI_Scatterv(all_data.mShares[0].data(), sendcounts_alpha.data(), displs_alpha.data(), MPI_INT64_T, data_to_be_search.mShares[0].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(all_data.mShares[1].data(), sendcounts_alpha.data(), displs_alpha.data(), MPI_INT64_T, data_to_be_search.mShares[1].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);

        // data.resize(0);
        all_data.resize(0, 0);

        sbMatrix all_keyset(beta*alpha, bitsize);
        for(size_t i=0; i<beta; i++){
            for(size_t j=0; j<alpha; j++){
                all_keyset.mShares[0](i*alpha + j, 0) = keyset[j*beta + i].mData[0];
                all_keyset.mShares[1](i*alpha + j, 0) = keyset[j*beta + i].mData[1];
            }
        }
        // scatter the data.
        MPI_Scatterv(all_keyset.mShares[0].data(), sendcounts_alpha.data(), displs_alpha.data(), MPI_INT64_T, keyset_to_be_search.mShares[0].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(all_keyset.mShares[1].data(), sendcounts_alpha.data(), displs_alpha.data(), MPI_INT64_T, keyset_to_be_search.mShares[1].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);

        // keyset.resize(0);
        all_keyset.resize(0, 0);
    }
    else{
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT64_T, data_to_be_search.mShares[0].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT64_T, data_to_be_search.mShares[1].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT64_T, keyset_to_be_search.mShares[0].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT64_T, keyset_to_be_search.mShares[1].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
    }

    // each task independetly compute the result.
    std::vector<si64Matrix> data_lowers(local_beta);
    std::vector<si64Matrix> tag_lowers(local_beta);
    si64Matrix data_res(local_beta, 1);
    for(size_t i=0; i<local_beta; i++){
        tag_lowers[i] = res_upper_a;
        data_lowers[i].resize(alpha, 1);
        for(size_t j=0; j<alpha; j++){
            data_lowers[i].mShares[0](j, 0) = data_to_be_search.mShares[0](i*alpha + j, 0);
            data_lowers[i].mShares[1](j, 0) = data_to_be_search.mShares[1](i*alpha + j, 0);
        }
    }
    constant_sint_dot(pIdx, data_lowers, tag_lowers, data_res, enc, eval, runtime);

    // constant inner product for key.
    std::vector<sbMatrix> key_lowers(local_beta);
    std::vector<sbMatrix> keytag_lowers(local_beta);
    sbMatrix key_res(local_beta, bitsize);
    for(size_t i=0; i<local_beta; i++){
        key_lowers[i].resize(alpha, bitsize);
        keytag_lowers[i].resize(alpha, bitsize);
        for(size_t j=0; j<alpha; j++){
            key_lowers[i].mShares[0](j, 0) = keyset_to_be_search.mShares[0](i*alpha + j, 0);
            key_lowers[i].mShares[1](j, 0) = keyset_to_be_search.mShares[1](i*alpha + j, 0);
            keytag_lowers[i].mShares[0](j, 0) = (res_upper.mShares[0](j, 0) == 1) ? -1 : 0;
            keytag_lowers[i].mShares[1](j, 0) = (res_upper.mShares[1](j, 0) == 1) ? -1 : 0;
        }
    }
    constant_bool_dot(pIdx, key_lowers, keytag_lowers, key_res, enc, eval, runtime);

    // gather the result and recurse
    sbMatrix key_res_all;
    si64Matrix data_res_all;
    std::vector<si64> data_to_be_search_all;
    std::vector<sb64> keyset_to_be_search_all;
    if(rank == 0){
        key_res_all.resize(beta, bitsize);
        data_res_all.resize(beta, 1);

        std::copy(key_res.mShares[0].data(), key_res.mShares[0].data() + local_beta, key_res_all.mShares[0].data());
        std::copy(key_res.mShares[1].data(), key_res.mShares[1].data() + local_beta, key_res_all.mShares[1].data());
        std::copy(data_res.mShares[0].data(), data_res.mShares[0].data() + local_beta, data_res_all.mShares[0].data());
        std::copy(data_res.mShares[1].data(), data_res.mShares[1].data() + local_beta, data_res_all.mShares[1].data());

        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, key_res_all.mShares[0].data(), sendcounts.data(), displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, key_res_all.mShares[1].data(), sendcounts.data(), displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, data_res_all.mShares[0].data(), sendcounts.data(), displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, data_res_all.mShares[1].data(), sendcounts.data(), displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);

        data_to_be_search_all.resize(beta);
        keyset_to_be_search_all.resize(beta);
        for(size_t i=0; i<beta; i++){
            // si64Matrix data_tmp(1, 1);
            // sbMatrix key_tmp(1, bitsize);
            // data_tmp.mShares[0](0, 0) = data_res_all.mShares[0](i, 0);
            // data_tmp.mShares[1](0, 0) = data_res_all.mShares[1](i, 0);
            // key_tmp.mShares[0](0, 0) = key_res_all.mShares[0](i, 0);
            // key_tmp.mShares[1](0, 0) = key_res_all.mShares[1](i, 0);
            // data_to_be_search_all.push_back(data_tmp);
            // keyset_to_be_search_all.push_back(key_tmp);
            data_to_be_search_all[i].mData[0] = data_res_all.mShares[0](i, 0);
            data_to_be_search_all[i].mData[1] = data_res_all.mShares[1](i, 0);
            keyset_to_be_search_all[i].mData[0] = key_res_all.mShares[0](i, 0);
            keyset_to_be_search_all[i].mData[1] = key_res_all.mShares[1](i, 0);
        }
    }
    else{
        MPI_Gatherv(key_res.mShares[0].data(), local_beta, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(key_res.mShares[1].data(), local_beta, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(data_res.mShares[0].data(), local_beta, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(data_res.mShares[1].data(), local_beta, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
    }

    return subHBS_with_TGL(data_to_be_search_all, keyset_to_be_search_all, key, res, pIdx, enc, eval, runtime, alpha, threshold);
}

int mtagBS(std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    int n = keyset.size();
    // int bitsize = keyset[0].bitCount();
    int bitsize = key.bitCount();

    if(checkPowerOfTwo(n) == false){
        THROW_RUNTIME_ERROR("Currently, the size of data should be power of 2!");
    }
    int logn = log2(n);

    // get the mid point.
    int p = n/2 - 1;
    sbMatrix c(1, 1);
    boolShare q(1, pIdx);
    sbMatrix _key(1, bitsize);
    _key.mShares[0](0, 0) = keyset[p].mData[0];
    _key.mShares[1](0, 0) = keyset[p].mData[1];
    bool_cipher_lt(pIdx, _key, key, c, enc, eval, runtime);
    bool_cipher_not(pIdx, c, c);

    aby3::sbMatrix tag(1, 1);
    aby3::sbMatrix tag_next;
    tag.mShares[0](0, 0) = q.bshares[0]; tag.mShares[1](0, 0) = q.bshares[1];

    // log-round binary search.
    for(size_t i=0; i<logn; i++){
        int tag_size = (1 << i);

        aby3::sbMatrix keymat(tag_size, bitsize);
        aby3::sbMatrix keytag(tag_size, bitsize);
        aby3::sbMatrix res_key_tmp(tag_size, bitsize);
        for(size_t j=0; j<tag_size; j++){
            int tag_ind = (((2*(j+1) - 1) * n) / (2*tag_size)) - 1;
            keymat.mShares[0](j, 0) = keyset[tag_ind].mData[0];
            keymat.mShares[1](j, 0) = keyset[tag_ind].mData[1];
            keytag.mShares[0](j, 0) = (tag.mShares[0](j, 0) == 1) ? -1 : 0;
            keytag.mShares[1](j, 0) = (tag.mShares[1](j, 0) == 1) ? -1 : 0;
        }

        aby3::sbMatrix akey(1, bitsize);

        std::vector<sbMatrix> enc_key_mat = {keymat};
        std::vector<sbMatrix> enc_key_tag = {keytag};
        constant_bool_dot(pIdx, enc_key_mat, enc_key_tag, akey, enc, eval, runtime);

        bool_cipher_lt(pIdx, akey, key, c, enc, eval, runtime);
        bool_cipher_not(pIdx, c, c);

        aby3::sbMatrix notC(1, 1);
        bool_cipher_not(pIdx, c, notC);

        // compute the tag for the next layer.
        int next_size = (1 << (i+1));
        tag_next.resize(next_size, 1);
        aby3::sbMatrix tag_expand(next_size, 1);
        aby3::sbMatrix comp_and_inv_comp(next_size, 1);

        for(size_t j=0; j<next_size; j++){
            tag_expand.mShares[0](j, 0) = tag.mShares[0](j/2, 0);
            tag_expand.mShares[1](j, 0) = tag.mShares[1](j/2, 0);
            if(j%2 == 0){
                comp_and_inv_comp.mShares[0](j, 0) = c.mShares[0](0, 0);
                comp_and_inv_comp.mShares[1](j, 0) = c.mShares[1](0, 0);
            }
            else{
                comp_and_inv_comp.mShares[0](j, 0) = notC.mShares[0](0, 0);
                comp_and_inv_comp.mShares[1](j, 0) = notC.mShares[1](0, 0);
            }
        }

        bool_cipher_and(pIdx, tag_expand, comp_and_inv_comp, tag_next, enc, eval, runtime);
        tag = tag_next;
    }
    
    res = tag;

    return 0;
}

int tagBS(std::vector<aby3::si64> &data, std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    int n = data.size();
    // int bitsize = keyset[0].bitCount();
    int bitsize = key.bitCount();

    if(checkPowerOfTwo(n) == false){
        THROW_RUNTIME_ERROR("Currently, the size of data should be power of 2!");
    }
    int logn = log2(n);
    // get the mid point.
    int p = n/2 - 1;
    sbMatrix c(1, 1);
    boolShare q(1, pIdx);
    sbMatrix _key(1, bitsize);
    _key.mShares[0](0, 0) = keyset[p].mData[0];
    _key.mShares[1](0, 0) = keyset[p].mData[1];
    bool_cipher_lt(pIdx, _key, key, c, enc, eval, runtime);
    bool_cipher_not(pIdx, c, c);

    si64Matrix pivot(1, 1);
    pivot.mShares[0](0, 0) = data[p].mData[0] - data[n-1].mData[0];
    pivot.mShares[1](0, 0) = data[p].mData[1] - data[n-1].mData[1];
    si64Matrix _pivot(1, 1);
    arith_bool_mul(pIdx, pivot, c, _pivot, enc, eval, runtime);
    pivot.mShares[0](0, 0) = _pivot.mShares[0](0, 0) + data[n-1].mData[0];
    pivot.mShares[1](0, 0) = _pivot.mShares[1](0, 0) + data[n-1].mData[1];

    aby3::sbMatrix prev_tag(1, 1);
    prev_tag.mShares[0](0, 0) = q.bshares[0];
    prev_tag.mShares[1](0, 0) = q.bshares[1];

    // log-round binary search.
    for(size_t i=1; i<logn; i++){
        int tag_size = (1 << i);
        aby3::sbMatrix tag(tag_size, 1);
        aby3::sbMatrix comp_inv(1, 1); bool_cipher_not(pIdx, c, comp_inv);

        aby3::sbMatrix prev_tag_expand(tag_size, 1);
        aby3::sbMatrix prev_comp_and_inv_comp(tag_size, 1);
    
        for(size_t j=0; j<tag_size; j++){
            prev_tag_expand.mShares[0](j, 0) = prev_tag.mShares[0](j/2, 0);
            prev_tag_expand.mShares[1](j, 0) = prev_tag.mShares[1](j/2, 0);
            if(j%2 == 0){
                prev_comp_and_inv_comp.mShares[0](j, 0) = c.mShares[0](0, 0);
                prev_comp_and_inv_comp.mShares[1](j, 0) = c.mShares[1](0, 0);
            }
            else{
                prev_comp_and_inv_comp.mShares[0](j, 0) = comp_inv.mShares[0](0, 0);
                prev_comp_and_inv_comp.mShares[1](j, 0) = comp_inv.mShares[1](0, 0);
            }
        }

        // compute the tag.
        bool_cipher_and(pIdx, prev_tag_expand, prev_comp_and_inv_comp, tag, enc, eval, runtime);

        prev_tag = tag;

        // aggregate this level.
        aby3::si64Matrix data_mat(tag_size, 1);
        aby3::sbMatrix key_mat(tag_size, bitsize);
        aby3::sbMatrix keytag(tag_size, bitsize);
        for(size_t j=0; j<tag_size; j++){
            int tag_ind = (((2*(j+1) - 1) * n) / (2*tag_size)) - 1;
            data_mat.mShares[0](j, 0) = data[tag_ind].mData[0];
            data_mat.mShares[1](j, 0) = data[tag_ind].mData[1];
            key_mat.mShares[0](j, 0) = keyset[tag_ind].mData[0];
            key_mat.mShares[1](j, 0) = keyset[tag_ind].mData[1];
            keytag.mShares[0](j, 0) = (tag.mShares[0](j, 0) == 1) ? -1 : 0;
            keytag.mShares[1](j, 0) = (tag.mShares[1](j, 0) == 1) ? -1 : 0;
        }

        // compute the result.
        aby3::si64Matrix res_tmp(tag_size, 1);
        aby3::sbMatrix res_key_tmp(tag_size, bitsize);
        arith_bool_mul(pIdx, data_mat, tag, res_tmp, enc, eval, runtime);
        
        aby3::si64Matrix a(1, 1);
        aby3::sbMatrix akey(1, bitsize);
        a.mShares[0](0, 0) = res_tmp.mShares[0](0, 0);
        a.mShares[1](0, 0) = res_tmp.mShares[1](0, 0);
        for(size_t j=1; j<tag_size; j++){
            a.mShares[0](0, 0) += res_tmp.mShares[0](j, 0);
            a.mShares[1](0, 0) += res_tmp.mShares[1](j, 0);
        }

        std::vector<sbMatrix> enc_key_mat = {key_mat};
        std::vector<sbMatrix> enc_key_tag = {keytag};
        constant_bool_dot(pIdx, enc_key_mat, enc_key_tag, akey, enc, eval, runtime);

        // update the pivot and the comp to the next level.
        bool_cipher_lt(pIdx, akey, key, c, enc, eval, runtime);
        bool_cipher_not(pIdx, c, c);

        aby3::si64Matrix pivot_tmp(1, 1);
        pivot_tmp.mShares[0](0, 0) = a.mShares[0](0, 0) - pivot.mShares[0](0, 0);
        pivot_tmp.mShares[1](0, 0) = a.mShares[1](0, 0) - pivot.mShares[1](0, 0);
        arith_bool_mul(pIdx, pivot_tmp, c, a, enc, eval, runtime);
        pivot.mShares[0](0, 0) += a.mShares[0](0, 0);
        pivot.mShares[1](0, 0) += a.mShares[1](0, 0);
    }

    // get the result.
    res = pivot;

    return 0;
}

int subHBS(std::vector<aby3::si64> &data, std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, int alpha, int threshold){
    int n = data.size();
    // int bitsize = keyset[0].bitCount();
    int bitsize = key.bitCount();
    int beta = n / alpha;
    if(alpha * beta != n) THROW_RUNTIME_ERROR("The size of data should be alpha * beta!");

    if(n < threshold){
        return compBS(data, keyset, key, res, pIdx, enc, eval, runtime);
    }

    // construct the upper-level tree.
    std::vector<sb64> keyset_upper;
    for(size_t i=0; i<alpha; i++){
        keyset_upper.push_back(keyset[i*beta + (beta-1)]);
    }

    // call MBS.
    aby3::sbMatrix res_upper(alpha, 1); // alpha * 1
    mcompBS(keyset_upper, key, res_upper, pIdx, enc, eval, runtime);
    // b2a
    aby3::si64Matrix all_ones(alpha, 1);
    init_ones(pIdx, enc, runtime, all_ones, alpha);
    aby3::si64Matrix res_upper_a(alpha, 1);
    arith_bool_mul(pIdx, all_ones, res_upper, res_upper_a, enc, eval, runtime);

    // constant inner product for data.
    std::vector<si64Matrix> data_lowers(beta);
    std::vector<si64Matrix> tag_lowers(beta);
    si64Matrix data_res(beta, 1);
    for(size_t i=0; i<beta; i++){
        tag_lowers[i] = res_upper_a;
        data_lowers[i].resize(alpha, 1);
        for(size_t j=0; j<alpha; j++){
            data_lowers[i].mShares[0](j, 0) = data[j*beta + i].mData[0];
            data_lowers[i].mShares[1](j, 0) = data[j*beta + i].mData[1];
        }
    }
    constant_sint_dot(pIdx, data_lowers, tag_lowers, data_res, enc, eval, runtime);


    // constant inner product for key.
    std::vector<sbMatrix> key_lowers(beta);
    std::vector<sbMatrix> keytag_lowers(beta);
    sbMatrix key_res(beta, bitsize);
    for(size_t i=0; i<beta; i++){
        key_lowers[i].resize(alpha, bitsize);
        keytag_lowers[i].resize(alpha, bitsize);
        for(size_t j=0; j<alpha; j++){
            key_lowers[i].mShares[0](j, 0) = keyset[j*beta + i].mData[0];
            key_lowers[i].mShares[1](j, 0) = keyset[j*beta + i].mData[1];
            keytag_lowers[i].mShares[0](j, 0) = (res_upper.mShares[0](j, 0) == 1) ? -1 : 0;
            keytag_lowers[i].mShares[1](j, 0) = (res_upper.mShares[1](j, 0) == 1) ? -1 : 0;
        }
    }
    constant_bool_dot(pIdx, key_lowers, keytag_lowers, key_res, enc, eval, runtime);

    // call BS.
    std::vector<si64> data_to_be_search(beta);
    std::vector<sb64> keyset_to_be_search(beta);

    for(size_t i=0; i<beta; i++){
        data_to_be_search[i].mData[0] = data_res.mShares[0](i, 0);
        data_to_be_search[i].mData[1] = data_res.mShares[1](i, 0);
        keyset_to_be_search[i].mData[0] = key_res.mShares[0](i, 0);
        keyset_to_be_search[i].mData[1] = key_res.mShares[1](i, 0);
    }

    // int sub_alpha = sqrtToPowerOfTwo(beta);
    return subHBS(data_to_be_search, keyset_to_be_search, key, res, pIdx, enc, eval, runtime, alpha, threshold);
}

int ptaBS(std::vector<aby3::si64> &data, std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, ABY3MPITask<aby3::sb64, aby3::sb64, aby3::si64, aby3::si64, PtABS>* ptaTask){

    int n = 1, m = data.size();
    ptaTask->set_default_value(GET_ZERO_SHARE);
    ptaTask->set_lookahead(1, 0);
    std::vector<si64> _res(1);

    Timer& timer = Timer::getInstance();

    sb64* _inputx_ptr;
    sb64* _inputy_ptr;
    si64* _selectv_ptr;

    ptaTask->signal_data_size(n, m);
    ptaTask->circuit_construct({n}, {m});
    ptaTask->subTask->set_key_bitsize(32);

    std::vector<si64> _data_vec = data;
    std::vector<sb64> _keyset_vec = keyset;
    std::vector<sb64> _target_key(1);

#ifdef FUNCTIONAL
    if(ptaTask->rank == 0){
        std::reverse(_data_vec.begin(), _data_vec.end());
        std::reverse(_keyset_vec.begin(), _keyset_vec.end());
    }
#endif

    if(ptaTask->rank == 0){
        _target_key[0].mData[0] = key.mShares[0](0, 0);
        _target_key[0].mData[1] = key.mShares[1](0, 0);
    }

    int m_size = ptaTask->subTask->get_partial_m_lens();
    if(ptaTask->rank == 0){
        m_size = m;
        _inputx_ptr = _target_key.data();
        _inputy_ptr = _keyset_vec.data();
        _selectv_ptr = _data_vec.data();
    }
    else{
        _inputx_ptr = new sb64[1];
        _inputy_ptr = new sb64[m_size];
        _selectv_ptr = new si64[m_size];
    }

    ptaTask->data_sharing<sb64>(_inputx_ptr, 1, 0);
    ptaTask->data_sharing<sb64>(_inputy_ptr, m_size, 1);
    ptaTask->data_sharing<si64>(_selectv_ptr, m_size, 1);

    if(ptaTask->rank == 0){
        m_size = ptaTask->subTask->get_partial_m_lens();
        _data_vec.resize(m_size);
        _keyset_vec.resize(m_size);
        _data_vec.shrink_to_fit();
        _keyset_vec.shrink_to_fit();
        _inputy_ptr = _keyset_vec.data();
        _selectv_ptr = _data_vec.data();
    }

    _res[0].mData[0] = 0; _res[0].mData[1] = 0;
    ptaTask->set_selective_value(_selectv_ptr, 0);

    if(ptaTask->rank == 0){
        ptaTask->circuit_evaluate(_inputx_ptr, _keyset_vec.data(), _data_vec.data(), _res.data());
    }
    else{
        ptaTask->circuit_evaluate(_inputx_ptr, _inputy_ptr, _selectv_ptr, _res.data());
    }

    if(ptaTask->rank == 0){
        res.resize(1, 1);
        res.mShares[0](0, 0) = _res[0].mData[0]; res.mShares[1](0, 0) = _res[0].mData[1];
    }

    return 0;
}

int ptaMBS(std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, ABY3MPIPairOnlyTask<sb64, sb64, sb64, sb64, PtAMBS>* ptaTask){
    int n = 1, m = keyset.size();
    ptaTask->set_lookahead(1, 0);
    ptaTask->signal_data_size(n, m);
    ptaTask->circuit_construct({n}, {m});

    ptaTask->subTask->set_key_bitsize(32);
    std::vector<sb64> _res;

    if(ptaTask->rank == 0){
        std::vector<sb64> _keyset_vec = keyset;
#ifdef FUNCTIONAL
        std::reverse(_keyset_vec.begin(), _keyset_vec.end());
#endif
        std::vector<sb64> _target_key(1);
        _target_key[0].mData[0] = key.mShares[0](0, 0);
        _target_key[0].mData[1] = key.mShares[1](0, 0);

        // synchronize the data and keyset with others.
        sb64* _inputx_ptr = _target_key.data();
        sb64* _inputy_ptr = _keyset_vec.data();
        ptaTask->data_sharing<sb64>(_inputx_ptr, 1, 0, 1);
        ptaTask->data_sharing<sb64>(_inputy_ptr, m, 1);
        int m_size = ptaTask->subTask->get_partial_m_lens();
        _keyset_vec.resize(m_size);
        ptaTask->circuit_evaluate(_target_key.data(), _keyset_vec.data(), nullptr, _res.data());
    }
    else{
        int m_size = ptaTask->subTask->get_partial_m_lens();
        sb64* _inputx_ptr = new sb64[1];
        sb64* _inputy_ptr = new sb64[m_size];
        ptaTask->data_sharing<sb64>(_inputx_ptr, 1, 0, 1);
        ptaTask->data_sharing<sb64>(_inputy_ptr, m_size, 1);
        ptaTask->circuit_evaluate(_inputx_ptr, _inputy_ptr, nullptr, _res.data());
        return 0;
    }

    // reorganize the results.
    res.resize(m, 1);
    for(int i=0; i<m; i++){
        res.mShares[0](i, 0) = ptaTask->table[m-i-1].mData[0];
        res.mShares[1](i, 0) = ptaTask->table[m-i-1].mData[1];
    }

    return 0;
}

int ptaMBSAB(std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::sbMatrix &res, aby3::si64Matrix& res_a, ABY3MPIPairOnlyTask<sb64, sb64, std::pair<si64, sb64>, std::pair<si64, sb64>, PtAMBSAB>* ptaTask){
    int n = 1, m = keyset.size();
    ptaTask->set_lookahead(1, 0);
    ptaTask->signal_data_size(n, m);
    ptaTask->circuit_construct({n}, {m});

    ptaTask->subTask->set_key_bitsize(32);
    std::vector<std::pair<si64, sb64>> _res;

    if(ptaTask->rank == 0){
        // std::vector<sb64> _keyset_vec(m);
        std::vector<sb64> _keyset_vec = keyset;
#ifdef FUNCTIONAL
        std::reverse(_keyset_vec.begin(), _keyset_vec.end());
#endif
        std::vector<sb64> _target_key(1);
        _target_key[0].mData[0] = key.mShares[0](0, 0);
        _target_key[0].mData[1] = key.mShares[1](0, 0);

        // synchronize the data and keyset with others.
        sb64* _inputx_ptr = _target_key.data();
        sb64* _inputy_ptr = _keyset_vec.data();
        ptaTask->data_sharing<sb64>(_inputx_ptr, 1, 0, 1);
        ptaTask->data_sharing<sb64>(_inputy_ptr, m, 1);
        int m_size = ptaTask->subTask->get_partial_m_lens();
        _keyset_vec.resize(m_size);
        ptaTask->circuit_evaluate(_target_key.data(), _keyset_vec.data(), nullptr, _res.data());
    }
    else{
        int m_size = ptaTask->subTask->get_partial_m_lens();
        sb64* _inputx_ptr = new sb64[1];
        sb64* _inputy_ptr = new sb64[m_size];
        ptaTask->data_sharing<sb64>(_inputx_ptr, 1, 0, 1);
        ptaTask->data_sharing<sb64>(_inputy_ptr, m_size, 1);
        ptaTask->circuit_evaluate(_inputx_ptr, _inputy_ptr, nullptr, _res.data());
        return 0;
    }

    // reorganize the results.
    res.resize(m, 1);
    res_a.resize(m, 1);
    for(int i=0; i<m; i++){
        res.mShares[0](i, 0) = ptaTask->table[m-i-1].second.mData[0];
        res.mShares[1](i, 0) = ptaTask->table[m-i-1].second.mData[1];
        res_a.mShares[0](i, 0) = ptaTask->table[m-i-1].first.mData[0];
        res_a.mShares[1](i, 0) = ptaTask->table[m-i-1].first.mData[1];
    }

    return 0;
}

int subHBS_with_PtA(std::vector<aby3::si64> &data, std::vector<aby3::sb64> &keyset, aby3::sbMatrix &key, aby3::si64Matrix &res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, int alpha, int threshold, ABY3MPITask<aby3::sb64, aby3::sb64, aby3::si64, aby3::si64, PtABS>* ptaBSTask, ABY3MPIPairOnlyTask<aby3::sb64, aby3::sb64, std::pair<aby3::si64, aby3::sb64>, std::pair<aby3::si64, aby3::sb64>, PtAMBSAB>* ptaMBSTask){  
    int n = data.size();
    int rank = ptaBSTask->rank;

    // only for rank0, the size is the correct data size, therefore broadcast.
    ptaBSTask->signal_data_size(n);
    if(n < threshold){
        // for too small input size, only rank0 compute is fine.
        if(n < ptaBSTask->total_tasks * PER_PARTY_MAX){
            if(ptaBSTask->rank == 0){
                return compBS(data, keyset, key, res, pIdx, enc, eval, runtime);
            }
            else{
                return 0;
            }
        }
        // otherwise call the ptaBS.
        return ptaBS(data, keyset, key, res, ptaBSTask);
    }

    int bitsize = key.bitCount();
    int beta = n / alpha;
    if(alpha * beta != n) THROW_RUNTIME_ERROR("The size of data should be alpha * beta, while now alpha = " + std::to_string(alpha) + " - beta = " + std::to_string(beta) + " - n = " + std::to_string(n) + " However alpha * beta = " + std::to_string(alpha * beta));

    std::vector<sb64> keyset_upper;
    aby3::sbMatrix res_upper(alpha, 1);
    aby3::si64Matrix res_upper_a(alpha, 1);
    if(rank == 0){
        // construct the upper-level tree.
        for(size_t i=0; i<alpha; i++){
            keyset_upper.push_back(keyset[i*beta + (beta-1)]);
        }
    }
    ptaMBSAB(keyset_upper, key, res_upper, res_upper_a, ptaMBSTask);

    // the p0 task perform the following computation till recursive.
    // parallelize the computation by beta. 
    // bcast res_upper to all parties.
    MPI_Bcast(res_upper.mShares[0].data(), (int)res_upper.i64Size(), MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(res_upper.mShares[1].data(), (int)res_upper.i64Size(), MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(res_upper_a.mShares[0].data(), alpha, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(res_upper_a.mShares[1].data(), alpha, MPI_INT64_T, 0, MPI_COMM_WORLD);

    // scatterv the data to all the parties for parallel computation.
    int local_beta = beta / ptaBSTask->total_tasks;
    if(rank == ptaBSTask->total_tasks - 1){
        local_beta = beta - local_beta * (ptaBSTask->total_tasks - 1);
    }
    std::vector<int> sendcounts(ptaBSTask->total_tasks);
    std::vector<int> displs(ptaBSTask->total_tasks);
    std::vector<int> sendcounts_alpha(ptaBSTask->total_tasks);
    std::vector<int> displs_alpha(ptaBSTask->total_tasks);
    if(rank == 0){
        for(int i=0; i<ptaBSTask->total_tasks - 1; i++){
            sendcounts[i] = local_beta;
            displs[i] = i * local_beta;
            sendcounts_alpha[i] = local_beta * alpha;
            displs_alpha[i] = i * local_beta * alpha;
        }
        sendcounts[ptaBSTask->total_tasks - 1] = (beta - (ptaBSTask->total_tasks - 1) * local_beta);
        sendcounts_alpha[ptaBSTask->total_tasks - 1] = (beta - (ptaBSTask->total_tasks - 1) * local_beta) * alpha;
        displs[ptaBSTask->total_tasks - 1] = (ptaBSTask->total_tasks - 1) * local_beta;
        displs_alpha[ptaBSTask->total_tasks - 1] = (ptaBSTask->total_tasks - 1) * local_beta * alpha;
    }

    si64Matrix data_to_be_search(local_beta*alpha, 1);
    sbMatrix keyset_to_be_search(local_beta*alpha, bitsize);
    if(rank == 0){
        si64Matrix all_data(beta*alpha, 1);
        for(size_t i=0; i<beta; i++){
            for(size_t j=0; j<alpha; j++){
                all_data.mShares[0](i*alpha + j, 0) = data[j*beta + i].mData[0];
                all_data.mShares[1](i*alpha + j, 0) = data[j*beta + i].mData[1];
            }
        }

        // scatter the data.
        MPI_Scatterv(all_data.mShares[0].data(), sendcounts_alpha.data(), displs_alpha.data(), MPI_INT64_T, data_to_be_search.mShares[0].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(all_data.mShares[1].data(), sendcounts_alpha.data(), displs_alpha.data(), MPI_INT64_T, data_to_be_search.mShares[1].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);

        data.resize(0);
        all_data.resize(0, 0);
        data.shrink_to_fit();

        sbMatrix all_keyset(beta*alpha, bitsize);
        for(size_t i=0; i<beta; i++){
            for(size_t j=0; j<alpha; j++){
                all_keyset.mShares[0](i*alpha + j, 0) = keyset[j*beta + i].mData[0];
                all_keyset.mShares[1](i*alpha + j, 0) = keyset[j*beta + i].mData[1];
            }
        }
        // scatter the data.
        MPI_Scatterv(all_keyset.mShares[0].data(), sendcounts_alpha.data(), displs_alpha.data(), MPI_INT64_T, keyset_to_be_search.mShares[0].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(all_keyset.mShares[1].data(), sendcounts_alpha.data(), displs_alpha.data(), MPI_INT64_T, keyset_to_be_search.mShares[1].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);

        keyset.resize(0);
        all_keyset.resize(0, 0);
        keyset.shrink_to_fit();
    }
    else{
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT64_T, data_to_be_search.mShares[0].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT64_T, data_to_be_search.mShares[1].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT64_T, keyset_to_be_search.mShares[0].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT64_T, keyset_to_be_search.mShares[1].data(), alpha*local_beta, MPI_INT64_T, 0, MPI_COMM_WORLD);
    }

    // each task independetly compute the result.
    std::vector<si64Matrix> data_lowers(local_beta);
    std::vector<si64Matrix> tag_lowers(local_beta);
    si64Matrix data_res(local_beta, 1);
    for(size_t i=0; i<local_beta; i++){
        tag_lowers[i] = res_upper_a;
        data_lowers[i].resize(alpha, 1);
        for(size_t j=0; j<alpha; j++){
            data_lowers[i].mShares[0](j, 0) = data_to_be_search.mShares[0](i*alpha + j, 0);
            data_lowers[i].mShares[1](j, 0) = data_to_be_search.mShares[1](i*alpha + j, 0);
        }
    }
    constant_sint_dot(pIdx, data_lowers, tag_lowers, data_res, enc, eval, runtime);

    // constant inner product for key.
    std::vector<sbMatrix> key_lowers(local_beta);
    std::vector<sbMatrix> keytag_lowers(local_beta);
    sbMatrix key_res(local_beta, bitsize);
    for(size_t i=0; i<local_beta; i++){
        key_lowers[i].resize(alpha, bitsize);
        keytag_lowers[i].resize(alpha, bitsize);
        for(size_t j=0; j<alpha; j++){
            key_lowers[i].mShares[0](j, 0) = keyset_to_be_search.mShares[0](i*alpha + j, 0);
            key_lowers[i].mShares[1](j, 0) = keyset_to_be_search.mShares[1](i*alpha + j, 0);
            keytag_lowers[i].mShares[0](j, 0) = (res_upper.mShares[0](j, 0) == 1) ? -1 : 0;
            keytag_lowers[i].mShares[1](j, 0) = (res_upper.mShares[1](j, 0) == 1) ? -1 : 0;
        }
    }
    constant_bool_dot(pIdx, key_lowers, keytag_lowers, key_res, enc, eval, runtime);

    // combine together to the task0.
    sbMatrix key_res_all;
    si64Matrix data_res_all;
    std::vector<si64> data_to_be_search_all;
    std::vector<sb64> keyset_to_be_search_all;
    if(rank == 0){
        key_res_all.resize(beta, bitsize);
        data_res_all.resize(beta, 1);

        std::copy(key_res.mShares[0].data(), key_res.mShares[0].data() + local_beta, key_res_all.mShares[0].data());
        std::copy(key_res.mShares[1].data(), key_res.mShares[1].data() + local_beta, key_res_all.mShares[1].data());
        std::copy(data_res.mShares[0].data(), data_res.mShares[0].data() + local_beta, data_res_all.mShares[0].data());
        std::copy(data_res.mShares[1].data(), data_res.mShares[1].data() + local_beta, data_res_all.mShares[1].data());

        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, key_res_all.mShares[0].data(), sendcounts.data(), displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, key_res_all.mShares[1].data(), sendcounts.data(), displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, data_res_all.mShares[0].data(), sendcounts.data(), displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_INT64_T, data_res_all.mShares[1].data(), sendcounts.data(), displs.data(), MPI_INT64_T, 0, MPI_COMM_WORLD);

        data_to_be_search_all.resize(beta);
        keyset_to_be_search_all.resize(beta);
        for(size_t i=0; i<beta; i++){
            data_to_be_search_all[i].mData[0] = data_res_all.mShares[0](i, 0);
            data_to_be_search_all[i].mData[1] = data_res_all.mShares[1](i, 0);
            keyset_to_be_search_all[i].mData[0] = key_res_all.mShares[0](i, 0);
            keyset_to_be_search_all[i].mData[1] = key_res_all.mShares[1](i, 0);
        }
    }
    else{
        MPI_Gatherv(key_res.mShares[0].data(), local_beta, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(key_res.mShares[1].data(), local_beta, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(data_res.mShares[0].data(), local_beta, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Gatherv(data_res.mShares[1].data(), local_beta, MPI_INT64_T, nullptr, nullptr, nullptr, MPI_INT64_T, 0, MPI_COMM_WORLD);
    }

    return subHBS_with_PtA(data_to_be_search_all, keyset_to_be_search_all, key, res, pIdx, enc, eval, runtime, alpha, threshold, ptaBSTask, ptaMBSTask);
}