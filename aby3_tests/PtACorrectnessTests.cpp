#include "Test.h"
#include <mpi.h>

#include "../aby3-RTR/PtATests.h"
#include "../aby3-RTR/GASTest.h"
#include "../aby3-Basic/timer.h"
#include "../aby3-Basic/Basics.h"

using namespace oc;
using namespace aby3;


int correctness_cipher_index_pta(oc::CLP& cmd){
    size_t n = 1, m=1<<10, optB=128;
    int task_num;
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);

    SETUP_PROCESS
    // construct the task.
    auto ptaTask = new ABY3MPITask<si64, int, si64, si64, CipherIndex>(task_num, optB, role, enc, runtime, eval);
    ptaTask->circuit_construct({n}, {m});
    ptaTask->set_default_value(GET_ZERO_SHARE);

    // data loading.
    std::vector<si64> inputX(n); std::vector<int> inputY; 
    std::vector<si64> inputV; std::vector<si64> res(n);

    // generate the test data.
    size_t partial_len = ptaTask->subTask->get_partial_m_lens();
    i64Matrix pinput_x(n, 1); i64Matrix pinput_v(partial_len, 1);
    pinput_x(0, 0) = 1;
    inputY.resize(partial_len); inputV.resize(partial_len);
    for(size_t i=0; i<partial_len; i++){
        inputY[i] = i + ptaTask->m_start;
        pinput_v(i, 0) = i + ptaTask->m_start;
    }

    aby3::si64Matrix sinputX(n, 1);
    aby3::si64Matrix sinputV(partial_len, 1);
    if(role == 0){
        enc.localIntMatrix(runtime, pinput_x, sinputX).get();
        enc.localIntMatrix(runtime, pinput_v, sinputV).get();
    }
    else{
        enc.remoteIntMatrix(runtime, sinputX).get();
        enc.remoteIntMatrix(runtime, sinputV).get();
    }
    inputX[0].mData[0] = sinputX.mShares[0](0, 0); 
    inputX[0].mData[1] = sinputX.mShares[1](0, 0);
    
    res[0] = GET_ZERO_SHARE;
    for(size_t i=0; i<partial_len; i++){
        // inputV[i] = sinputV(i, 0);
        inputV[i].mData[0] = sinputV.mShares[0](i, 0);
        inputV[i].mData[1] = sinputV.mShares[1](i, 0);
    }

    // evaluate the task.
    ptaTask->set_selective_value(inputV.data(), 0);
    ptaTask->circuit_evaluate(inputX.data(), inputY.data(), inputV.data(), res.data());

    // check the result.   
    aby3::i64Matrix test_res(1, 1);
    aby3::si64Matrix runtime_res(1, 1);
    // runtime_res(0, 0) = res[0];
    runtime_res.mShares[0](0, 0) = res[0].mData[0];
    runtime_res.mShares[1](0, 0) = res[0].mData[1];
    enc.revealAll(runtime, runtime_res, test_res).get();

    if(rank == 0 && role == 0){
        check_result("pta cipher_index", test_res, pinput_x);
    }
    return 0;
}

int correctness_sort_pta(oc::CLP& cmd){
    size_t n = 50, m=50, optB=128;
    int task_num;
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);
    MPI_Barrier(MPI_COMM_WORLD);
    SETUP_PROCESS
    // construct the task.
    auto ptaTask_step1 = new ABY3MPITask<si64, si64, si64, si64, Rank>(task_num, optB, role, enc, runtime, eval);
    ptaTask_step1->set_default_value(GET_ZERO_SHARE);
    ptaTask_step1->circuit_construct({n}, {m});
    auto ptaTask_step2 = new ABY3MPITask<int, si64, si64, si64, CipherIndex>(task_num, optB, role, enc, runtime, eval);
    ptaTask_step2->set_default_value(GET_ZERO_SHARE);
    ptaTask_step2->circuit_construct({n}, {m});

    // data_construction.
    i64Matrix plain_data(n, 1);
    std::vector<int> indices(n);
    std::vector<si64> middle_res(n);
    si64Matrix init_zero(n, 1);
    init_zeros(role, enc, runtime, init_zero, n);
    std::vector<si64> res(n);

    for(size_t i=0; i<n; i++){
        indices[i] = n-i-1;
        plain_data(i, 0) = i;
        middle_res[i].mData[0] = init_zero.mShares[0](i, 0);
        middle_res[i].mData[1] = init_zero.mShares[1](i, 0);
        res[i].mData[0] = init_zero.mShares[0](i, 0);
        res[i].mData[1] = init_zero.mShares[1](i, 0);
    }

    si64Matrix sdata(n, 1);
    if(role == 0){
        enc.localIntMatrix(runtime, plain_data, sdata).get();
    }
    else{
        enc.remoteIntMatrix(runtime, sdata).get();
    }

    std::vector<si64> inputX(n);
    for(size_t i=0; i<n; i++){
        inputX[i].mData[0] = sdata.mShares[0](i, 0);
        inputX[i].mData[1] = sdata.mShares[1](i, 0);
    }

    size_t partial_len = ptaTask_step1->subTask->get_partial_m_lens();
    size_t m_start = ptaTask_step1->m_start;

    std::vector<si64> partialX(partial_len);
    for(size_t i=0; i<partial_len; i++){
        partialX[i].mData[0] = sdata.mShares[0](m_start + i, 0);
        partialX[i].mData[1] = sdata.mShares[1](m_start + i, 0);
    }

    // evaluate the task (1 step)
    ptaTask_step1->circuit_evaluate(inputX.data(), partialX.data(), nullptr, middle_res.data());

    ptaTask_step2->set_selective_value(partialX.data(), 0);

    // resharing the result and prepare for the next step.
    size_t sharing_length;
    aby3::si64* partial_data;
    if(rank == 0){
        sharing_length = ptaTask_step1->m;
        partial_data = middle_res.data();
    }
    else{
        sharing_length = (ptaTask_step1->m_end - ptaTask_step1->m_start) + 1;
        partial_data = new aby3::si64[sharing_length];
    }
    ptaTask_step2->data_sharing<aby3::si64>(partial_data, sharing_length, 1);

    // debug
    std::vector<aby3::si64> _partial_data(sharing_length);
    for(size_t i=0; i<sharing_length; i++){
        _partial_data[i].mData[0] = partial_data[i].mData[0];
        _partial_data[i].mData[1] = partial_data[i].mData[1];
    }

    ptaTask_step2->circuit_evaluate(indices.data(), partial_data, partialX.data(), res.data());

    // check the result.
    // debug_output_vector_mpi(role, rank, res, runtime, enc, "sort res");
    aby3::i64Matrix test_res(n, 1);
    for(size_t i=0; i<n; i++){
        test_res(i, 0) = i;
    }
    aby3::si64Matrix runtime_res(n, 1);
    aby3::i64Matrix res_mat(n, 1);
    vec2mat(res, runtime_res);
    enc.revealAll(runtime, runtime_res, res_mat).get();

    if(rank == 0 && role == 0){
        check_result("pta sort", test_res, res_mat);
    }

    return 0;
}

int correctness_sum_pta(oc::CLP& cmd){
    size_t n = 1, m=1<<10, optB=128;
    int task_num;
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);
    MPI_Barrier(MPI_COMM_WORLD);
    SETUP_PROCESS

    // construct the task.
    auto ptaTask = new ABY3MPITask<si64, si64, si64, si64, Sum>(task_num, optB, role, enc, runtime, eval);
    ptaTask->set_default_value(GET_ZERO_SHARE);
    ptaTask->circuit_construct({n}, {m});

    // data_construction.
    i64Matrix plain_data(m, 1);
    i64Matrix sum(1, 1);
    sum(0, 0) = 0;
    for(size_t i=0; i<m; i++){
        plain_data(i, 0) = i;
        sum(0, 0) += i;
    }
    si64Matrix sdata(m, 1);
    if(role == 0){
        enc.localIntMatrix(runtime, plain_data, sdata).get();
    }
    else{
        enc.remoteIntMatrix(runtime, sdata).get();
    }

    std::vector<si64> inputX(n);
    std::vector<si64> res(n);
    inputX[0].mData[0] = 0;
    inputX[0].mData[1] = 0;
    res[0].mData[0] = 0;
    res[0].mData[1] = 0;

    size_t partial_len = ptaTask->subTask->get_partial_m_lens();
    size_t m_start = ptaTask->m_start;

    std::vector<si64> partialX(partial_len);
    for(size_t i=0; i<partial_len; i++){
        partialX[i].mData[0] = sdata.mShares[0](m_start + i, 0);
        partialX[i].mData[1] = sdata.mShares[1](m_start + i, 0);
    }

    // evaluate the task.
    ptaTask->circuit_evaluate(inputX.data(), partialX.data(), nullptr, res.data());

    // check the result.
    aby3::i64Matrix test_res(1, 1);
    aby3::si64Matrix runtime_res(1, 1);
    runtime_res.mShares[0](0, 0) = res[0].mData[0];
    runtime_res.mShares[1](0, 0) = res[0].mData[1];
    enc.revealAll(runtime, runtime_res, test_res).get();

    if(rank == 0 && role == 0){
        check_result("pta sum", test_res, sum);
    }

    return 0;
}


int correctness_max_pta(oc::CLP& cmd){
    size_t n = 1, m=1<<10, optB=128;
    int task_num;
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);
    MPI_Barrier(MPI_COMM_WORLD);
    SETUP_PROCESS

    // construct the task.
    auto ptaTask = new ABY3MPITask<si64, si64, si64, si64, Max>(task_num, optB, role, enc, runtime, eval);
    ptaTask->set_default_value(GET_ZERO_SHARE);
    ptaTask->circuit_construct({n}, {m});

    // data_construction.
    i64Matrix plain_data(m, 1);
    i64Matrix max(1, 1);
    max(0, 0) = m-1;
    for(size_t i=0; i<m; i++){
        plain_data(i, 0) = i;
    }
    si64Matrix sdata(m, 1);
    if(role == 0){
        enc.localIntMatrix(runtime, plain_data, sdata).get();
    }
    else{
        enc.remoteIntMatrix(runtime, sdata).get();
    }

    std::vector<si64> inputX(n);
    std::vector<si64> res(n);
    inputX[0].mData[0] = 0;
    inputX[0].mData[1] = 0;
    res[0].mData[0] = 0;
    res[0].mData[1] = 0;

    size_t partial_len = ptaTask->subTask->get_partial_m_lens();
    size_t m_start = ptaTask->m_start;

    std::vector<si64> partialX(partial_len);
    for(size_t i=0; i<partial_len; i++){
        partialX[i].mData[0] = sdata.mShares[0](m_start + i, 0);
        partialX[i].mData[1] = sdata.mShares[1](m_start + i, 0);
    }

    // evaluate the task.
    ptaTask->circuit_evaluate(inputX.data(), partialX.data(), nullptr, res.data());

    // check the result.
    aby3::i64Matrix test_res(1, 1);
    aby3::si64Matrix runtime_res(1, 1);
    runtime_res.mShares[0](0, 0) = res[0].mData[0];
    runtime_res.mShares[1](0, 0) = res[0].mData[1];
    enc.revealAll(runtime, runtime_res, test_res).get();

    if(rank == 0 && role == 0){
        check_result("pta max", test_res, max);
    }

    return 0;
}

int correctness_metric_pta(oc::CLP& cmd){

    size_t n = 1, m=1<<10, optB=128, k=1;
    int task_num;
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);
    MPI_Barrier(MPI_COMM_WORLD);
    SETUP_PROCESS

    // construct the task.
    auto ptaTask = new ABY3MPITask<std::vector<si64>, std::vector<si64>, si64, si64, BioMetric>(task_num, optB, role, enc, runtime, eval);
    ptaTask->set_default_value(get_share(role, m+1));
    ptaTask->circuit_construct({n}, {m});

    // data_construction.
    i64Matrix plain_data(m, 1);
    i64Matrix min_dis(1, 1);
    min_dis(0, 0) = 1;
    i64Matrix target(1, 1);
    target(0, 0) = 0;
    for(size_t i=0; i<m; i++){
        plain_data(i, 0) = i + 1;
    }
    si64Matrix sdata(m, 1);
    si64Matrix starget(1, 1);
    if(role == 0){
        enc.localIntMatrix(runtime, plain_data, sdata).get();
        enc.localIntMatrix(runtime, target, starget).get();
    }
    else{
        enc.remoteIntMatrix(runtime, sdata).get();
        enc.remoteIntMatrix(runtime, starget).get();
    }

    std::vector<std::vector<si64>> inputX(n, std::vector<si64>(1));
    std::vector<si64> res(n);
    inputX[0][0].mData[0] = starget.mShares[0](0, 0);
    inputX[0][0].mData[1] = starget.mShares[1](0, 0);
    res[0].mData[0] = 0;
    res[0].mData[1] = 0;

    size_t partial_len = ptaTask->subTask->get_partial_m_lens();
    size_t m_start = ptaTask->m_start;

    std::vector<std::vector<si64>> partialX(partial_len, std::vector<si64>(1));
    for(size_t i=0; i<partial_len; i++){
        partialX[i][0].mData[0] = sdata.mShares[0](m_start + i, 0);
        partialX[i][0].mData[1] = sdata.mShares[1](m_start + i, 0);
    }

    // evaluate the task.
    ptaTask->circuit_evaluate(inputX.data(), partialX.data(), nullptr, res.data());

    // check the result.
    aby3::i64Matrix test_res(1, 1);
    aby3::si64Matrix runtime_res(1, 1);
    runtime_res.mShares[0](0, 0) = res[0].mData[0];
    runtime_res.mShares[1](0, 0) = res[0].mData[1];
    enc.revealAll(runtime, runtime_res, test_res).get();

    if(rank == 0 && role == 0){
        check_result("pta metric", test_res, min_dis);
    }

    return 0;
}

int correctness_search_pta(oc::CLP& cmd){
    int n=1, m=1<<10, optB=128;

    SETUP_PROCESS

    // construct the task.
    auto ptaTask = new ABY3MPITask<sb64, sb64, si64, si64, PtABS>(size, optB, role, enc, runtime, eval);
    ptaTask->set_default_value(GET_ZERO_SHARE);

    // prepare the data pointers.
    std::vector<si64> res(n);
    i64Matrix target(1, 1); target(0, 0) = m-2;

    // p0 prepares the data.
    if(rank == 0){
        // data loading.
        i64Matrix plain_data(m, 1);
        for(size_t i=0; i<m; i++){
            plain_data(i, 0) = i;
        }
        si64Matrix sdata(m, 1);
        sbMatrix skeyset(m, 1);
        sbMatrix starget(1, 1);

        if(role == 0){
            enc.localIntMatrix(runtime, plain_data, sdata).get();
            enc.localBinMatrix(runtime, plain_data, skeyset).get();
            enc.localBinMatrix(runtime, target, starget).get();
        }
        else{
            enc.remoteIntMatrix(runtime, sdata).get();
            enc.remoteBinMatrix(runtime, skeyset).get();
            enc.remoteBinMatrix(runtime, starget).get();
        }

        std::vector<sb64> inputX(n); std::vector<sb64> inputY(m); std::vector<si64> selectV(m);
        // std::vector<si64> res(n);
        inputX[0].mData[0] = starget.mShares[0](0, 0);
        inputX[0].mData[1] = starget.mShares[1](0, 0);
        res[0].mData[0] = 0; res[0].mData[1] = 0;

        for(size_t i=0; i<m; i++){
            inputY[i].mData[0] = skeyset.mShares[0](i, 0);
            inputY[i].mData[1] = skeyset.mShares[1](i, 0);
            selectV[i].mData[0] = sdata.mShares[0](i, 0);
            selectV[i].mData[1] = sdata.mShares[1](i, 0);
        }

        // synchronize with others.
        ptaTask->signal_data_size(n, m);
        ptaTask->set_lookahead(1, 0);
        ptaTask->circuit_construct({n}, {m});
        ptaTask->subTask->set_key_bitsize(32);

        sb64* _inputx_ptr = inputX.data();
        sb64* _inputy_ptr = inputY.data();
        si64* _selectv_ptr = selectV.data();
        ptaTask->data_sharing<sb64>(_inputx_ptr, 1, 0);
        ptaTask->data_sharing<sb64>(_inputy_ptr, m, 1);
        ptaTask->data_sharing<si64>(_selectv_ptr, m, 1);
        
        ptaTask->set_selective_value(selectV.data(), 0);
        ptaTask->circuit_evaluate(inputX.data(), inputY.data(), selectV.data(), res.data());
    }
    else{
        // all the other tasks get the data:
        ptaTask->signal_data_size(n, m);
        ptaTask->set_lookahead(1, 0);
        ptaTask->circuit_construct({n}, {m});
        ptaTask->subTask->set_key_bitsize(32);
        int m_size = ptaTask->subTask->get_partial_m_lens();
        sb64* inputX_ptr = new sb64[n];
        sb64* inputY_ptr = new sb64[m_size];
        si64* selectV_ptr = new si64[m_size];
        // std::vector<si64> res(n);
        res[0].mData[0] = 0; res[0].mData[1] = 0;
        ptaTask->data_sharing<sb64>(inputX_ptr, 1, 0);
        ptaTask->data_sharing<sb64>(inputY_ptr, m_size, 1);
        ptaTask->data_sharing<si64>(selectV_ptr, m_size, 1);

        ptaTask->set_selective_value(selectV_ptr, 0);
        ptaTask->circuit_evaluate(inputX_ptr, inputY_ptr, selectV_ptr, res.data());
    }

    // p0 get the result.
    aby3::i64Matrix test_res(1, 1);
    aby3::si64Matrix runtime_res(1, 1);
    runtime_res.mShares[0](0, 0) = res[0].mData[0];
    runtime_res.mShares[1](0, 0) = res[0].mData[1];
    enc.revealAll(runtime, runtime_res, test_res).get();

    if(rank == 0 && role == 0){
        check_result("pta search", test_res, target);
    }

    return 0;
}