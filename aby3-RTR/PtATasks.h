#pragma once
#include "GeneralPTA.h"
#include "DataGenerator.h"
#include "../aby3-Basic/Basics.h"
// #define PROFILE


const aby3::si64 GET_ZERO_SHARE = [] {
    aby3::si64 dval;
    dval.mData[0] = 0;
    dval.mData[1] = 0;
    return dval;
}();

aby3::si64 get_share(int pIdx, int target_val);

template<typename T>
struct indexData{
  T value;
  aby3::si64 index;
};

// functional class override.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class CipherIndex : public SubTask<NUMX, NUMY, NUMT, NUMR> { 

public:
    // aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    double eq_time = 0;
    double mul_time = 0;
    double other_time = 0;

    std::string logging_file = "";

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    CipherIndex(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = true;
    }

    void set_logging_file(std::string log_file){
        this->logging_file = log_file;
    }

    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        
        auto stamp1 = std::chrono::high_resolution_clock::now();

        aby3::u64 block_length = binfo->block_len;
        // pairwise comparison
        aby3::sbMatrix partTable(block_length, 1); 
        vector_cipher_eq(this->pIdx, expandX, expandY, partTable, *(this->eval), *(this->runtime));

        auto stamp2 = std::chrono::high_resolution_clock::now();

        // pairwise abmul.
        aby3::si64Matrix expandV(block_length, 1), p_v_after_mul(block_length, 1);
        for(size_t i=0; i<block_length; i++) {
            expandV.mShares[0](i, 0) = this->selectV[binfo->t_start + i].mData[0];
            expandV.mShares[1](i, 0) = this->selectV[binfo->t_start + i].mData[1];
        }

        auto stamp3 = std::chrono::high_resolution_clock::now();

        cipher_mul(this->pIdx, expandV, partTable, p_v_after_mul, *this->eval, *this->enc, *this->runtime);

        auto stamp4 = std::chrono::high_resolution_clock::now();

        //trans to local_table.
        for (size_t i = 0; i < block_length; i++) {
            local_table[i].mData[0] = p_v_after_mul.mShares[0](i, 0);
            local_table[i].mData[1] = p_v_after_mul.mShares[1](i, 0);
        }

        auto stamp5 = std::chrono::high_resolution_clock::now();

#ifdef PROFILE
        eq_time = std::chrono::duration_cast<std::chrono::milliseconds>(stamp2 - stamp1).count();
        mul_time = std::chrono::duration_cast<std::chrono::milliseconds>(stamp4 - stamp3).count();
        other_time = std::chrono::duration_cast<std::chrono::milliseconds>(stamp5 - stamp4).count();

        if(this->task_id == 0){
            std::ofstream ofs(this->logging_file, ios::app);
            ofs << "eq_time: " << eq_time << " ms" << endl;
            ofs << "eq processing time: " << comp_process_time.count() / 1000 << " ms" << endl;
            ofs << "eq comm time: " << comp_comm_time.count() / 1000 << " ms" << endl;
            ofs << "mul_time: " << mul_time << " ms" << endl;
            ofs << "other_time: " << other_time << " ms" << endl;
            ofs << "----------" << endl;
            ofs.close();
        }

#endif
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        for (size_t i = 0; i < resLeft.size(); i++) {
            local_res[i] = resLeft[i] + resRight[i];
        }
    }

    std::tuple<std::vector<aby3::si64>, std::vector<int>, std::vector<aby3::si64>> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::tuple<std::vector<aby3::si64>, std::vector<int>, std::vector<aby3::si64>> data_loading(){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecX = generate_vector_si64(this->n, this->pIdx, *this->enc, *this->runtime);
        std::vector<int> vecY = generate_vector_int(partial_len);
        std::vector<aby3::si64> vecV = generate_vector_si64(partial_len,  this->pIdx, *this->enc, *this->runtime);
        return std::make_tuple(vecX, vecY, vecV);
    }
};

// Max function class
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class Max : public SubTask<NUMX, NUMY, NUMT, NUMR> {

public:
    // aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    Max(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
    }

    // functional code.
    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        local_table = expandY;
        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        aby3::sbMatrix comp_res;
        vector_cipher_gt(this->pIdx, resLeft, resRight, comp_res, *(this->eval), *(this->enc), *(this->runtime));
        aby3::si64Matrix mul_mat(resLeft.size(), 1), mul_res(resLeft.size(), 1);
        for(size_t i=0; i<resLeft.size(); i++) {
            mul_mat.mShares[0](i, 0) = resLeft[i].mData[0] - resRight[i].mData[0];
            mul_mat.mShares[1](i, 0) = resLeft[i].mData[1] - resRight[i].mData[1];
        }
        cipher_mul_seq(this->pIdx, mul_mat, comp_res, mul_res, *(this->eval), *(this->enc), *(this->runtime));
        for(size_t i=0; i<resLeft.size(); i++) {
            local_res[i].mData[0] = resRight[i].mData[0] + mul_res.mShares[0](i, 0);
            local_res[i].mData[1] = resRight[i].mData[1] + mul_res.mShares[1](i, 0);
        }
        return;
    }

    std::vector<aby3::si64> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::vector<aby3::si64> data_loading(){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecY = generate_vector_si64(partial_len, this->pIdx, *this->enc, *this->runtime);
        return vecY;
    }

};

// Rank function class.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class Rank : public SubTask<NUMX, NUMY, NUMT, NUMR> {

public:
    //aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    Rank(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
    }

    // functional code.
    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        aby3::u64 block_length = binfo->block_len;
        aby3::sbMatrix partTable(block_length, 1);
        vector_cipher_gt(this->pIdx, expandY, expandX, local_table, *(this->eval), *(this->enc), *(this->runtime));
        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        for(size_t i=0; i<resLeft.size(); i++) local_res[i] = resLeft[i] + resRight[i];
        return;
    }

    std::vector<aby3::si64> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::pair<std::vector<aby3::si64>, std::vector<aby3::si64>> data_loading(){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecX = generate_vector_si64(this->n, this->pIdx, *this->enc, *this->runtime);
        std::vector<aby3::si64> vecY = generate_vector_si64(partial_len, this->pIdx, *this->enc, *this->runtime);
        return std::make_pair(vecX, vecY);
    }

};

// Sum function class.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class Sum : public SubTask<NUMX, NUMY, NUMT, NUMR> {

public:
    //aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    Sum(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
    }

    // functional code.
    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        local_table = expandY;
        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        for(size_t i=0; i<resLeft.size(); i++) local_res[i] = resLeft[i] + resRight[i];
        return;
    }

    std::vector<aby3::si64> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::vector<aby3::si64> data_loading(){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecY = generate_vector_si64(partial_len, this->pIdx, *this->enc, *this->runtime);
        return vecY;
    }
};

// BioMetric function class.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class BioMetric : public SubTask<NUMX, NUMY, NUMT, NUMR> {
public:
    //aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    BioMetric(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
    }

    // functional code.
    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        aby3::u64 block_length = binfo->block_len;
        size_t k = expandX[0].size();

        // flat the two-dimensional inputs.
        size_t exp_len = expandX.size()*k;
        std::vector<typename NUMX::value_type> flatX, flatY;
        for (const auto& innerVec : expandX){
            for (const auto& element : innerVec) flatX.push_back(element);
        }
        for (const auto& innerVec : expandY){
            for (const auto& element : innerVec) flatY.push_back(element);
        }

        // vector mul
        if(std::is_same<typename NUMX::value_type, aby3::si64>::value){
            vector_mean_square(this->pIdx, flatX, flatY, flatX, *(this->eval), *(this->enc), *(this->runtime));
        }

        // delete the space of flatY.
        flatY.clear();
        
        // reduce to local table
        for(int i=0; i<expandX.size(); i++){
            local_table[i] = GET_ZERO_SHARE;
            for(int j=0; j<k; j++){
                local_table[i].mData[0] += flatX[i*k+j].mData[0];
                local_table[i].mData[1] += flatX[i*k+j].mData[1];
            }
        }

        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        aby3::sbMatrix comp_res;

        vector_cipher_gt(this->pIdx, resLeft, resRight, comp_res, *(this->eval), *(this->enc), *(this->runtime));

        // multiply for value extraction.
        aby3::si64Matrix mat_mul(resLeft.size(), 1), mat_res(resLeft.size(), 1);
        for(size_t i=0; i<resLeft.size(); i++) {
            // mat_mul(i, 0, resRight[i] - resLeft[i]);
            mat_mul.mShares[0](i, 0) = resRight[i].mData[0] - resLeft[i].mData[0];
            mat_mul.mShares[1](i, 0) = resRight[i].mData[1] - resLeft[i].mData[1];
        }
        cipher_mul_seq(this->pIdx, mat_mul, comp_res, mat_res, *(this->eval), *(this->enc), *(this->runtime));
        
        // free the space of mat_mul and comp_res.
        comp_res.resize(0, 0);
        mat_mul.resize(0, 0);

        // compute the final result
        for(int i=0; i<resLeft.size(); i++){
            local_res[i].mData[0] = resLeft[i].mData[0] + mat_res.mShares[0](i, 0);
            local_res[i].mData[1] = resLeft[i].mData[1] + mat_res.mShares[1](i, 0);            
        }
        return;
    }

    std::pair<std::vector<std::vector<aby3::si64>>, std::vector<std::vector<aby3::si64>>> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::pair<std::vector<std::vector<aby3::si64>>, std::vector<std::vector<aby3::si64>>> data_loading(size_t k=1){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecX = generate_vector_si64(this->n * k, this->pIdx, *this->enc, *this->runtime);
        std::vector<aby3::si64> vecY = generate_vector_si64(partial_len * k, this->pIdx, *this->enc, *this->runtime);

        std::vector<std::vector<aby3::si64>> vecX2d(this->n, std::vector<aby3::si64>(k));
        std::vector<std::vector<aby3::si64>> vecY2d(partial_len, std::vector<aby3::si64>(k));
        for(size_t i=0; i<this->n; i++){
            for(size_t j=0; j<k; j++){
                vecX2d[i][j] = vecX[i*k + j];
            }
        }
        for(size_t i=0; i<partial_len; i++){
            for(size_t j=0; j<k; j++){
                vecY2d[i][j] = vecY[i*k + j];
            }
        }
        return std::make_pair(vecX2d, vecY2d);
    }
};

// Search function class.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class PtABS : public SubTask<NUMX, NUMY, NUMT, NUMR> {
public:
    //aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;
    int key_bitsize = 32;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    PtABS(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = true;
            this->lookahead = 1;
            this->lookahead_axis = 0;
    }

    void set_key_bitsize(int key_bitsize){
        if(key_bitsize < 1 || key_bitsize > 53) THROW_RUNTIME_ERROR("Invalid key bitsize.");
        this->key_bitsize = key_bitsize;
    }

    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        aby3::u64 block_length = binfo->block_len;
        aby3::u64 valid_length = binfo->block_len + this->lookahead * this->n;
        if(binfo->t_start + valid_length >= (this->m * this->n)) valid_length = (this->n * this->m) - binfo->t_start;

        // compute the local table.
        aby3::sbMatrix one_more_pieces_table(valid_length, 1);
        aby3::sbMatrix expand_key(valid_length, this->key_bitsize);
        aby3::sbMatrix expand_keyset(valid_length, this->key_bitsize);
        for(size_t i=0; i<valid_length; i++){
            expand_key.mShares[0](i, 0) = expandX[i].mData[0];
            expand_key.mShares[1](i, 0) = expandX[i].mData[1];
            expand_keyset.mShares[0](i, 0) = expandY[i].mData[0];
            expand_keyset.mShares[1](i, 0) = expandY[i].mData[1];
        }

        bool_cipher_lt(this->pIdx, expand_key, expand_keyset, one_more_pieces_table, *(this->enc), *(this->eval), *(this->runtime));
        bool_cipher_not(this->pIdx, one_more_pieces_table, one_more_pieces_table);

        // differential substraction!
        aby3::sbMatrix diff_table(binfo->block_len, 1);
        if(binfo->t_start < (this->m - 1) * this->n){
            for(int i=0; i<valid_length - this->n; i++){
                diff_table.mShares[0](i, 0) = one_more_pieces_table.mShares[0](i+this->n, 0) ^ one_more_pieces_table.mShares[0](i, 0);
                diff_table.mShares[1](i, 0) = one_more_pieces_table.mShares[1](i+this->n, 0) ^ one_more_pieces_table.mShares[1](i, 0);
            }
            for(int i=valid_length - this->n; i<binfo->block_len; i++){
                diff_table.mShares[0](i, 0) = one_more_pieces_table.mShares[0](i, 0);
                diff_table.mShares[1](i, 0) = one_more_pieces_table.mShares[1](i, 0);
            }
        }
        else{
            for(int i=0; i<binfo->block_len; i++){
                diff_table.mShares[0](i, 0) = one_more_pieces_table.mShares[0](i+this->n, 0) ^ one_more_pieces_table.mShares[0](i, 0);
                diff_table.mShares[1](i, 0) = one_more_pieces_table.mShares[1](i+this->n, 0) ^ one_more_pieces_table.mShares[1](i, 0);
            }
        }

        // reduce to local table.
        aby3::si64Matrix expandV(binfo->block_len, 1);
        for(size_t i=0; i<binfo->block_len; i++){
            expandV.mShares[0](i, 0) = this->selectV[binfo->t_start + i].mData[0];
            expandV.mShares[1](i, 0) = this->selectV[binfo->t_start + i].mData[1];
        }

        // cipher mul.
        aby3::si64Matrix table_res(binfo->block_len, 1);
        (*this->eval).asyncMul((*this->runtime), expandV, diff_table, table_res).get();
        for(size_t i=0; i<binfo->block_len; i++){
            local_table[i].mData[0] = table_res.mShares[0](i, 0);
            local_table[i].mData[1] = table_res.mShares[1](i, 0);
        }

        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        for(size_t i=0; i<resLeft.size(); i++){
            local_res[i] = resLeft[i] + resRight[i];
        }
        return;
    }
};


template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class PtAMBS : public PairOnlySubTask<NUMX, NUMY, NUMT, NUMR> {
public:
    //aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;
    int key_bitsize = 32;

    using PairOnlySubTask<NUMX, NUMY, NUMT, NUMR>::PairOnlySubTask;

    PtAMBS(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        PairOnlySubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
            this->lookahead = 1;
            this->lookahead_axis = 0;
    }

    void set_key_bitsize(int key_bitsize){
        if(key_bitsize < 1 || key_bitsize > 53) THROW_RUNTIME_ERROR("Invalid key bitsize.");
        this->key_bitsize = key_bitsize;
    }

    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {

        aby3::u64 block_length = binfo->block_len;
        aby3::u64 valid_length = binfo->block_len + this->lookahead * this->n;
        if(binfo->t_start + valid_length >= (this->m * this->n)) valid_length = (this->n * this->m) - binfo->t_start;

        // compute the local table.
        aby3::sbMatrix one_more_pieces_table(valid_length, 1);
        aby3::sbMatrix expand_key(valid_length, this->key_bitsize);
        aby3::sbMatrix expand_keyset(valid_length, this->key_bitsize);
        for(size_t i=0; i<valid_length; i++){
            expand_key.mShares[0](i, 0) = expandX[i].mData[0];
            expand_key.mShares[1](i, 0) = expandX[i].mData[1];
            expand_keyset.mShares[0](i, 0) = expandY[i].mData[0];
            expand_keyset.mShares[1](i, 0) = expandY[i].mData[1];
        }

        bool_cipher_lt(this->pIdx, expand_key, expand_keyset, one_more_pieces_table, *(this->enc), *(this->eval), *(this->runtime));
        bool_cipher_not(this->pIdx, one_more_pieces_table, one_more_pieces_table);

        // differential substraction!
        aby3::sbMatrix diff_table(binfo->block_len, 1);
        if(binfo->t_start < (this->m - 1) * this->n){
            for(int i=0; i<valid_length - this->n; i++){
                diff_table.mShares[0](i, 0) = one_more_pieces_table.mShares[0](i+this->n, 0) ^ one_more_pieces_table.mShares[0](i, 0);
                diff_table.mShares[1](i, 0) = one_more_pieces_table.mShares[1](i+this->n, 0) ^ one_more_pieces_table.mShares[1](i, 0);
            }
            for(int i=valid_length - this->n; i<binfo->block_len; i++){
                diff_table.mShares[0](i, 0) = one_more_pieces_table.mShares[0](i, 0);
                diff_table.mShares[1](i, 0) = one_more_pieces_table.mShares[1](i, 0);
            }
        }
        else{
            for(int i=0; i<binfo->block_len; i++){
                diff_table.mShares[0](i, 0) = one_more_pieces_table.mShares[0](i+this->n, 0) ^ one_more_pieces_table.mShares[0](i, 0);
                diff_table.mShares[1](i, 0) = one_more_pieces_table.mShares[1](i+this->n, 0) ^ one_more_pieces_table.mShares[1](i, 0);
            }
        }

        for(size_t i=0; i<binfo->block_len; i++){
            local_table[i].mData[0] = diff_table.mShares[0](i, 0);
            local_table[i].mData[1] = diff_table.mShares[1](i, 0);
        }

        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        return;
    }
};

