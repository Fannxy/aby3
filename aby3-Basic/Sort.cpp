#include "Sort.h"
#include <numeric>

static const size_t MAX_SENDING_SIZE = 1 << 25;

using namespace oc;
using namespace aby3;

int bc_sort_different(std::vector<aby3::sbMatrix> &data, std::vector<size_t> &lows, std::vector<size_t> &highs, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t max_size){

    size_t BITSIZE = data[0].bitCount();

    size_t num_interval = lows.size();
    size_t total_size = 0;
    for(size_t i=0; i<num_interval; i++){
        total_size += (highs[i] - lows[i]) * (highs[i] - lows[i]);
    }
    size_t valid_size = std::min(total_size, max_size);


    // Prepare the tiling & repeating indices.
    std::vector<size_t> index_tile(valid_size, 0);
    std::vector<size_t> index_repeat(valid_size, 0);

    // Setup the recorders.
    size_t start = 0;
    size_t total_current = 0;
    size_t pos_current = 0;

    // Prepare the data location exchanging vectors.
    size_t pos_num = std::accumulate(highs.begin(), highs.end(), 0) - std::accumulate(lows.begin(), lows.end(), 0);

    std::vector<size_t> data_pos(pos_num, 0);
    std::vector<size_t> origion_pos(pos_num, 0);

    while(start < num_interval){
        size_t current = total_current;
        for(size_t i=0; i<(num_interval - start); i++){
            size_t low = lows[start + i], high = highs[start + i];
            size_t after = current + (high - low) * (high - low);

            // if(pIdx == 0) debug_info("current after size = " + std::to_string((after - total_current)) + " | after = " + std::to_string(after) + " | total_current = " + std::to_string(total_current));
            // if(pIdx == 0) debug_info("i = " + std::to_string(i) + " | low = " + std::to_string(low) + " | high = " + std::to_string(high) + " | current = " + std::to_string(current));

            if((after - total_current) > valid_size) break;

            // tiling the indices.
            std::vector<size_t> ranging_indices(high - low, 0);
            std::iota(ranging_indices.begin(), ranging_indices.end(), low);

            // std::string ranging_indices_str = "";
            // for(size_t j=0; j<(high - low); j++){
            //     ranging_indices_str += std::to_string(ranging_indices[j]) + " ";
            // }    

            vector_tile<size_t>(ranging_indices, (high - low), index_tile, current - total_current);
            vector_repeat<size_t>(ranging_indices, (high - low), index_repeat, current - total_current);
            current = after;
        }

        // if(pIdx == 0){
        //     std::string ranging_indices_str = "";
        //     for(size_t i=0; i<index_tile.size(); i++){
        //         ranging_indices_str += std::to_string(index_tile[i]) + " ";
        //     }
        // }

        // if(pIdx == 0) debug_info("bc sorts - before tiling & repeating data construction.");

        // vectorized comparison.
        sbMatrix data_tile(current - total_current, BITSIZE);
        sbMatrix data_repeat(current - total_current, BITSIZE);
        sbMatrix data_comp(current - total_current, 1);
        for(size_t i=0; i<(current - total_current); i++){
            data_tile.mShares[0](i, 0) = data[index_tile[i]].mShares[0](0, 0);
            data_tile.mShares[1](i, 0) = data[index_tile[i]].mShares[1](0, 0);
            data_repeat.mShares[0](i, 0) = data[index_repeat[i]].mShares[0](0, 0);
            data_repeat.mShares[1](i, 0) = data[index_repeat[i]].mShares[1](0, 0);
        }

        bool_cipher_lt(pIdx, data_tile, data_repeat, data_comp, enc, eval, runtime);
        i64Matrix plain_comp(current - total_current, 1);
        enc.revealAll(runtime, data_comp, plain_comp).get();

        // compute the location after sorting and construct the data.
        current = total_current;
        size_t next_start = start;

        for(size_t i=0; i<(num_interval - start); i++){
            size_t low = lows[start + i], high = highs[start + i];
            size_t after = current + (high - low) * (high - low);

            if((after - total_current) > valid_size){
                next_start = start + i;
                break;
            }

            next_start += 1;
            size_t pos_after = pos_current + high - low;
            std::vector<size_t> rela_pos(high - low, low);
            for(size_t i=0; i<high - low; i++){
                for(size_t j=0; j<high - low; j++){
                    rela_pos[i] += (size_t) plain_comp(current - total_current + i * (high - low) + j, 0);
                }
            }

            for(size_t i=pos_current; i<pos_after; i++){
                data_pos[i] = rela_pos[i - pos_current];
                origion_pos[i] = i - pos_current + low;
            }

            pos_current = pos_after;
            current = after;
        }

        start = next_start;
        total_current = current;
    }

    // sort the data.
    std::vector<sbMatrix> tmp_data = data;
    for(size_t i=0; i<pos_num; i++){
        data[data_pos[i]] = tmp_data[origion_pos[i]];
    }   

    return 0;
}

int bc_sort_different(aby3::si64Matrix &data, std::vector<size_t> &lows, std::vector<size_t> &highs, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t max_size){

    size_t num_interval = lows.size();
    size_t total_size = 0;
    for(size_t i=0; i<num_interval; i++){
        total_size += (highs[i] - lows[i]) * (highs[i] - lows[i]);
    }
    size_t valid_size = std::min(total_size, max_size);

    // Prepare the tiling & repeating indices.
    std::vector<size_t> index_tile(valid_size, 0);
    std::vector<size_t> index_repeat(valid_size, 0);

    // Setup the recorders.
    size_t start = 0;
    size_t total_current = 0;
    size_t pos_current = 0;

    // Prepare the data location exchanging vectors.
    size_t pos_num = std::accumulate(highs.begin(), highs.end(), 0) - std::accumulate(lows.begin(), lows.end(), 0);

    std::vector<size_t> data_pos(pos_num, 0);
    std::vector<size_t> origion_pos(pos_num, 0);

    while(start < num_interval){
        size_t current = total_current;
        for(size_t i=0; i<(num_interval - start); i++){
            size_t low = lows[start + i], high = highs[start + i];
            size_t after = current + (high - low) * (high - low);

            if((after - total_current) > valid_size) break;

            // tiling the indices.
            std::vector<size_t> ranging_indices(high - low, 0);
            std::iota(ranging_indices.begin(), ranging_indices.end(), low);

            vector_tile<size_t>(ranging_indices, (high - low), index_tile, current - total_current);
            vector_repeat<size_t>(ranging_indices, (high - low), index_repeat, current - total_current);
            current = after;
        }

        // vectorized comparison.
        si64Matrix data_tile(current - total_current, 1);
        si64Matrix data_repeat(current - total_current, 1);
        sbMatrix data_comp(current - total_current, 1);
        for(size_t i=0; i<(current - total_current); i++){
            data_tile.mShares[0](i, 0) = data.mShares[0](index_tile[i], 0);
            data_tile.mShares[1](i, 0) = data.mShares[1](index_tile[i], 0);
            data_repeat.mShares[0](i, 0) = data.mShares[0](index_repeat[i], 0);
            data_repeat.mShares[1](i, 0) = data.mShares[1](index_repeat[i], 0);
        }

        arith_cipher_lt(pIdx, data_tile, data_repeat, data_comp, enc, eval, runtime);
        i64Matrix plain_comp(current - total_current, 1);
        enc.revealAll(runtime, data_comp, plain_comp).get();

        // compute the location after sorting and construct the data.
        current = total_current;
        size_t next_start = start;

        for(size_t i=0; i<(num_interval - start); i++){
            size_t low = lows[start + i], high = highs[start + i];
            size_t after = current + (high - low) * (high - low);

            if((after - total_current) > valid_size){
                next_start = start + i;
                break;
            }

            next_start += 1;
            size_t pos_after = pos_current + high - low;
            std::vector<size_t> rela_pos(high - low, low);
            for(size_t i=0; i<high - low; i++){
                for(size_t j=0; j<high - low; j++){
                    rela_pos[i] += (size_t) plain_comp(current - total_current + i * (high - low) + j, 0);
                }
            }

            for(size_t i=pos_current; i<pos_after; i++){
                data_pos[i] = rela_pos[i - pos_current];
                origion_pos[i] = i - pos_current + low;
            }

            pos_current = pos_after;
            current = after;
        }

        start = next_start;
        total_current = current;

    }

    // sort the data.
    si64Matrix tmp_data = data;
    for(size_t i=0; i<pos_num; i++){
        data.mShares[0](data_pos[i], 0) = tmp_data.mShares[0](origion_pos[i], 0);
        data.mShares[1](data_pos[i], 0) = tmp_data.mShares[1](origion_pos[i], 0);
    }

    return 0;
}


int quick_sort_different(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size){
    size_t BITSIZE = data[0].bitCount();
    size_t len = data.size();

    std::vector<size_t> lows = {0};
    std::vector<size_t> highs = {len};
    size_t num_interval = 1;

    std::vector<size_t> bc_lows;
    std::vector<size_t> bc_highs;

    if(min_size < 4) min_size = 4;

    while(num_interval > 0){
        std::vector<size_t> mids(lows.size());
        for(size_t i=0; i<lows.size(); i++){
            mids[i] = ((highs[i] - 1 - lows[i]) >> 1) + lows[i];
        }
        
        // Select out the pivot candidates for each 
        std::vector<sbMatrix> pivot_candidates(num_interval * 3);
        for(size_t i=0; i<num_interval; i++){
            pivot_candidates[i * 3] = data[lows[i]];
            pivot_candidates[i * 3 + 1] = data[mids[i]];
            pivot_candidates[i * 3 + 2] = data[highs[i] - 1];
        }

        sbMatrix tile_candidate_pivots(num_interval * 3 * 3, BITSIZE);
        sbMatrix repeat_candidate_pivots(num_interval * 3 * 3, BITSIZE);

        for(size_t i=0; i < num_interval; i++){
            size_t start_ind = i * 3 * 3;
            for(size_t j=0; j<3; j++){
                for(size_t k=0; k<3; k++){
                    tile_candidate_pivots.mShares[0](start_ind + j*3 + k, 0) = pivot_candidates[i * 3 + k].mShares[0](0, 0);
                    tile_candidate_pivots.mShares[1](start_ind + j*3 + k, 0) = pivot_candidates[i * 3 + k].mShares[1](0, 0);
                    repeat_candidate_pivots.mShares[0](start_ind + j*3 + k, 0) = pivot_candidates[i * 3 + j].mShares[0](0, 0);
                    repeat_candidate_pivots.mShares[1](start_ind + j*3 + k, 0) = pivot_candidates[i * 3 + j].mShares[1](0, 0);
                }
            }
        }

        // compare the pivots in each interval through one vectorized lt and then reveal it.
        sbMatrix comp_candidate_pivots(num_interval * 3 * 3, 1);
        bool_cipher_lt(pIdx, tile_candidate_pivots, repeat_candidate_pivots, comp_candidate_pivots, enc, eval, runtime);
        i64Matrix plain_comp_candidates(num_interval * 3 * 3, 1);
        enc.revealAll(runtime, comp_candidate_pivots, plain_comp_candidates).get();

        // compute the indices of each interval pivots.
        std::vector<std::vector<int>> comp_sum(num_interval);
        for(size_t i=0; i<num_interval; i++){
            comp_sum[i].resize(3);
            std::fill(comp_sum[i].begin(), comp_sum[i].begin() + 3, 0);
            for(size_t j=0; j<3; j++){
                for(size_t k=0; k<3; k++) comp_sum[i][j] += (int) plain_comp_candidates(i * 3*3 + j * 3 + k, 0);
            }
        }

        // compute the indices of each interval pivots.
        std::vector<size_t> which_pivots = argwhere(comp_sum, 1);
        
        // swicth the elements of each interval, put the pivot in the first.
        std::vector<size_t> switch_mid = argwhere(which_pivots, 1);
        for(size_t i=0; i<switch_mid.size(); i++){
            size_t si = switch_mid[i];
            data[lows[si]] = pivot_candidates[si * 3 + 1];
            data[mids[si]] = pivot_candidates[si * 3];
        }

        std::vector<size_t> switch_high = argwhere(which_pivots, 2);
        for(size_t i=0; i<switch_high.size(); i++){
            size_t si = switch_high[i];
            data[lows[si]] = pivot_candidates[si * 3 + 2];
            data[highs[si] - 1] = pivot_candidates[si * 3];
        }

        // vectorized compare the elements with the midpoint of each interval.
        size_t total_cmp_size = std::accumulate(highs.begin(), highs.end(), 0) - std::accumulate(lows.begin(), lows.end(), 0) - num_interval;
        std::vector<size_t> to_cmp_with_low_pos(total_cmp_size, 0);
        std::vector<size_t> to_cmp_low_pos(total_cmp_size, 0);

        size_t current = 0;
        for(size_t i=0; i<num_interval; i++){
            size_t low = lows[i], high = highs[i];
            size_t after = current + high - low - 1;

            for(size_t j=current; j<after; j++){
                to_cmp_with_low_pos[j] = low + 1 + j - current;
                to_cmp_low_pos[j] = low;
            }
            current = after;
        }


        sbMatrix to_cmp_low(total_cmp_size, BITSIZE);
        sbMatrix to_cmp_with_low(total_cmp_size, BITSIZE);
        sbMatrix tmp_cmp(total_cmp_size, 1);
        i64Matrix plain_cmp(total_cmp_size, 1);

        for(size_t i=0; i<total_cmp_size; i++){
            to_cmp_low.mShares[0](i, 0) = data[to_cmp_low_pos[i]].mShares[0](0, 0);
            to_cmp_low.mShares[1](i, 0) = data[to_cmp_low_pos[i]].mShares[1](0, 0);
            to_cmp_with_low.mShares[0](i, 0) = data[to_cmp_with_low_pos[i]].mShares[0](0, 0);
            to_cmp_with_low.mShares[1](i, 0) = data[to_cmp_with_low_pos[i]].mShares[1](0, 0);
        }
        bool_cipher_lt(pIdx, to_cmp_low, to_cmp_with_low, tmp_cmp, enc, eval, runtime);
        enc.revealAll(runtime, tmp_cmp, plain_cmp).get();

        // update the intervals.
        std::vector<size_t> new_lows;
        std::vector<size_t> new_highs;
        std::vector<size_t> data_pos(current + num_interval, 0);
        std::vector<size_t> origion_pos(current + num_interval, 0);

        // construct the corresponsing positions.
        current = 0;
        size_t new_current = 0;
        for(size_t i=0; i<num_interval; i++){
            size_t low = lows[i], high = highs[i];
            size_t after = current + high - low - 1;
            size_t new_after = new_current + high - low;

            size_t large_num = 0;
            for(size_t j=current; j<after; j++){
                large_num += (size_t) plain_cmp(j, 0);
            }
            size_t low_pos = high - 1 - large_num;
            for(size_t j=new_current; j<new_after; j++) data_pos[j] = low + j - new_current;

            aby3::i64Matrix cmp_interval = plain_cmp.block(current, 0, after - current, 1);

            std::vector<size_t> tmp_lower_args = argwhere(cmp_interval, 0);
            std::vector<size_t> tmp_higher_args = argwhere(cmp_interval, 1);

            for(size_t j=new_current; j<new_after; j++){
                if(j < new_after - large_num - 1){
                    origion_pos[j] = low + 1 + tmp_lower_args[j - new_current];
                }
                else if(j > new_after - large_num - 1){
                    origion_pos[j] = low + 1 + tmp_higher_args[j - (new_after - large_num)];
                }
                else{
                    origion_pos[j] = low;
                }
            }

            // construct the following interval end-points.
            if((low_pos - low) > min_size){
                new_lows.push_back(low);
                new_highs.push_back(low_pos);
            }
            else{
                if((low_pos - low) > 1){
                    bc_lows.push_back(low);
                    bc_highs.push_back(low_pos);
                }
            }

            if((high - low_pos - 1) > min_size){
                new_lows.push_back(low_pos + 1);
                new_highs.push_back(high);
            }
            else{
                if((high - low_pos - 1) > 1) {
                    bc_lows.push_back(low_pos + 1);
                    bc_highs.push_back(high);
                }
            }

            current = after;
            new_current = new_after;
        }
        
        std::vector<sbMatrix> tmp_data = data;
        for(size_t i=0; i<current + num_interval; i++){
            data[data_pos[i]] = tmp_data[origion_pos[i]];
        }

        lows = new_lows;
        highs = new_highs;
        num_interval = lows.size();
    }

    bc_sort_different(data, bc_lows, bc_highs, pIdx, enc, eval, runtime, 1048576);
    return 0;
}

int quick_sort(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size){

    tag_append(pIdx, data);
    efficient_shuffle(data, pIdx, data, enc, eval, runtime);
    quick_sort_different(data, pIdx, enc, eval, runtime, min_size);
    size_t tag_size = std::ceil(std::log2(data.size()));
    tag_remove(pIdx, tag_size, data);

    return 0;
}

int quick_sort_different(aby3::si64Matrix& data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size){

    size_t len = data.rows(); // data is in one-dimensional.
    std::vector<size_t> lows = {0};
    std::vector<size_t> highs = {len};
    size_t num_interval = 1;

    std::vector<size_t> bc_lows;
    std::vector<size_t> bc_highs;

    if(min_size < 4) min_size = 4;

    while(num_interval > 0){
        std::vector<size_t> mids(lows.size());
        for(size_t i=0; i<lows.size(); i++){
            mids[i] = ((highs[i] - 1 - lows[i]) >> 1) + lows[i];
        }

        // Select out the pivot candidates for each
        si64Matrix pivot_candidates(num_interval * 3, 1);
        for(size_t i=0; i<num_interval; i++){
            // pivot_candidates(i * 3, 0) = data(lows[i], 0);
            pivot_candidates.mShares[0](i * 3, 0) = data.mShares[0](lows[i], 0);
            pivot_candidates.mShares[1](i * 3, 0) = data.mShares[1](lows[i], 0);
            // pivot_candidates(i * 3 + 1, 0) = data(mids[i], 0);
            pivot_candidates.mShares[0](i * 3 + 1, 0) = data.mShares[0](mids[i], 0);
            pivot_candidates.mShares[1](i * 3 + 1, 0) = data.mShares[1](mids[i], 0);
            // pivot_candidates(i * 3 + 2, 0) = data(highs[i] - 1, 0);
            pivot_candidates.mShares[0](i * 3 + 2, 0) = data.mShares[0](highs[i] - 1, 0);
            pivot_candidates.mShares[1](i * 3 + 2, 0) = data.mShares[1](highs[i] - 1, 0);
        }

        si64Matrix tile_candidate_pivots(num_interval * 3 * 3, 1);
        si64Matrix repeat_candidate_pivots(num_interval * 3 * 3, 1);

        for(size_t i=0; i < num_interval; i++){
            size_t start_ind = i * 3 * 3;
            for(size_t j=0; j<3; j++){
                for(size_t k=0; k<3; k++){
                    tile_candidate_pivots.mShares[0](start_ind + j*3 + k, 0) = pivot_candidates.mShares[0](i * 3 + k, 0);
                    tile_candidate_pivots.mShares[1](start_ind + j*3 + k, 0) = pivot_candidates.mShares[1](i * 3 + k, 0);
                    repeat_candidate_pivots.mShares[0](start_ind + j*3 + k, 0) = pivot_candidates.mShares[0](i * 3 + j, 0);
                    repeat_candidate_pivots.mShares[1](start_ind + j*3 + k, 0) = pivot_candidates.mShares[1](i * 3 + j, 0);
                }
            }
        }

        // compare the pivots in each interval through one vectorized lt and then reveal it.
        sbMatrix comp_candidate_pivots(num_interval * 3 * 3, 1);
        arith_cipher_lt(pIdx, repeat_candidate_pivots, tile_candidate_pivots, comp_candidate_pivots, enc, eval, runtime);
        i64Matrix plain_comp_candidates(num_interval * 3 * 3, 1);
        enc.revealAll(runtime, comp_candidate_pivots, plain_comp_candidates).get();

        // compute the indices of each interval pivots.
        std::vector<std::vector<int>> comp_sum(num_interval);
        for(size_t i=0; i<num_interval; i++){
            comp_sum[i].resize(3);
            std::fill(comp_sum[i].begin(), comp_sum[i].begin() + 3, 0);
            for(size_t j=0; j<3; j++){
                for(size_t k=0; k<3; k++) comp_sum[i][j] += (int) plain_comp_candidates(i * 3*3 + j * 3 + k, 0);
            }
        }

        // compute the indices of each interval pivots.
        std::vector<size_t> which_pivots = argwhere(comp_sum, 1);
        
        // swicth the elements of each interval, put the pivot in the first.
        std::vector<size_t> switch_mid = argwhere(which_pivots, 1);
        for(size_t i=0; i<switch_mid.size(); i++){
            size_t si = switch_mid[i];
            data.mShares[0](lows[si], 0) = pivot_candidates.mShares[0](si * 3 + 1, 0);
            data.mShares[1](lows[si], 0) = pivot_candidates.mShares[1](si * 3 + 1, 0);
            data.mShares[0](mids[si], 0) = pivot_candidates.mShares[0](si * 3, 0);
            data.mShares[1](mids[si], 0) = pivot_candidates.mShares[1](si * 3, 0);
            // data[lows[si]] = pivot_candidates[si * 3 + 1];
            // data[mids[si]] = pivot_candidates[si * 3];
        }

        std::vector<size_t> switch_high = argwhere(which_pivots, 2);
        for(size_t i=0; i<switch_high.size(); i++){
            size_t si = switch_high[i];
            // data[lows[si]] = pivot_candidates[si * 3 + 2];
            // data[highs[si] - 1] = pivot_candidates[si * 3];
            data.mShares[0](lows[si], 0) = pivot_candidates.mShares[0](si * 3 + 2, 0);
            data.mShares[1](lows[si], 0) = pivot_candidates.mShares[1](si * 3 + 2, 0);
            data.mShares[0](highs[si] - 1, 0) = pivot_candidates.mShares[0](si * 3, 0);
            data.mShares[1](highs[si] - 1, 0) = pivot_candidates.mShares[1](si * 3, 0);
        }

        // vectorized compare the elements with the midpoint of each interval.
        size_t total_cmp_size = std::accumulate(highs.begin(), highs.end(), 0) - std::accumulate(lows.begin(), lows.end(), 0) - num_interval;
        std::vector<size_t> to_cmp_with_low_pos(total_cmp_size, 0);
        std::vector<size_t> to_cmp_low_pos(total_cmp_size, 0);

        size_t current = 0;
        for(size_t i=0; i<num_interval; i++){
            size_t low = lows[i], high = highs[i];
            size_t after = current + high - low - 1;

            for(size_t j=current; j<after; j++){
                to_cmp_with_low_pos[j] = low + 1 + j - current;
                to_cmp_low_pos[j] = low;
            }
            current = after;
        }

        si64Matrix to_cmp_low(total_cmp_size, 1);
        si64Matrix to_cmp_with_low(total_cmp_size, 1);
        sbMatrix tmp_cmp(total_cmp_size, 1);
        i64Matrix plain_cmp(total_cmp_size, 1);
        for(size_t i=0; i<total_cmp_size; i++){
            to_cmp_low.mShares[0](i, 0) = data.mShares[0](to_cmp_low_pos[i], 0);
            to_cmp_low.mShares[1](i, 0) = data.mShares[1](to_cmp_low_pos[i], 0);
            to_cmp_with_low.mShares[0](i, 0) = data.mShares[0](to_cmp_with_low_pos[i], 0);
            to_cmp_with_low.mShares[1](i, 0) = data.mShares[1](to_cmp_with_low_pos[i], 0);
        }
        arith_cipher_lt(pIdx, to_cmp_low, to_cmp_with_low, tmp_cmp, enc, eval, runtime);
        enc.revealAll(runtime, tmp_cmp, plain_cmp).get();

        // update the intervals.
        std::vector<size_t> new_lows;
        std::vector<size_t> new_highs;
        std::vector<size_t> data_pos(current + num_interval, 0);
        std::vector<size_t> origion_pos(current + num_interval, 0);

        // construct the corresponsing positions.
        current = 0;
        size_t new_current = 0;
        for(size_t i=0; i<num_interval; i++){
            size_t low = lows[i], high = highs[i];
            size_t after = current + high - low - 1;
            size_t new_after = new_current + high - low;

            size_t large_num = 0;
            for(size_t j=current; j<after; j++){
                large_num += (size_t) plain_cmp(j, 0);
            }
            size_t low_pos = high - 1 - large_num;
            for(size_t j=new_current; j<new_after; j++) data_pos[j] = low + j - new_current;

            aby3::i64Matrix cmp_interval = plain_cmp.block(current, 0, after - current, 1);

            std::vector<size_t> tmp_lower_args = argwhere(cmp_interval, 0);
            std::vector<size_t> tmp_higher_args = argwhere(cmp_interval, 1);

            for(size_t j=new_current; j<new_after; j++){
                if(j < new_after - large_num - 1){
                    origion_pos[j] = low + 1 + tmp_lower_args[j - new_current];
                }
                else if(j > new_after - large_num - 1){
                    origion_pos[j] = low + 1 + tmp_higher_args[j - (new_after - large_num)];
                }
                else{
                    origion_pos[j] = low;
                }
            }

            // construct the following interval end-points.
            if((low_pos - low) > min_size){
                new_lows.push_back(low);
                new_highs.push_back(low_pos);
            }
            else{
                if((low_pos - low) > 1){
                    bc_lows.push_back(low);
                    bc_highs.push_back(low_pos);
                }
            }

            if((high - low_pos - 1) > min_size){
                new_lows.push_back(low_pos + 1);
                new_highs.push_back(high);
            }
            else{
                if((high - low_pos - 1) > 1) {
                    bc_lows.push_back(low_pos + 1);
                    bc_highs.push_back(high);
                }
            }

            current = after;
            new_current = new_after;
    }

    si64Matrix tmp_data = data;
    for(size_t i=0; i<current + num_interval; i++){
        data.mShares[0](data_pos[i], 0) = tmp_data.mShares[0](origion_pos[i], 0);
        data.mShares[1](data_pos[i], 0) = tmp_data.mShares[1](origion_pos[i], 0);
    }

    lows = new_lows;
    highs = new_highs;
    num_interval = lows.size();
    }

    bc_sort_different(data, bc_lows, bc_highs, pIdx, enc, eval, runtime, 1048576);

    return 0;
}

int quick_sort(aby3::si64Matrix& data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size){

    aby3::sbMatrix bool_data(data.rows(), 64); // fixed data bit size.
    arith2bool(pIdx, data, bool_data, enc, eval, runtime);
    efficient_shuffle(bool_data, pIdx, bool_data, enc, eval, runtime);
    bool2arith(pIdx, bool_data, data, enc, eval, runtime);
    quick_sort_different(data, pIdx, enc, eval, runtime, min_size);

    return 0;
}

int odd_even_merge(aby3::sbMatrix& data1, aby3::sbMatrix& data2, aby3::sbMatrix& res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    int arr1_length = data1.rows();
    int arr2_length = data2.rows();
    int length = std::max(arr1_length, arr2_length);

    // init the result.
    sbMatrix result(length*2, BITSIZE);
    sbMatrix max_ele, max1, max2;
    boolIndex max_arr1, max_arr2;
    max_arr1.indexShares[0] = data1.mShares[0](arr1_length - 1, 0);
    max_arr1.indexShares[1] = data1.mShares[1](arr1_length - 1, 0);
    max_arr2.indexShares[0] = data2.mShares[0](arr2_length - 1, 0);
    max_arr2.indexShares[1] = data2.mShares[1](arr2_length - 1, 0);

    max1 = max_arr1.to_matrix();
    max2 = max_arr2.to_matrix();

    bool_cipher_max(pIdx, max1, max2, max_ele, enc, eval, runtime);
    
    std::fill(result.mShares[0].data(), result.mShares[0].data() + result.mShares[0].size(), max_ele.mShares[0](0, 0));
    std::fill(result.mShares[1].data(), result.mShares[1].data() + result.mShares[1].size(), max_ele.mShares[1](0, 0));

    // organize the result
    for(size_t i=0; i<arr1_length; i++){
        result.mShares[0](i*2, 0) = data1.mShares[0](i, 0);
        result.mShares[1](i*2, 0) = data1.mShares[1](i, 0);
    }
    for(size_t i=0; i<arr2_length; i++){
        result.mShares[0](i*2 + 1, 0) = data2.mShares[0](i, 0);
        result.mShares[1](i*2 + 1, 0) = data2.mShares[1](i, 0);
    }

    // begin the odd_even merge
    size_t t = (size_t) std::ceil((std::log2((double)length)) + 1);
    size_t q = std::pow(2, t-1);
    size_t d = 1;
    size_t r = 0;

    while(d > 0){
        // get the to-be-compared indexing mask.
        std::vector<int> x_mask, y_mask;
        for(int i=r; i<length*2 - d; i+=2){
            x_mask.push_back(i);
            y_mask.push_back(i + d);
        }

        // get the result into the sbMatrix.
        sbMatrix x_mask_mat(x_mask.size(), BITSIZE);
        sbMatrix y_mask_mat(y_mask.size(), BITSIZE);
        for(int i=0; i<x_mask.size(); i++){
            x_mask_mat.mShares[0](i, 0) = result.mShares[0](x_mask[i], 0);
            x_mask_mat.mShares[1](i, 0) = result.mShares[1](x_mask[i], 0);
            y_mask_mat.mShares[0](i, 0) = result.mShares[0](y_mask[i], 0);
            y_mask_mat.mShares[1](i, 0) = result.mShares[1](y_mask[i], 0);
        }
        sbMatrix max_mat, min_mat;
        bool_cipher_max_min_split(pIdx, x_mask_mat, y_mask_mat, max_mat, min_mat, enc, eval, runtime);
        
        // update the result.
        for(int i=0; i<x_mask.size(); i++){
            result.mShares[0](x_mask[i], 0) = min_mat.mShares[0](i, 0);
            result.mShares[1](x_mask[i], 0) = min_mat.mShares[1](i, 0);
            result.mShares[0](y_mask[i], 0) = max_mat.mShares[0](i, 0);
            result.mShares[1](y_mask[i], 0) = max_mat.mShares[1](i, 0);
        }

        // update the d, r.
        d = (size_t) (q - 1);
        q = q >> 1;
        r = 1;
    }
    res.resize(arr1_length+arr2_length, BITSIZE);
    for(size_t i=0; i<arr1_length+arr2_length; i++){
        res.mShares[0](i, 0) = result.mShares[0](i, 0);
        res.mShares[1](i, 0) = result.mShares[1](i, 0);
    }

    return 0;
}

int odd_even_merge(aby3::si64Matrix& data1, aby3::si64Matrix& data2, aby3::si64Matrix& res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    

    return 0;
}

/*
    The odd_even_multi_merge function is used to merge the data sets in the data vector.
    The data vector contains the data sets to be merged.
    The sorted_res is the result of the merge.
*/
int odd_even_multi_merge(std::vector<aby3::sbMatrix> &data, aby3::sbMatrix& sorted_res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    int k = data.size();

    while(k!=1){
        if(k%2 != 0){
            aby3::sbMatrix res;
            odd_even_merge(data[k-2], data[k-1], res, pIdx, enc, eval, runtime);
            data[k-2] = res;
            k -=1;
        }
        else{
            for(int i=0; i<k; i+=2){
                aby3::sbMatrix res;
                odd_even_merge(data[i], data[i+1], res, pIdx, enc, eval, runtime);
                data[i/2] = res;
            }
            k = k>>1;
        }
    }

    sorted_res = data[0];

    return 0;
}

int high_dimensional_odd_even_merge(std::vector<aby3::sbMatrix> &data1, std::vector<aby3::sbMatrix>& data2, std::vector<aby3::sbMatrix>& sorted_res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    size_t dim = data1.size();
    if(dim != data2.size()){
        THROW_RUNTIME_ERROR("The dimensions of the two data sets are not equal!");
    }

    int arr1_length = data1[0].rows();
    int arr2_length = data2[0].rows();
    int length = std::max(arr1_length, arr2_length);

    // init the result.
    std::vector<sbMatrix> result(dim);
    for(size_t i=0; i<dim; i++){
        result[i].resize(length*2, BITSIZE);
    }

    sbMatrix max_ele, max1, max2;
    max1.resize(dim, BITSIZE);
    max2.resize(dim, BITSIZE);
    max_ele.resize(dim, BITSIZE);
    for(size_t i=0; i<dim; i++){
        max1.mShares[0](i, 0) = data1[i].mShares[0](arr1_length - 1, 0);
        max1.mShares[1](i, 0) = data1[i].mShares[1](arr1_length - 1, 0);
        max2.mShares[0](i, 0) = data2[i].mShares[0](arr2_length - 1, 0);
        max2.mShares[1](i, 0) = data2[i].mShares[1](arr2_length - 1, 0);
    }

    size_t round = (size_t) ceil(dim / (double) MAX_SENDING_SIZE);
    size_t last_len = dim - (round - 1) * MAX_SENDING_SIZE;

    // if(pIdx == 0) debug_info("round = " + std::to_string(round) + " | last_len = " + std::to_string(last_len));
    for(size_t i=0; i<round; i++){
        size_t unit_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
        sbMatrix max1_part(unit_len, BITSIZE);
        sbMatrix max2_part(unit_len, BITSIZE);
        sbMatrix max_ele_part(unit_len, BITSIZE);

        std::memcpy(max1_part.mShares[0].data(), max1.mShares[0].data() + i * MAX_SENDING_SIZE, unit_len * sizeof(max1.mShares[0](0, 0)));
        std::memcpy(max1_part.mShares[1].data(), max1.mShares[1].data() + i * MAX_SENDING_SIZE, unit_len * sizeof(max1.mShares[1](0, 0)));
        std::memcpy(max2_part.mShares[0].data(), max2.mShares[0].data() + i * MAX_SENDING_SIZE, unit_len * sizeof(max2.mShares[0](0, 0)));
        std::memcpy(max2_part.mShares[1].data(), max2.mShares[1].data() + i * MAX_SENDING_SIZE, unit_len * sizeof(max2.mShares[1](0, 0)));

        bool_cipher_max(pIdx, max1_part, max2_part, max_ele_part, enc, eval, runtime);

        std::memcpy(max_ele.mShares[0].data() + i * MAX_SENDING_SIZE, max_ele_part.mShares[0].data(), unit_len * sizeof(max_ele_part.mShares[0](0, 0)));
        std::memcpy(max_ele.mShares[1].data() + i * MAX_SENDING_SIZE, max_ele_part.mShares[1].data(), unit_len * sizeof(max_ele_part.mShares[1](0, 0)));
    }
    // if(pIdx == 0) debug_info("after bool_max");
    // bool_cipher_max(pIdx, max1, max2, max_ele, enc, eval, runtime);

    for(size_t i=0; i<dim; i++){
        std::fill(result[i].mShares[0].data(), result[i].mShares[0].data() + result[i].mShares[0].size(), max_ele.mShares[0](i, 0));
        std::fill(result[i].mShares[1].data(), result[i].mShares[1].data() + result[i].mShares[1].size(), max_ele.mShares[1](i, 0));
    }

    // organize the result
    for(size_t i=0; i<dim; i++){
        for(size_t k=0; k<arr1_length; k++){
            result[i].mShares[0](k*2, 0) = data1[i].mShares[0](k, 0);
            result[i].mShares[1](k*2, 0) = data1[i].mShares[1](k, 0);
        }
        for(size_t k=0; k<arr2_length; k++){
            result[i].mShares[0](k*2 + 1, 0) = data2[i].mShares[0](k, 0);
            result[i].mShares[1](k*2 + 1, 0) = data2[i].mShares[1](k, 0);
        }
    }

    // begin the high-dimensional odd_even merge sort.
    size_t t = (size_t) std::ceil((std::log2((double)length)) + 1);
    size_t q = std::pow(2, t-1);
    size_t d = 1;
    size_t r = 0;

    while(d > 0){
        std::vector<int> x_mask, y_mask;
        for(int i=r; i<length*2 - d; i+=2){
            x_mask.push_back(i);
            y_mask.push_back(i + d);
        }
        size_t unit_mask_len = x_mask.size();

        // get the result into the sbMatrix.
        sbMatrix x_mask_mat(unit_mask_len * dim, BITSIZE);
        sbMatrix y_mask_mat(unit_mask_len * dim, BITSIZE);
        for(size_t i=0; i<dim; i++){
            for(size_t j=0; j<unit_mask_len; j++){
                x_mask_mat.mShares[0](i*unit_mask_len + j, 0) = result[i].mShares[0](x_mask[j], 0);
                x_mask_mat.mShares[1](i*unit_mask_len + j, 0) = result[i].mShares[1](x_mask[j], 0);
                y_mask_mat.mShares[0](i*unit_mask_len + j, 0) = result[i].mShares[0](y_mask[j], 0);
                y_mask_mat.mShares[1](i*unit_mask_len + j, 0) = result[i].mShares[1](y_mask[j], 0);
            }
        }
        // sbMatrix max_mat, min_mat;
        sbMatrix max_mat(unit_mask_len * dim, BITSIZE);
        sbMatrix min_mat(unit_mask_len * dim, BITSIZE);
        // bool_cipher_max_min_split(pIdx, x_mask_mat, y_mask_mat, max_mat, min_mat, enc, eval, runtime);
        round = (size_t) ceil((unit_mask_len * dim) / (double)MAX_SENDING_SIZE);
        last_len = unit_mask_len * dim - (round - 1) * MAX_SENDING_SIZE;
        for(size_t i=0; i<round; i++){
            size_t unit_len_ = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
            sbMatrix x_mask_part(unit_len_, BITSIZE);
            sbMatrix y_mask_part(unit_len_, BITSIZE);
            sbMatrix max_part(unit_len_, BITSIZE);
            sbMatrix min_part(unit_len_, BITSIZE);

            std::memcpy(x_mask_part.mShares[0].data(), x_mask_mat.mShares[0].data() + i * MAX_SENDING_SIZE, unit_len_ * sizeof(x_mask_mat.mShares[0](0, 0)));
            std::memcpy(x_mask_part.mShares[1].data(), x_mask_mat.mShares[1].data() + i * MAX_SENDING_SIZE, unit_len_ * sizeof(x_mask_mat.mShares[1](0, 0)));
            std::memcpy(y_mask_part.mShares[0].data(), y_mask_mat.mShares[0].data() + i * MAX_SENDING_SIZE, unit_len_ * sizeof(y_mask_mat.mShares[0](0, 0)));
            std::memcpy(y_mask_part.mShares[1].data(), y_mask_mat.mShares[1].data() + i * MAX_SENDING_SIZE, unit_len_ * sizeof(y_mask_mat.mShares[1](0, 0)));

            bool_cipher_max_min_split(pIdx, x_mask_part, y_mask_part, max_part, min_part, enc, eval, runtime);

            std::memcpy(max_mat.mShares[0].data() + i * MAX_SENDING_SIZE, max_part.mShares[0].data(), unit_len_ * sizeof(max_part.mShares[0](0, 0)));
            std::memcpy(max_mat.mShares[1].data() + i * MAX_SENDING_SIZE, max_part.mShares[1].data(), unit_len_ * sizeof(max_part.mShares[1](0, 0)));
            std::memcpy(min_mat.mShares[0].data() + i * MAX_SENDING_SIZE, min_part.mShares[0].data(), unit_len_ * sizeof(min_part.mShares[0](0, 0)));
            std::memcpy(min_mat.mShares[1].data() + i * MAX_SENDING_SIZE, min_part.mShares[1].data(), unit_len_ * sizeof(min_part.mShares[1](0, 0)));
        }

        // update the result.
        for(size_t i=0; i<dim; i++){
            for(size_t j=0; j<unit_mask_len; j++){
                result[i].mShares[0](x_mask[j], 0) = min_mat.mShares[0](i*unit_mask_len + j, 0);
                result[i].mShares[1](x_mask[j], 0) = min_mat.mShares[1](i*unit_mask_len + j, 0);
                result[i].mShares[0](y_mask[j], 0) = max_mat.mShares[0](i*unit_mask_len + j, 0);
                result[i].mShares[1](y_mask[j], 0) = max_mat.mShares[1](i*unit_mask_len + j, 0);
            }
        }

        // update the d, r.
        d = (size_t) (q - 1);
        q = q >> 1;
        r = 1;
    }

    sorted_res.resize(dim);
    for(size_t i=0; i<dim; i++){
        sorted_res[i].resize(arr1_length+arr2_length, BITSIZE);
        for(size_t j=0; j<arr1_length+arr2_length; j++){
            sorted_res[i].mShares[0](j, 0) = result[i].mShares[0](j, 0);
            sorted_res[i].mShares[1](j, 0) = result[i].mShares[1](j, 0);
        }
    }

    return 0;
}

int high_dimensional_odd_even_multi_merge(std::vector<std::vector<aby3::sbMatrix>> &data,std::vector<aby3::sbMatrix>& sorted_res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    size_t dim = data.size();
    size_t k = data[0].size(); // the number of data sets in each dimension.

    while(k!=1){
        if(k%2 != 0){
            std::vector<aby3::sbMatrix> res;
            std::vector<aby3::sbMatrix> data1(dim), data2(dim);
            for(size_t i=0; i<dim; i++){
                data1[i] = data[i][k-2];
                data2[i] = data[i][k-1];
            }
            high_dimensional_odd_even_merge(data1, data2, res, pIdx, enc, eval, runtime);
            for(size_t i=0; i<dim; i++){
                data[i][k-2] = res[i];
            }
            k -=1;
        }
        else{
            for(size_t i=0; i<k; i+=2){
                std::vector<aby3::sbMatrix> res;
                std::vector<aby3::sbMatrix> data1(dim), data2(dim);
                for(size_t j=0; j<dim; j++){
                    data1[j] = data[j][i];
                    data2[j] = data[j][i+1];
                }
                high_dimensional_odd_even_merge(data1, data2, res, pIdx, enc, eval, runtime);
                for(size_t j=0; j<dim; j++){
                    data[j][i/2] = res[j];
                }
            }
            k = k>>1;
        }
    }

    // sorted_res = data1[0];
    sorted_res.resize(dim);
    for(size_t i=0; i<dim; i++){
        sorted_res[i] = data[i][0];
    }

    return 0;
}