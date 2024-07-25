#pragma once
#include "debug.h"
#include "./Oram/include/oram.h"
#include "Basics.h"
#include "Shuffle.h"
#include "assert.h"

template <typename CIPHERDATA>
class ABY3StashElement : public StashElement<CIPHERDATA> {
   public:
    CIPHERDATA data;
    boolIndex logicalIndex;
    ABY3StashElement(CIPHERDATA data, boolIndex logicalIndex)
        : data(data), logicalIndex(logicalIndex) {}
};

class ABY3PackedIndex : public PackedIndex<boolIndex> {
   public:
    ABY3PackedIndex(size_t pack, boolIndex logicalIndex)
        : PackedIndex<boolIndex>(pack, logicalIndex) {}

    aby3::sbMatrix pack_to_single_matrix() {
        // the size of the wrap_indices is pack + 1.
        aby3::sbMatrix wrap_indices;
        wrap_indices.resize(this->pack + 1, BITSIZE);
        // the first element is the logical index.
        wrap_indices.mShares[0](0, 0) = logicalIndex.indexShares[0];
        wrap_indices.mShares[1](0, 0) = logicalIndex.indexShares[1];
        // all the other elements are the contents.
        for (size_t i = 0; i < this->pack; i++) {
            wrap_indices.mShares[0](i + 1, 0) = packedIndices[i].indexShares[0];
            wrap_indices.mShares[1](i + 1, 0) = packedIndices[i].indexShares[1];
        }
        return wrap_indices;
    }

    void unpack_from_single_matrix(aby3::sbMatrix &wrap_indices) {
        // the size of the wrap_indices is pack + 1.
        this->pack = wrap_indices.rows() - 1;
        this->packedIndices.resize(this->pack);
        // the first element is the logical index.
        logicalIndex.indexShares[0] = wrap_indices.mShares[0](0, 0);
        logicalIndex.indexShares[1] = wrap_indices.mShares[1](0, 0);
        // all the other elements are the contents.
        for (size_t i = 0; i < this->pack; i++) {
            packedIndices[i].indexShares[0] = wrap_indices.mShares[0](i + 1, 0);
            packedIndices[i].indexShares[1] = wrap_indices.mShares[1](i + 1, 0);
        }
    }

    void to_string(std::ofstream &ofs, int pIdx, aby3::Sh3Encryptor &enc,
                   aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime) {
        aby3::sbMatrix data = this->pack_to_single_matrix();
        aby3::i64Matrix plain_data = back2plain(pIdx, data, enc, eval, runtime);
        ofs << std::to_string(plain_data(0, 0)) + ": {";
        for (size_t i = 1; i < this->pack + 1; i++) {
            ofs << std::to_string(plain_data(i, 0)) + ", ";
        }
        ofs << " } ";
    }
};

class ABY3PosMap : public PosMap<boolIndex> {
   public:
    std::vector<boolShare> usage_map;
    ABY3PosMap *subPosMap;
    std::vector<ABY3PackedIndex> packed_index;
    std::vector<ABY3PackedIndex> stash;
    int pIdx;
    aby3::Sh3Encryptor *enc;
    aby3::Sh3Evaluator *eval;
    aby3::Sh3Runtime *runtime;

    ABY3PosMap(size_t len, size_t pack, size_t S,
               std::vector<boolIndex> &permutation, int pIdx,
               aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
               aby3::Sh3Runtime &runtime)
        : PosMap<boolIndex>(len, pack, S),
          pIdx(pIdx),
          enc(&enc),
          runtime(&runtime),
          eval(&eval) {
        if (!checkPowerOfTwo(pack))
            THROW_RUNTIME_ERROR("pack = " + std::to_string(pack) +
                                " must be a power of 2.");

        if (linear) {  // construct the linear map.
            this->usage_map.resize(this->n);
            for (int i = 0; i < this->n; i++)
                bool_init_false(this->pIdx, this->usage_map[i]);
            this->permutation = permutation;
        } else {  // recursively constructe the posMap.

            // 1. construct the packed posMap.
            for (size_t i = 0; i < this->map_len; i++) {
                std::vector<boolIndex> tmp;
                for (size_t j = 0; j < this->pack; j++) {
                    tmp.push_back(permutation[i * this->pack + j]);
                }
                boolIndex tmp_index(i, this->pIdx);
                ABY3PackedIndex tmp_packed_index(this->pack, tmp_index);
                tmp_packed_index.set_packed_indices(tmp);
                this->packed_index.push_back(tmp_packed_index);
            }

            // 2. shuffle the packed posMap, store the permutaion into Pi.
            std::vector<aby3::sbMatrix> _permutation(this->map_len);
            for (size_t i = 0; i < this->map_len; i++) {
                _permutation[i] = this->packed_index[i].pack_to_single_matrix();
            }

            std::vector<aby3::si64> Pi(this->map_len);
            efficient_shuffle_with_random_permutation(
                _permutation, this->pIdx, _permutation, Pi, *(this->enc),
                *(this->eval), *(this->runtime));

            for (size_t i = 0; i < this->map_len; i++) {
                this->packed_index[i].unpack_from_single_matrix(
                    _permutation[i]);
            }

            // 3. construct the subPosMap.
            std::vector<boolIndex> sub_permutation(this->map_len);
            for (size_t i = 0; i < this->map_len; i++) {
                sub_permutation[i] = boolIndex(Pi[i].mData[0], Pi[i].mData[1]);
            }
            this->subPosMap = new ABY3PosMap(
                this->map_len, this->pack, this->S, sub_permutation, this->pIdx,
                *(this->enc), *(this->eval), *(this->runtime));
        }
    }

    aby3::i64 access(boolIndex index, boolShare fake) {
        boolIndex physical_index;

        if (linear) {
            // construct the s1_map.
            aby3::sbMatrix s1_map(this->n, 1);
            aby3::sbMatrix _tmp_index(this->n, 1);
            aby3::i64Matrix _tmp_range(this->n, 1);

            for (size_t i = 0; i < this->n; i++) {
                _tmp_index.mShares[0](i, 0) = index.indexShares[0];
                _tmp_index.mShares[1](i, 0) = index.indexShares[1];
                _tmp_range(i, 0) = i;
            }
            _tmp_index.resize(this->n, BITSIZE);
            bool_cipher_eq(this->pIdx, _tmp_index, _tmp_range, s1_map,
                           *this->enc, *this->eval, *this->runtime);

            // construct the s2_map.
            aby3::sbMatrix s2_map(this->n, 1);
            bool_get_first_zero_mask(this->pIdx, this->usage_map, s2_map,
                                     *this->enc, *this->eval, *this->runtime);

            // expand the length of s1_map and s2_map.
            s1_map.resize(this->n, BITSIZE);
            s2_map.resize(this->n, BITSIZE);

            for (size_t i = 0; i < this->n; i++) {
                s1_map.mShares[0](i, 0) =
                    (s1_map.mShares[0](i, 0) == 1) ? -1 : 0;
                s1_map.mShares[1](i, 0) =
                    (s1_map.mShares[1](i, 0) == 1) ? -1 : 0;
                s2_map.mShares[0](i, 0) =
                    (s2_map.mShares[0](i, 0) == 1) ? -1 : 0;
                s2_map.mShares[1](i, 0) =
                    (s2_map.mShares[1](i, 0) == 1) ? -1 : 0;
            }

            aby3::sbMatrix _permutation =
                vecBoolIndices(this->permutation).to_matrix();

            // get the secret shares of the physical index.
            aby3::sbMatrix _res_s1(1, 1);
            aby3::sbMatrix _res_s2(1, 1);

            _permutation.resize(this->n, BITSIZE);
            bool_cipher_dot(this->pIdx, _permutation, s1_map, _res_s1,
                            *this->enc, *this->eval, *this->runtime);
            bool_cipher_dot(this->pIdx, _permutation, s2_map, _res_s2,
                            *this->enc, *this->eval, *this->runtime);

            // get the physical index.
            aby3::sbMatrix _res(1, BITSIZE);
            bool_cipher_selector(this->pIdx, fake, _res_s2, _res_s1, _res,
                                 *this->enc, *this->eval, *this->runtime);
            physical_index.indexShares[0] = _res.mShares[0](0, 0);
            physical_index.indexShares[1] = _res.mShares[1](0, 0);

            // update use_map, OR(this->use_map, s1_map, this->use_map);
            aby3::sbMatrix _use_map(this->n, 1);
            for (size_t i = 0; i < this->n; i++) {
                _use_map.mShares[0](i, 0) = this->usage_map[i].bshares[0];
                _use_map.mShares[1](i, 0) = this->usage_map[i].bshares[1];
            }
            s1_map.resize(this->n, 1);
            bool_cipher_or(this->pIdx, s1_map, _use_map, _use_map, *this->enc,
                           *this->eval, *this->runtime);
            for (size_t i = 0; i < this->n; i++) {
                this->usage_map[i].from_matrix(_use_map.mShares[0](i, 0),
                                               _use_map.mShares[1](i, 0));
            }
        } else {  // get the recursive indices.

            boolIndex h, l;
            bool_shift_and_left(this->pIdx, index, log2(pack), h, l);

            // boolShare found(false, this->pIdx);
            boolShare found = fake;
            aby3::sbMatrix _phyIndex_in_stash(1, BITSIZE);
            _phyIndex_in_stash = boolIndex(-1, this->pIdx).to_matrix();

            if (this->t > 0) {
                // 1. check whether the stash contains the target element <h>.
                // construct the mask indicating whether the stash contains the
                // target element.
                aby3::sbMatrix expand_h(this->t, BITSIZE);
                aby3::sbMatrix stash_index(this->t, BITSIZE);
                aby3::sbMatrix stash_hit_mask(this->t, 1);
                for (size_t i = 0; i < this->t; i++) {
                    expand_h.mShares[0](i, 0) = h.indexShares[0];
                    expand_h.mShares[1](i, 0) = h.indexShares[1];
                    stash_index.mShares[0](i, 0) =
                        this->stash[i].logicalIndex.indexShares[0];
                    stash_index.mShares[1](i, 0) =
                        this->stash[i].logicalIndex.indexShares[1];
                }
                bool_cipher_eq(this->pIdx, expand_h, stash_index,
                               stash_hit_mask, *this->enc, *this->eval,
                               *this->runtime);

                // update found flag using the above hit_mask.
                boolShare tmp_found =
                    boolShare((bool)(stash_hit_mask.mShares[0](0, 0) & 1),
                              (bool)(stash_hit_mask.mShares[1](0, 0) & 1));
                for (size_t i = 1; i < this->t; i++) {
                    tmp_found.bshares[0] ^=
                        (stash_hit_mask.mShares[0](i, 0) & 1);
                    tmp_found.bshares[1] ^=
                        (stash_hit_mask.mShares[1](i, 0) & 1);
                }
                bool_cipher_or(this->pIdx, found, tmp_found, found, *this->enc,
                               *this->eval, *this->runtime);

                // 2. get the 'physical_index' from the stash.
                // 1) construct the per-pack map, indicating which block
                // contains the target element.
                aby3::sbMatrix expand_l(this->pack, BITSIZE);
                aby3::i64Matrix range_pack(this->pack, 1);
                aby3::sbMatrix per_stash_mask(this->pack, 1);
                for (size_t i = 0; i < this->pack; i++) {
                    expand_l.mShares[0](i, 0) = l.indexShares[0];
                    expand_l.mShares[1](i, 0) = l.indexShares[1];
                    range_pack(i, 0) = i;
                }
                bool_cipher_eq(this->pIdx, expand_l, range_pack, per_stash_mask,
                               *this->enc, *this->eval, *this->runtime);

                // 2) expand the per-pack map and the stash-hit map to (this->t
                // * this->pack) for vectorization.
                expand_h.resize(this->t * this->pack, BITSIZE);
                expand_l.resize(this->t * this->pack, BITSIZE);
                stash_index.resize(this->t * this->pack, BITSIZE);

                for (size_t i = 0; i < this->t; i++) {
                    for (size_t j = 0; j < this->pack; j++) {
                        // whether hit the ith stash element
                        expand_h.mShares[0](i * this->pack + j, 0) =
                            (stash_hit_mask.mShares[0](i) == 1) ? -1 : 0;
                        expand_h.mShares[1](i * this->pack + j, 0) =
                            (stash_hit_mask.mShares[1](i) == 1) ? -1 : 0;

                        // whether hit the jth block in the stash element
                        expand_l.mShares[0](i * this->pack + j, 0) =
                            (per_stash_mask.mShares[0](j) == 1) ? -1 : 0;
                        expand_l.mShares[1](i * this->pack + j, 0) =
                            (per_stash_mask.mShares[1](j) == 1) ? -1 : 0;

                        // construct the expanded stash elements
                        stash_index.mShares[0](i * this->pack + j, 0) =
                            this->stash[i].packedIndices[j].indexShares[0];
                        stash_index.mShares[1](i * this->pack + j, 0) =
                            this->stash[i].packedIndices[j].indexShares[1];
                    }
                }

                // 3) construct the target index map through AND gate.
                aby3::sbMatrix target_index_map(this->t * this->pack, BITSIZE);
                bool_cipher_and(this->pIdx, expand_h, expand_l,
                                target_index_map, *this->enc, *this->eval,
                                *this->runtime);

                // 4) fetch the target index through DOT gate.
                // aby3::sbMatrix _phyIndex_in_stash(1, BITSIZE);
                bool_cipher_dot(this->pIdx, target_index_map, stash_index,
                                _phyIndex_in_stash, *this->enc, *this->eval,
                                *this->runtime);
                boolIndex phyIndex_in_stash(_phyIndex_in_stash);
            }

            // 3. get the physical index in the next level posMap.
            aby3::i64 phyIndex_in_nextMap = this->subPosMap->access(h, found);

            // 4. select the physical index from the stash and the next level
            // posMap.
            boolShare main_flag;
            bool_cipher_not(this->pIdx, found, main_flag);
            bool_cipher_or(this->pIdx, main_flag, fake, main_flag, *this->enc,
                           *this->eval, *this->runtime);

            // 5. update the stash.
            this->stash.push_back(this->packed_index[phyIndex_in_nextMap]);
            this->t += 1;

            // 6. get the l-th element from the stash element.
            std::vector<aby3::sbMatrix> stashed_elements(this->pack);
            for (size_t i = 0; i < this->pack; i++) {
                stashed_elements[i].resize(1, BITSIZE);
                stashed_elements[i].mShares[0](0, 0) =
                    this->packed_index[phyIndex_in_nextMap]
                        .packedIndices[i]
                        .indexShares[0];
                stashed_elements[i].mShares[1](0, 0) =
                    this->packed_index[phyIndex_in_nextMap]
                        .packedIndices[i]
                        .indexShares[1];
            }

            aby3::sbMatrix _res_stashed_element(1, BITSIZE);
            this->linear_ram(stashed_elements, l, _res_stashed_element);

            aby3::sbMatrix _res(1, BITSIZE);
            bool_cipher_selector(this->pIdx, main_flag, _res_stashed_element,
                                 _phyIndex_in_stash, _res, *this->enc,
                                 *this->eval, *this->runtime);
            physical_index = boolIndex(_res);
        }
        aby3::i64 plain_index =
            back2plain(this->pIdx, physical_index, *this->enc, *this->eval,
                       *this->runtime);
        return plain_index;
    }

    void linear_ram(std::vector<aby3::sbMatrix> &data, boolIndex &index,
                    aby3::sbMatrix &res) {
        size_t len = data.size();

        // construct the s1_map.
        aby3::sbMatrix s1_map(len, 1);
        aby3::sbMatrix _tmp_index(len, 1);
        aby3::i64Matrix _tmp_range(len, 1);

        for (size_t i = 0; i < len; i++) {
            _tmp_index.mShares[0](i, 0) = index.indexShares[0];
            _tmp_index.mShares[1](i, 0) = index.indexShares[1];
            _tmp_range(i, 0) = i;
        }
        _tmp_index.resize(len, BITSIZE);
        bool_cipher_eq(this->pIdx, _tmp_index, _tmp_range, s1_map, *this->enc,
                       *this->eval, *this->runtime);

        // using dot to get the result.
        bool_cipher_dot(this->pIdx, data, s1_map, res, *this->enc, *this->eval,
                        *this->runtime);

        return;
    }
};

class ABY3SqrtOram : public SqrtOram<aby3::sbMatrix> {
   public:
    ABY3PosMap *posMap;
    int pIdx;
    aby3::Sh3Encryptor *enc;
    aby3::Sh3Evaluator *eval;
    aby3::Sh3Runtime *runtime;
    std::vector<ABY3StashElement<aby3::sbMatrix>> stash;

    ABY3SqrtOram(int n, int S, int pack, int pIdx, aby3::Sh3Encryptor &enc,
                 aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime)
        : SqrtOram<aby3::sbMatrix>(n, S, pack),
          pIdx(pIdx),
          enc(&enc),
          runtime(&runtime),
          eval(&eval) {
            int sqrt_n = (int) ceil(sqrt(n));
            if(S > sqrt_n) S = sqrt_n;
          }

    void initiate(std::vector<aby3::sbMatrix> &data) {
        std::vector<aby3::si64> Pi(this->n);
        efficient_shuffle_with_random_permutation(
            data, this->pIdx, this->shuffle_mem, Pi, *(this->enc),
            *(this->eval), *(this->runtime));
        std::vector<boolIndex> permutation(this->n);
        for (int i = 0; i < this->n; i++) {
            permutation[i] = boolIndex(Pi[i].mData[0], Pi[i].mData[1]);
        }
        this->posMap = new ABY3PosMap(this->n, this->pack, this->S, permutation,
                                      this->pIdx, *(this->enc), *(this->eval),
                                      *(this->runtime));
    }

    aby3::sbMatrix access(boolIndex index) {
        // first look at the stash.
        boolShare found_flag(false, this->pIdx);

        aby3::sbMatrix _res;

        if (this->t > 0) {
            aby3::sbMatrix _stash_index(this->t, BITSIZE);
            aby3::sbMatrix expand_index(this->t, BITSIZE);
            std::vector<aby3::sbMatrix> expand_val(this->t);
            for (size_t i = 0; i < this->t; i++) {
                _stash_index.mShares[0](i, 0) =
                    this->stash[i].logicalIndex.indexShares[0];
                _stash_index.mShares[1](i, 0) =
                    this->stash[i].logicalIndex.indexShares[1];
                expand_index.mShares[0](i, 0) = index.indexShares[0];
                expand_index.mShares[1](i, 0) = index.indexShares[1];
                expand_val[i] = this->stash[i].data;
            }
            bool_cipher_eq(this->pIdx, _stash_index, expand_index, _stash_index,
                           *(this->enc), *(this->eval), *(this->runtime));
            for (size_t i = 0; i < this->t;
                 i++) {  // exist one 1 will turn found flag to 1.
                found_flag.bshares[0] ^= _stash_index.mShares[0](i, 0);
                found_flag.bshares[1] ^= _stash_index.mShares[1](i, 0);
            }

            // the dot should support any length sbMatrix.
            size_t block_len = this->shuffle_mem[0].rows();
            bool_cipher_dot(this->pIdx, expand_val, _stash_index, _res,
                            *(this->enc), *(this->eval), *(this->runtime));
        }

        // get the physical index from the posMap.
        aby3::i64 phy_index = this->posMap->access(index, found_flag);
        // get the real result.
        aby3::sbMatrix res = this->shuffle_mem[phy_index];

        if (this->t > 0) {
            bool_cipher_selector(this->pIdx, found_flag, _res, res, res,
                                 *(this->enc), *(this->eval), *(this->runtime));
        }
        this->stash.push_back(ABY3StashElement<aby3::sbMatrix>(res, index));

        return res;
    }
};