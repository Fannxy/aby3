#include "Shuffle.h"

using namespace oc;
using namespace aby3;

// #define DEBUG_SHUFFLE2

static size_t MAX_COMM_SIZE = 1 << 25;

/**
 * The shuffle protocol accoring to
 * https://dl.acm.org/doi/10.1145/3460120.3484560 (advanced version).
 */
int efficient_shuffle(std::vector<sbMatrix>& T, int pIdx,
                      std::vector<sbMatrix>& Tres, Sh3Encryptor& enc,
                      Sh3Evaluator& eval, Sh3Runtime& runtime) {
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = T.size();
    size_t unit_len = T[0].i64Size();

    // generate the prev, next - correlated randomness.
    // 1 - generate the permutations.
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    // 2 - generate the random masks Z.
    std::vector<i64Matrix> prev_maskZ(len);
    std::vector<i64Matrix> next_maskZ(len);
    for (size_t i = 0; i < len; i++) {
        prev_maskZ[i].resize(unit_len, T[0].bitCount());
        next_maskZ[i].resize(unit_len, T[0].bitCount());
        get_random_mask(pIdx, prev_maskZ[i], prevSeed);
        get_random_mask(pIdx, next_maskZ[i], nextSeed);
    }

    // get the random permutation.
    if (pIdx == 0) {
        // pre-generate the randomness.
        std::vector<i64Matrix> maskB(len);
        std::vector<i64Matrix> maskA(len);
        for (size_t i = 0; i < len; i++) {
            maskB[i].resize(unit_len, T[0].bitCount());
            maskA[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskB[i], nextSeed);
            get_random_mask(pIdx, maskA[i], prevSeed);
        }

        // get the sharedX1.
        std::vector<i64Matrix> sharedX1(len);
        for (size_t i = 0; i < len; i++) {
            sharedX1[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedX1[i](j, 0) =
                    T[i].mShares[0](j) ^ T[i].mShares[1](j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedX1);

        // get the sharedX2.
        std::vector<i64Matrix> sharedX2(len);
        for (size_t i = 0; i < len; i++) {
            sharedX2[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedX2[i](j, 0) = sharedX1[i](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, sharedX2);

        // send the sharedX2 to P1.
        i64Matrix sharedX2_matrix(len * unit_len, T[0].bitCount());
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                sharedX2_matrix(i * unit_len + j, 0) = sharedX2[i](j, 0);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(sharedX2_matrix.data(),
                                          sharedX2_matrix.size());

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskA[i](j);
                Tres[i].mShares[0](j, 0) = maskB[i](j);
            }
        }
    }
    if (pIdx == 1) {
        // pre-generate the randomness.
        std::vector<i64Matrix> maskB(len);
        for (size_t i = 0; i < len; i++) {
            maskB[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskB[i], prevSeed);
        }

        // compute the sharedY1 and send to the next party.
        std::vector<i64Matrix> shared_Y1(len);
        for (size_t i = 0; i < len; i++) {
            shared_Y1[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                shared_Y1[i](j, 0) = T[i].mShares[0](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, shared_Y1);

        i64Matrix sharedY1_matrix(len * unit_len, T[0].bitCount());
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                sharedY1_matrix(i * unit_len + j, 0) = shared_Y1[i](j, 0);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(sharedY1_matrix.data(),
                                          sharedY1_matrix.size());

        i64Matrix sharedX2_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(sharedX2_matrix.data(),
                                 sharedX2_matrix.size());

        // compute the sharedX3.
        std::vector<i64Matrix> sharedX3(len);
        for (size_t i = 0; i < len; i++) {
            sharedX3[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedX3[i](j, 0) =
                    sharedX2_matrix(i * unit_len + j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedX3);

        // compute the masked C1.
        // std::vector<i64Matrix> maskedC1(len);
        i64Matrix maskedC1_matrix(len * unit_len, T[0].bitCount());
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                maskedC1_matrix(i * unit_len + j, 0) =
                    sharedX3[i](j) ^ maskB[i](j);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(maskedC1_matrix.data(),
                                          maskedC1_matrix.size());

        i64Matrix maskedC2_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mNext.recv(maskedC2_matrix.data(),
                                 maskedC2_matrix.size());

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskB[i](j);
                Tres[i].mShares[0](j, 0) = maskedC1_matrix(i * unit_len + j) ^
                                           maskedC2_matrix(i * unit_len + j);
            }
        }
    }
    if (pIdx == 2) {
        // pre-generate the randomness.
        std::vector<i64Matrix> maskA(len);
        for (size_t i = 0; i < len; i++) {
            maskA[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskA[i], nextSeed);
        }

        // receive the sharedY1 from the previous party.
        i64Matrix sharedY1_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(sharedY1_matrix.data(),
                                 sharedY1_matrix.size());

        // compute the sharedY2.
        std::vector<i64Matrix> sharedY2(len);
        for (size_t i = 0; i < len; i++) {
            sharedY2[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedY2[i](j, 0) =
                    sharedY1_matrix(i * unit_len + j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedY2);

        std::vector<i64Matrix> sharedY3(len);
        for (size_t i = 0; i < len; i++) {
            sharedY3[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedY3[i](j, 0) = sharedY2[i](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, sharedY3);

        // compute the masked C2.
        i64Matrix maskedC2_matrix(len * unit_len, T[0].bitCount());
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                maskedC2_matrix(i * unit_len + j, 0) =
                    sharedY3[i](j) ^ maskA[i](j);
            }
        }
        runtime.mComm.mPrev.asyncSendCopy(maskedC2_matrix.data(),
                                          maskedC2_matrix.size());
        i64Matrix maskedC1_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(maskedC1_matrix.data(),
                                 maskedC1_matrix.size());

        std::vector<i64Matrix> maskC(len);
        for (size_t i = 0; i < len; i++) {
            maskC[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                maskC[i](j, 0) = maskedC1_matrix(i * unit_len + j) ^
                                 maskedC2_matrix(i * unit_len + j);
            }
        }

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskC[i](j);
                Tres[i].mShares[0](j, 0) = maskA[i](j);
            }
        }
    }
    return 0;
}


int efficient_shuffle(aby3::sbMatrix &T, int pIdx, aby3::sbMatrix &Tres, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = T.rows();
    size_t unit_size = (T.bitCount() + 63) / 64;

    // generate the prev, next - correlated randomness.
    // 1 - generate the permutations.
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    // 2 - generate the random masks Z.
    i64Matrix prev_maskZ(len, unit_size);
    i64Matrix next_maskZ(len, unit_size);
    get_random_mask(pIdx, prev_maskZ, prevSeed);
    get_random_mask(pIdx, next_maskZ, nextSeed);

    // get the random permutation.
    if(pIdx == 0){
        // pre-generate the randomness.
        i64Matrix maskB(len, unit_size);
        i64Matrix maskA(len, unit_size);
        get_random_mask(pIdx, maskB, nextSeed);
        get_random_mask(pIdx, maskA, prevSeed);

        // get the sharedX1.
        i64Matrix sharedX1(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX1(i, 0) = T.mShares[0](i, 0) ^ T.mShares[1](i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedX1);

        // get the sharedX2.
        i64Matrix sharedX2(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX2(i, 0) = sharedX1(i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedX2);

        // send the sharedX2 to P1.
        // runtime.mComm.mNext.asyncSendCopy(sharedX2.data(), sharedX2.size());
        large_data_sending(pIdx, sharedX2, runtime, true);

        // compute the final shares.
        Tres.resize(len, unit_size);
        for(size_t i=0; i<len; i++){
            Tres.mShares[1](i, 0) = maskA(i, 0);
            Tres.mShares[0](i, 0) = maskB(i, 0);
        }
    }
    if(pIdx == 1){
        // pre-generate the randomness.
        i64Matrix maskB(len, unit_size);
        get_random_mask(pIdx, maskB, prevSeed);

        // compute the sharedY1 and send to the next party.
        i64Matrix sharedY1(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY1(i, 0) = T.mShares[0](i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedY1);

        // runtime.mComm.mNext.asyncSendCopy(sharedY1.data(), sharedY1.size());
        large_data_sending(pIdx, sharedY1, runtime, true);

        i64Matrix sharedX2(len, unit_size);
        // runtime.mComm.mPrev.recv(sharedX2.data(), sharedX2.size());
        large_data_receiving(pIdx, sharedX2, runtime, true);

        // compute the sharedX3.
        i64Matrix sharedX3(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedX3(i, 0) = sharedX2(i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedX3);

        // // compute the masked C1.
        // i64Matrix maskedC1(len, unit_size);
        // for(size_t i=0; i<len; i++){
        //     maskedC1(i, 0) = sharedX3(i, 0) ^ maskB(i, 0);
        // }
        // // runtime.mComm.mNext.asyncSendCopy(maskedC1.data(), maskedC1.size());
        // large_data_sending(pIdx, maskedC1, runtime, true);

        i64Matrix maskedC2(len, unit_size);
        // runtime.mComm.mNext.recv(maskedC2.data(), maskedC2.size());
        large_data_receiving(pIdx, maskedC2, runtime, false);

        // compute the masked C1.
        i64Matrix maskedC1(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskedC1(i, 0) = sharedX3(i, 0) ^ maskB(i, 0);
        }
        // runtime.mComm.mNext.asyncSendCopy(maskedC1.data(), maskedC1.size());
        large_data_sending(pIdx, maskedC1, runtime, true);

        // compute the final shares.
        Tres.resize(len, unit_size);
        for(size_t i=0; i<len; i++){
            Tres.mShares[1](i, 0) = maskB(i, 0);
            Tres.mShares[0](i, 0) = maskedC1(i, 0) ^ maskedC2(i, 0);
        }
    }
    if(pIdx == 2){
        // pre-generate the randomness.
        i64Matrix maskA(len, unit_size);
        get_random_mask(pIdx, maskA, nextSeed);

        // receive the sharedY1 from the previous party.
        i64Matrix sharedY1(len, unit_size);
        // runtime.mComm.mPrev.recv(sharedY1.data(), sharedY1.size());
        large_data_receiving(pIdx, sharedY1, runtime, true);

        // compute the sharedY2.
        i64Matrix sharedY2(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY2(i, 0) = sharedY1(i, 0) ^ next_maskZ(i, 0);
        }
        plain_permutate(next_permutation, sharedY2);

        i64Matrix sharedY3(len, unit_size);
        for(size_t i=0; i<len; i++){
            sharedY3(i, 0) = sharedY2(i, 0) ^ prev_maskZ(i, 0);
        }
        plain_permutate(prev_permutation, sharedY3);

        // compute the masked C2.
        i64Matrix maskedC2(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskedC2(i, 0) = sharedY3(i, 0) ^ maskA(i, 0);
        }
        // runtime.mComm.mPrev.asyncSendCopy(maskedC2.data(), maskedC2.size());
        large_data_sending(pIdx, maskedC2, runtime, false);

        i64Matrix maskedC1(len, unit_size);
        // runtime.mComm.mPrev.recv(maskedC1.data(), maskedC1.size());
        large_data_receiving(pIdx, maskedC1, runtime, true);

        i64Matrix maskC(len, unit_size);
        for(size_t i=0; i<len; i++){
            maskC(i, 0) = maskedC1(i, 0) ^ maskedC2(i, 0);
        }

        // compute the final shares.
        Tres.resize(len, unit_size);
        for(size_t i=0; i<len; i++){
            Tres.mShares[1](i, 0) = maskC(i, 0);
            Tres.mShares[0](i, 0) = maskA(i, 0);
        }
    }

    return 0;
}

// /**
//  * The shuffle protocol accoring to
//  * https://dl.acm.org/doi/10.1145/3460120.3484560 (advanced version).
//  */
// int efficient_shuffle_with_random_permutation(
//     std::vector<sbMatrix>& T, int pIdx, std::vector<sbMatrix>& Tres,
//     std::vector<aby3::si64>& Pi, aby3::Sh3Encryptor& enc,
//     aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime) {
//     // get the common randomness.
//     block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
//     block nextSeed = enc.mShareGen.mNextCommon.getSeed();
//     size_t len = T.size();
//     size_t unit_len = T[0].i64Size();
//     size_t bit_count = 1;  // currently, we only support 1-64 bit element representation.


//     // generate the prev, next - correlated randomness.
//     // 1 - generate the permutations and the inverse permutations.
//     std::vector<size_t> prev_permutation;
//     std::vector<size_t> next_permutation;
//     get_permutation(len, prev_permutation, prevSeed);
//     get_permutation(len, next_permutation, nextSeed);
//     std::vector<size_t> prev_inverse_permutation;
//     std::vector<size_t> next_inverse_permutation;
//     get_inverse_permutation(prev_permutation, prev_inverse_permutation);
//     get_inverse_permutation(next_permutation, next_inverse_permutation);

//     // init Pi with the plain shares of R.
//     Pi.resize(len);
//     for (size_t i = 0; i < len; i++) {
//         if (pIdx == 0) {
//             Pi[i].mData[0] = ~0;
//             Pi[i].mData[0] = i;
//         }
//         if (pIdx == 1) {
//             Pi[i].mData[0] = ~0;
//             Pi[i].mData[1] = ~0;
//         }
//         if (pIdx == 2) {
//             Pi[i].mData[0] = i;
//             Pi[i].mData[1] = ~0;
//         }
//     }

//     // generate the random masks Z and RZ.
//     std::vector<i64Matrix> prev_maskZ(len);
//     std::vector<i64Matrix> next_maskZ(len);
//     i64Matrix prev_maskRZ(len, bit_count);
//     i64Matrix next_maskRZ(len, bit_count);

//     for (size_t i = 0; i < len; i++) {
//         prev_maskZ[i].resize(unit_len, bit_count);
//         next_maskZ[i].resize(unit_len, bit_count);
//         get_random_mask(pIdx, prev_maskZ[i], prevSeed);
//         get_random_mask(pIdx, next_maskZ[i], nextSeed);
//     }
//     get_random_mask(pIdx, prev_maskRZ, prevSeed);
//     get_random_mask(pIdx, next_maskRZ, nextSeed);

// #ifdef DEBUG_SHUFFLE2
//     std::string party_debug_file = "/root/aby3/party-" + std::to_string(pIdx) + ".txt";
//     std::ofstream partyfs(party_debug_file, std::ios_base::app);
//     partyfs << "before shuffle ptotocol" << std::endl;
// #endif

//     // get the random permutation.
//     if (pIdx == 0) {
//         // pre-generate the randomness, masks and the maskRs.
//         std::vector<i64Matrix> maskB(len);
//         std::vector<i64Matrix> maskA(len);
//         for (size_t i = 0; i < len; i++) {
//             maskB[i].resize(unit_len, bit_count);
//             maskA[i].resize(unit_len, bit_count);
//             get_random_mask(pIdx, maskB[i], nextSeed);
//             get_random_mask(pIdx, maskA[i], prevSeed);
//         }
//         i64Matrix maskRA(len, bit_count);
//         get_random_mask(pIdx, maskRA, prevSeed);

//         // get the sharedX1.
//         std::vector<i64Matrix> sharedX1(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedX1[i].resize(unit_len, bit_count);
//             for (size_t j = 0; j < unit_len; j++) {
//                 sharedX1[i](j, 0) =
//                     T[i].mShares[0](j) ^ T[i].mShares[1](j) ^ next_maskZ[i](j);
//             }
//             next_maskZ[i].resize(0, 0);
//         }
//         plain_permutate(next_permutation, sharedX1);

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after sharedX1" << std::endl;
// #endif

//         // get the sharedX2.
//         std::vector<i64Matrix> sharedX2(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedX2[i].resize(unit_len, bit_count);
//             for (size_t j = 0; j < unit_len; j++) {
//                 sharedX2[i](j, 0) = sharedX1[i](j) ^ prev_maskZ[i](j);
//             }
//             prev_maskZ[i].resize(0, 0);
//             sharedX1[i].resize(0, 0);
//         }
//         plain_permutate(prev_permutation, sharedX2);

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after sharedX2" << std::endl;
// #endif

//         // send the sharedX2 to P1.
//         i64Matrix sharedX2_matrix(len * unit_len, bit_count);
//         for (size_t i = 0; i < len; i++) {
//             for (size_t j = 0; j < unit_len; j++) {
//                 sharedX2_matrix(i * unit_len + j, 0) = sharedX2[i](j, 0);
//             }
//             sharedX2[i].resize(0, 0);
//         }

//         size_t round = (size_t)ceil(sharedX2_matrix.size() / (double)MAX_COMM_SIZE);
//         size_t last_len = sharedX2_matrix.size() - (round - 1) * MAX_COMM_SIZE;

//         // per-round async communications.
//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;

//             aby3::i64Matrix _data = sharedX2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1);
//             auto sendFu = runtime.mComm.mNext.asyncSendFuture(_data.data(), _data.size());
//             sendFu.get();
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after first communication" << std::endl;
// #endif

//         // delete the space.
//         sharedX2_matrix.resize(0, 0);

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after space deletion of sharedX1 and sharedX2" << std::endl;
// #endif

//         // the following is for the shared permutation.
//         // receive the sharedRY1 from P1.
//         i64Matrix sharedRY1_matrix(len, bit_count);
//         // runtime.mComm.mNext.recv(sharedRY1_matrix.data(),
//                                 //  sharedRY1_matrix.size());
//         round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
//         last_len = len - (round - 1) * MAX_COMM_SIZE;
//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             aby3::i64Matrix _recv_buffer(_len, 1);

//             auto recvFu = runtime.mComm.mNext.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
//             recvFu.get();
//             sharedRY1_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "permutation share 1" << std::endl;
// #endif

//         // compute the sharedRY2.
//         std::vector<i64> sharedRY2(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedRY2[i] = sharedRY1_matrix(i, 0) ^ prev_maskRZ(i, 0);
//         }
//         plain_permutate(prev_inverse_permutation, sharedRY2);

//         // compute the sharedRY3.
//         std::vector<i64> sharedRY3(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedRY3[i] = sharedRY2[i] ^ next_maskRZ(i, 0);
//         }
//         plain_permutate(next_inverse_permutation, sharedRY3);

//         i64Matrix maskRB2(len, bit_count);
//         for (size_t i = 0; i < len; i++) {
//             maskRB2(i, 0) = sharedRY3[i] ^ maskRA(i, 0);
//         }

//         // runtime.mComm.mNext.asyncSendCopy(maskRB2.data(), maskRB2.size());

//         i64Matrix maskRB1(len, bit_count);
//         // runtime.mComm.mNext.recv(maskRB1.data(), maskRB1.size());

//         round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
//         last_len = len - (round - 1) * MAX_COMM_SIZE;
//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             aby3::i64Matrix _recv_buffer(_len, 1);
//             aby3::i64Matrix _data = maskRB2.block(i * MAX_COMM_SIZE, 0, _len, 1);

//             auto sendFu = runtime.mComm.mNext.asyncSendFuture(_data.data(), _data.size());
//             auto recvFu = runtime.mComm.mNext.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
//             sendFu.get();
//             recvFu.get();
//             maskRB1.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "permutation RB2" << std::endl;
// #endif

//         // compute the final shares.
//         for (size_t i = 0; i < len; i++) {
//             Tres[i].resize(unit_len, BITSIZE);
//             for (size_t j = 0; j < unit_len; j++) {
//                 Tres[i].mShares[1](j, 0) = maskA[i](j);
//                 Tres[i].mShares[0](j, 0) = maskB[i](j);
//             }
//         }

//         // compute the final shares of Pi.
//         for (size_t i = 0; i < len; i++) {
//             Pi[i].mData[1] = maskRA(i, 0);
//             Pi[i].mData[0] = maskRB1(i, 0) ^ maskRB2(i, 0);
//         }
//     }
//     if (pIdx == 1) {
//         // pre-generate the randomness.
//         std::vector<i64Matrix> maskB(len);
//         for (size_t i = 0; i < len; i++) {
//             maskB[i].resize(unit_len, bit_count);
//             get_random_mask(pIdx, maskB[i], prevSeed);
//         }
//         i64Matrix maskRC(len, bit_count);
//         get_random_mask(pIdx, maskRC, nextSeed);

//         // compute the sharedY1 and send to the next party.
//         std::vector<i64Matrix> shared_Y1(len);
//         for (size_t i = 0; i < len; i++) {
//             shared_Y1[i].resize(unit_len, bit_count);
//             for (size_t j = 0; j < unit_len; j++) {
//                 shared_Y1[i](j, 0) = T[i].mShares[0](j) ^ prev_maskZ[i](j);
//             }
//             prev_maskZ[i].resize(0, 0);
//         }
//         plain_permutate(prev_permutation, shared_Y1);

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after sharedY1" << std::endl;
// #endif

//         i64Matrix sharedY1_matrix(len * unit_len, bit_count);
//         for (size_t i = 0; i < len; i++) {
//             for (size_t j = 0; j < unit_len; j++) {
//                 sharedY1_matrix(i * unit_len + j, 0) = shared_Y1[i](j, 0);
//             }
//             // then delete the space of shared_Y1.
//             shared_Y1[i].resize(0, 0);
//         }
//         i64Matrix sharedX2_matrix(len * unit_len, bit_count);

//         size_t round = (size_t)ceil(sharedY1_matrix.size() / (double)MAX_COMM_SIZE);
//         size_t last_len = sharedY1_matrix.size() - (round - 1) * MAX_COMM_SIZE;

//         // send to P2 and recv from P0 in async manner.
//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             aby3::i64Matrix _data = sharedY1_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1);
//             aby3::i64Matrix _recv_buffer(_len, 1);

//             auto sendFu = runtime.mComm.mNext.asyncSendFuture(_data.data(), _data.size());
//             auto recvFu = runtime.mComm.mPrev.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
//             sendFu.get();
//             recvFu.get();
//             sharedX2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after first round communication" << std::endl;
// #endif

//         // delete the space of sharedY1.
//         sharedY1_matrix.resize(0, 0);

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after sharedY1_martix deletion" << std::endl;
// #endif

//         // compute the sharedX3.
//         std::vector<i64Matrix> sharedX3(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedX3[i].resize(unit_len, bit_count);
//             for (size_t j = 0; j < unit_len; j++) {
//                 sharedX3[i](j, 0) =
//                     sharedX2_matrix(i * unit_len + j) ^ next_maskZ[i](j);
//             }
//             next_maskZ[i].resize(0, 0);
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after sharedX3 deletion" << std::endl;
// #endif
        
//         // delete the space of sharedX2.
//         sharedX2_matrix.resize(0, 0);

//         plain_permutate(next_permutation, sharedX3);

//         // compute the masked C1.
//         // std::vector<i64Matrix> maskedC1(len);
//         i64Matrix maskedC1_matrix(len * unit_len, bit_count);
//         for (size_t i = 0; i < len; i++) {
//             for (size_t j = 0; j < unit_len; j++) {
//                 maskedC1_matrix(i * unit_len + j, 0) =
//                     sharedX3[i](j) ^ maskB[i](j);
//             }
//             sharedX3[i].resize(0, 0);
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after maskedC1" << std::endl;
// #endif

//         i64Matrix maskedC2_matrix(len * unit_len, bit_count);

//         // send maskC1 to the next and reav maskC2 from the next.
//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             aby3::i64Matrix _data = maskedC1_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1);
//             aby3::i64Matrix _recv_buffer(_len, 1);

//             auto sendFu = runtime.mComm.mNext.asyncSendFuture(_data.data(), _data.size());
//             auto recvFu = runtime.mComm.mNext.asyncRecv(_recv_buffer.data(), _recv_buffer.size());

//             sendFu.get();
//             recvFu.get();
//             maskedC2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after second rounf comm" << std::endl;
// #endif

//         // the following is for shared permutation.
//         // compute the sharedRY1
//         std::vector<i64> sharedRY1(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedRY1[i] = Pi[i].mData[1] ^ next_maskRZ(i, 0);
//         }
//         plain_permutate(next_inverse_permutation, sharedRY1);
//         // runtime.mComm.mPrev.asyncSendCopy(sharedRY1.data(), sharedRY1.size());

//         // receive the sharedRX2 from the next party.
//         i64Matrix sharedRX2_matrix(len, bit_count);
//         // runtime.mComm.mNext.recv(sharedRX2_matrix.data(),
//                                 //  sharedRX2_matrix.size());
//         round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
//         last_len = len - (round - 1) * MAX_COMM_SIZE;

//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             aby3::i64Matrix _recv_buffer(_len, 1);
//             // aby3::i64Matrix _data = sharedRY1.block(i * MAX_COMM_SIZE, 0, _len, 1);
//             std::vector<i64> _data(sharedRY1.begin() + i*MAX_COMM_SIZE, sharedRY1.begin() + i*MAX_COMM_SIZE + _len);
//             auto sendFu = runtime.mComm.mPrev.asyncSendFuture(_data.data(), _data.size());
//             auto recvFu = runtime.mComm.mNext.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
//             sendFu.get();
//             recvFu.get();
//             sharedRX2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "permutation sharedRX2 " << round << " length = " << last_len << std::endl;
// #endif


//         // compute the sharedRX3.
//         std::vector<i64> sharedRX3(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedRX3[i] = sharedRX2_matrix(i, 0) ^ prev_maskRZ(i, 0);
//         }
//         plain_permutate(prev_inverse_permutation, sharedRX3);

//         // compute the maskedRB1.
//         i64Matrix maskRB1(len, bit_count);
//         for (size_t i = 0; i < len; i++) {
//             maskRB1(i, 0) = sharedRX3[i] ^ maskRC(i, 0);
//         }
//         // runtime.mComm.mPrev.asyncSendCopy(maskRB1.data(), maskRB1.size());


//         i64Matrix maskRB2(len, bit_count);
//         // runtime.mComm.mPrev.recv(maskRB2.data(), maskRB2.size());
//         round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
//         last_len = len - (round - 1) * MAX_COMM_SIZE;
//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             aby3::i64Matrix _recv_buffer(_len, 1);
//             aby3::i64Matrix _data = maskRB1.block(i * MAX_COMM_SIZE, 0, _len, 1);

//             auto sendFu = runtime.mComm.mPrev.asyncSendFuture(_data.data(), _data.size());
//             auto recvFu = runtime.mComm.mPrev.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
//             sendFu.get();
//             recvFu.get();
//             maskRB2.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "permutation maskRB2 " << round << " length = " << last_len << std::endl;
// #endif

//         // compute the final shares.
//         for (size_t i = 0; i < len; i++) {
//             Tres[i].resize(unit_len, BITSIZE);
//             for (size_t j = 0; j < unit_len; j++) {
//                 Tres[i].mShares[1](j, 0) = maskB[i](j);
//                 Tres[i].mShares[0](j, 0) = maskedC1_matrix(i * unit_len + j) ^
//                                            maskedC2_matrix(i * unit_len + j);
//             }
//         }

//         // compute the final shares of Pi.
//         for (size_t i = 0; i < len; i++) {
//             Pi[i].mData[0] = maskRC(i, 0);
//             Pi[i].mData[1] = maskRB1(i, 0) ^ maskRB2(i, 0);
//         }
//     }
//     if (pIdx == 2) {
//         // pre-generate the randomness.
//         std::vector<i64Matrix> maskA(len);
//         for (size_t i = 0; i < len; i++) {
//             maskA[i].resize(unit_len, bit_count);
//             get_random_mask(pIdx, maskA[i], nextSeed);
//         }

//         i64Matrix maskRC(len, bit_count);
//         i64Matrix maskRA(len, bit_count);
//         get_random_mask(pIdx, maskRC, prevSeed);
//         get_random_mask(pIdx, maskRA, nextSeed);

//         // receive the sharedY1 from the previous party.
//         i64Matrix sharedY1_matrix(len * unit_len, bit_count);

//         size_t round = (size_t)ceil(sharedY1_matrix.size() / (double)MAX_COMM_SIZE);
//         size_t last_len = sharedY1_matrix.size() - (round - 1) * MAX_COMM_SIZE;
//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             aby3::i64Matrix _recv_buffer(_len, 1);
//             auto recvFu = runtime.mComm.mPrev.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
//             recvFu.get();
//             sharedY1_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
//         }


// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after first round comm" << std::endl;
// #endif

//         // compute the sharedY2.
//         std::vector<i64Matrix> sharedY2(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedY2[i].resize(unit_len, bit_count);
//             for (size_t j = 0; j < unit_len; j++) {
//                 sharedY2[i](j, 0) =
//                     sharedY1_matrix(i * unit_len + j) ^ next_maskZ[i](j);
//             }
//             next_maskZ[i].resize(0, 0);
//         }
//         plain_permutate(next_permutation, sharedY2);

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after sharedY2" << std::endl;
// #endif

//         std::vector<i64Matrix> sharedY3(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedY3[i].resize(unit_len, bit_count);
//             for (size_t j = 0; j < unit_len; j++) {
//                 sharedY3[i](j, 0) = sharedY2[i](j) ^ prev_maskZ[i](j);
//             }
//             prev_maskZ[i].resize(0, 0);
//         }
//         plain_permutate(prev_permutation, sharedY3);

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after sharedY3" << std::endl;
// #endif

//         // compute the masked C2.
//         i64Matrix maskedC2_matrix(len * unit_len, bit_count);
//         for (size_t i = 0; i < len; i++) {
//             for (size_t j = 0; j < unit_len; j++) {
//                 maskedC2_matrix(i * unit_len + j, 0) =
//                     sharedY3[i](j) ^ maskA[i](j);
//             }
//             sharedY3[i].resize(0, 0);
//         }
//         i64Matrix maskedC1_matrix(len * unit_len, bit_count);

//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             aby3::i64Matrix _send_buffer = maskedC2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1);
//             aby3::i64Matrix _recv_buffer(_len, 1);
//             auto sendFu = runtime.mComm.mPrev.asyncSendFuture(_send_buffer.data(), _send_buffer.size());
//             auto recvFu = runtime.mComm.mPrev.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
//             sendFu.get();
//             recvFu.get();
//             maskedC1_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after second round comm" << std::endl;
// #endif

//         std::vector<i64Matrix> maskC(len);
//         for (size_t i = 0; i < len; i++) {
//             maskC[i].resize(unit_len, bit_count);
//             for (size_t j = 0; j < unit_len; j++) {
//                 maskC[i](j, 0) = maskedC1_matrix(i * unit_len + j) ^
//                                  maskedC2_matrix(i * unit_len + j);
//             }
//         }

//         // delete the space of maskedC1_matrix and maskedC2_matrix.
//         maskedC1_matrix.resize(0, 0);
//         maskedC2_matrix.resize(0, 0);

//         // the following is for the shared permutation.

//         // compute the sharedRX1.
//         std::vector<i64> sharedRX1(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedRX1[i] = Pi[i].mData[0] ^ Pi[i].mData[1] ^ prev_maskRZ(i, 0);
//         }
//         plain_permutate(prev_inverse_permutation, sharedRX1);

//         // compute the sharedRX2.
//         std::vector<i64> sharedRX2(len);
//         for (size_t i = 0; i < len; i++) {
//             sharedRX2[i] = sharedRX1[i] ^ next_maskRZ(i, 0);
//         }
//         plain_permutate(next_inverse_permutation, sharedRX2);

//         // runtime.mComm.mPrev.asyncSendCopy(sharedRX2.data(), sharedRX2.size());
//         round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
//         last_len = len - (round - 1) * MAX_COMM_SIZE;
//         for(size_t i=0; i<round; i++){
//             size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
//             std::vector<i64> _data(sharedRX2.begin() + i*MAX_COMM_SIZE, sharedRX2.begin() + i*MAX_COMM_SIZE + _len);
//             auto sendFu = runtime.mComm.mPrev.asyncSendFuture(_data.data(), _data.size());
//             sendFu.get();
//         }

// #ifdef DEBUG_SHUFFLE2
//     partyfs << "after sharedRX2" << std::endl;
// #endif

//         // compute the final shares.
//         for (size_t i = 0; i < len; i++) {
//             Tres[i].resize(unit_len, BITSIZE);
//             for (size_t j = 0; j < unit_len; j++) {
//                 Tres[i].mShares[1](j, 0) = maskC[i](j);
//                 Tres[i].mShares[0](j, 0) = maskA[i](j);
//             }
//         }

//         // compute the final shares of Pi.
//         for (size_t i = 0; i < len; i++) {
//             Pi[i].mData[0] = maskRA(i, 0);
//             Pi[i].mData[1] = maskRC(i, 0);
//         }
//     }

//     return 0;
// }

int efficient_shuffle_with_random_permutation(
    std::vector<sbMatrix>& T, int pIdx, std::vector<sbMatrix>& Tres,
    std::vector<aby3::si64>& Pi, aby3::Sh3Encryptor& enc,
    aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime) {
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = T.size();
    size_t unit_len = T[0].i64Size();
    size_t bit_count = 1;  // currently, we only support 1-64 bit element representation.

    // generate the prev, next - correlated randomness.
    // 1 - generate the permutations and the inverse permutations.
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);
    std::vector<size_t> prev_inverse_permutation;
    std::vector<size_t> next_inverse_permutation;
    get_inverse_permutation(prev_permutation, prev_inverse_permutation);
    get_inverse_permutation(next_permutation, next_inverse_permutation);

    // init Pi with the plain shares of R.
    Pi.resize(len);
    for (size_t i = 0; i < len; i++) {
        if (pIdx == 0) {
            Pi[i].mData[0] = ~0;
            Pi[i].mData[0] = i;
        }
        if (pIdx == 1) {
            Pi[i].mData[0] = ~0;
            Pi[i].mData[1] = ~0;
        }
        if (pIdx == 2) {
            Pi[i].mData[0] = i;
            Pi[i].mData[1] = ~0;
        }
    }

    std::ofstream partyfs("/root/aby3/party-" + std::to_string(pIdx) + ".txt", std::ios::app);

#ifdef DEBUG_SHUFFLE2
    std::string party_debug_file = "/root/aby3/party-" + std::to_string(pIdx) + ".txt";
    std::ofstream partyfs(party_debug_file, std::ios_base::app);
    partyfs << "before shuffle ptotocol" << std::endl;
#endif

    // get the random permutation.
    if (pIdx == 0) {

        // get the sharedX1.
        std::vector<i64Matrix> sharedX1(len);
        for (size_t i = 0; i < len; i++) {
            i64Matrix next_maskZ_unit(unit_len, bit_count);
            get_random_mask(pIdx, next_maskZ_unit, nextSeed);
            sharedX1[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                sharedX1[i](j, 0) =
                    T[i].mShares[0](j) ^ T[i].mShares[1](j) ^ next_maskZ_unit(j, 0);
            }
        }
        plain_permutate(next_permutation, sharedX1);

        // get the sharedX2.
        std::vector<i64Matrix> sharedX2(len);
        for (size_t i = 0; i < len; i++) {
            i64Matrix prev_maskZ_unit(unit_len, bit_count);
            get_random_mask(pIdx, prev_maskZ_unit, prevSeed);
            sharedX2[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                // sharedX2[i](j, 0) = sharedX1[i](j) ^ prev_maskZ[i](j);
                sharedX2[i](j, 0) = sharedX1[i](j, 0) ^ prev_maskZ_unit(j, 0);
            }
            // prev_maskZ[i].resize(0, 0);
            sharedX1[i].resize(0, 0);
        }
        plain_permutate(prev_permutation, sharedX2);

        size_t round = (size_t)ceil((len*unit_len) / (double)MAX_COMM_SIZE);
        size_t last_len = (len*unit_len) - (round - 1) * MAX_COMM_SIZE;

        // per-round async communications.
        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;

            // aby3::i64Matrix _data = sharedX2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1);
            aby3::i64Matrix _data(_len, 1);
            for(size_t j=0; j<_len; j++){
                _data(j, 0) = sharedX2[(i*MAX_COMM_SIZE + j)/unit_len]((i*MAX_COMM_SIZE + j)%unit_len, 0);
            }
            auto sendFu = runtime.mComm.mNext.asyncSendFuture(_data.data(), _data.size());
            sendFu.get();
        }

#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 1" << std::endl;
#endif

        // delete the space.
        for(size_t i=0; i<len; i++){
            sharedX2[i].resize(0, 0);
        }
        // sharedX2_matrix.resize(0, 0);

        // the following is for the shared permutation.
        // receive the sharedRY1 from P1.
        // i64Matrix sharedRY1_matrix(len, bit_count);
        std::vector<i64> sharedRY2(len);
        // runtime.mComm.mNext.recv(sharedRY1_matrix.data(),
                                //  sharedRY1_matrix.size());
        round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
        last_len = len - (round - 1) * MAX_COMM_SIZE;
        i64Matrix prev_maskRZ(len, bit_count);
        get_random_mask(pIdx, prev_maskRZ, prevSeed);
        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix _recv_buffer(_len, 1);

            auto recvFu = runtime.mComm.mNext.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
            recvFu.get();
            for(size_t j=0; j<_len; j++){
                sharedRY2[i * MAX_COMM_SIZE + j] = _recv_buffer(j, 0) ^ prev_maskRZ(i * MAX_COMM_SIZE + j, 0);
            }
        }
        prev_maskRZ.resize(0, 0);

#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 2" << std::endl;
#endif

        // compute the sharedRY2.
        // std::vector<i64> sharedRY2(len);
        // for (size_t i = 0; i < len; i++) {
        //     sharedRY2[i] = sharedRY1_matrix(i, 0) ^ prev_maskRZ(i, 0);
        // }
        plain_permutate(prev_inverse_permutation, sharedRY2);

        // compute the sharedRY3.
        std::vector<i64> sharedRY3(len);
        i64Matrix next_maskRZ(len, bit_count);
        get_random_mask(pIdx, next_maskRZ, nextSeed);
        for (size_t i = 0; i < len; i++) {
            sharedRY3[i] = sharedRY2[i] ^ next_maskRZ(i, 0);
        }
        sharedRY2.resize(0);
        next_maskRZ.resize(0, 0);
        plain_permutate(next_inverse_permutation, sharedRY3);
        
        // std::vector<i64Matrix> maskA(len);
        // for (size_t i = 0; i < len; i++) {
        //     maskB[i].resize(unit_len, bit_count);
        //     maskA[i].resize(unit_len, bit_count);
        //     get_random_mask(pIdx, maskB[i], nextSeed);
        //     get_random_mask(pIdx, maskA[i], prevSeed);
        // }

        i64Matrix maskRA(len, bit_count);
        get_random_mask(pIdx, maskRA, prevSeed);
        i64Matrix maskRB2(len, bit_count);
        for (size_t i = 0; i < len; i++) {
            maskRB2(i, 0) = sharedRY3[i] ^ maskRA(i, 0);
            Pi[i].mData[1] = maskRA(i, 0);
        }
        maskRA.resize(0, 0);
        // runtime.mComm.mNext.asyncSendCopy(maskRB2.data(), maskRB2.size());

        // i64Matrix maskRB1(len, bit_count);
        // runtime.mComm.mNext.recv(maskRB1.data(), maskRB1.size());

        round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
        last_len = len - (round - 1) * MAX_COMM_SIZE;
        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix _recv_buffer(_len, 1);
            aby3::i64Matrix _data = maskRB2.block(i * MAX_COMM_SIZE, 0, _len, 1);

            auto sendFu = runtime.mComm.mNext.asyncSendFuture(_data.data(), _data.size());
            auto recvFu = runtime.mComm.mNext.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
            sendFu.get();
            recvFu.get();
            // maskRB1.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;

            for(size_t j=0; j<_len; j++){
                Pi[i * MAX_COMM_SIZE + j].mData[0] = _recv_buffer(j, 0) ^ maskRB2(i * MAX_COMM_SIZE + j, 0);
            }
        }

#ifdef DEBUG_SHUFFLE2
    partyfs << "permutation 3" << std::endl;
#endif


        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, BITSIZE);
            i64Matrix maskA_unit(unit_len, bit_count);
            i64Matrix maskB_unit(unit_len, bit_count);
            get_random_mask(pIdx, maskA_unit, prevSeed);
            get_random_mask(pIdx, maskB_unit, nextSeed);
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskA_unit(j, 0);
                Tres[i].mShares[0](j, 0) = maskB_unit(j, 0);
            }
        }

#ifdef DEBUG_SHUFFLE2
    partyfs << "finished!!!" << std::endl;
#endif

    }
    if (pIdx == 1) {

        // compute the sharedY1 and send to the next party.
        std::vector<i64Matrix> shared_Y1(len);
        for (size_t i = 0; i < len; i++) {
            i64Matrix perv_maskZ_unit(unit_len, bit_count);
            get_random_mask(pIdx, perv_maskZ_unit, prevSeed);
            shared_Y1[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                // shared_Y1[i](j, 0) = T[i].mShares[0](j) ^ prev_maskZ[i](j);
                shared_Y1[i](j, 0) = T[i].mShares[0](j) ^ perv_maskZ_unit(j, 0);
            }
            // prev_maskZ[i].resize(0, 0);
        }
        plain_permutate(prev_permutation, shared_Y1);


        size_t round = (size_t)ceil(len * unit_len / (double)MAX_COMM_SIZE);
        size_t last_len = len * unit_len - (round - 1) * MAX_COMM_SIZE;

        std::vector<i64Matrix> sharedX3(len);
        std::vector<i64Matrix> next_maskZ(len);
        for(size_t i=0; i<len; i++){
            sharedX3[i].resize(unit_len, bit_count);
            next_maskZ[i].resize(unit_len, bit_count);
            get_random_mask(pIdx, next_maskZ[i], nextSeed);
        }

        // send to P2 and recv from P0 in async manner.
        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            // aby3::i64Matrix _data = sharedY1_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1);
            aby3::i64Matrix _data(_len, 1);
            for(size_t j=0; j<_len; j++){
                _data(j, 0) = shared_Y1[(i*MAX_COMM_SIZE + j)/unit_len]((i*MAX_COMM_SIZE + j) % unit_len, 0);
            }
            aby3::i64Matrix _recv_buffer(_len, 1);

            auto sendFu = runtime.mComm.mNext.asyncSendFuture(_data.data(), _data.size());
            auto recvFu = runtime.mComm.mPrev.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
            sendFu.get();
            recvFu.get();
            for(size_t j=0; j<_len; j++){
                sharedX3[(i*MAX_COMM_SIZE + j) / unit_len]((i*MAX_COMM_SIZE + j) % unit_len, 0) = _recv_buffer(j, 0) ^ next_maskZ[(i*MAX_COMM_SIZE + j)/unit_len]((i*MAX_COMM_SIZE + j) % unit_len, 0);
            }
        }

#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 1" << std::endl;
#endif

        for (size_t i = 0; i < len; i++) {
            next_maskZ[i].resize(0, 0);
            shared_Y1[i].resize(0, 0);
        }
        plain_permutate(next_permutation, sharedX3);

        i64Matrix maskedC1_matrix(len * unit_len, bit_count);
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, BITSIZE);
            i64Matrix maskB_unit(unit_len, bit_count);
            get_random_mask(pIdx, maskB_unit, prevSeed);
            for (size_t j = 0; j < unit_len; j++) {
                sharedX3[i](j, 0) = sharedX3[i](j, 0) ^ maskB_unit(j, 0);
                Tres[i].mShares[1](j, 0) = maskB_unit(j, 0);
            }
        }

        // i64Matrix maskedC2_matrix(len * unit_len, bit_count);

        // send maskC1 to the next and reav maskC2 from the next.
        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix _data(_len, 1);
            for(size_t j=0; j<_len; j++){
                _data(j, 0) = sharedX3[(i*MAX_COMM_SIZE + j) / unit_len]((i*MAX_COMM_SIZE + j) % unit_len, 0);
            }
            aby3::i64Matrix _recv_buffer(_len, 1);

            auto sendFu = runtime.mComm.mNext.asyncSendFuture(_data.data(), _data.size());
            auto recvFu = runtime.mComm.mNext.asyncRecv(_recv_buffer.data(), _recv_buffer.size());

            sendFu.get();
            recvFu.get();
            // maskedC2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
            for(size_t j=0; j<_len; j++){
                Tres[(i*MAX_COMM_SIZE + j) / unit_len].mShares[0]((i*MAX_COMM_SIZE + j) % unit_len, 0) = _recv_buffer(j, 0) ^ _data(j, 0);
            }
        }

#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 2" << std::endl;
#endif

        // the following is for shared permutation.
        // compute the sharedRY1
        std::vector<i64> sharedRY1(len);
        i64Matrix next_maskRZ(len, bit_count);
        get_random_mask(pIdx, next_maskRZ, nextSeed);
        for (size_t i = 0; i < len; i++) {
            sharedX3[i].resize(0, 0);
            sharedRY1[i] = Pi[i].mData[1] ^ next_maskRZ(i, 0);
        }
        plain_permutate(next_inverse_permutation, sharedRY1);
        next_maskRZ.resize(0, 0);
        // runtime.mComm.mPrev.asyncSendCopy(sharedRY1.data(), sharedRY1.size());

        // receive the sharedRX2 from the next party.
        i64Matrix sharedRX2_matrix(len, bit_count);
        // runtime.mComm.mNext.recv(sharedRX2_matrix.data(),
                                //  sharedRX2_matrix.size());
        round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
        last_len = len - (round - 1) * MAX_COMM_SIZE;

        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix _recv_buffer(_len, 1);
            // aby3::i64Matrix _data = sharedRY1.block(i * MAX_COMM_SIZE, 0, _len, 1);
            std::vector<i64> _data(sharedRY1.begin() + i*MAX_COMM_SIZE, sharedRY1.begin() + i*MAX_COMM_SIZE + _len);
            auto sendFu = runtime.mComm.mPrev.asyncSendFuture(_data.data(), _data.size());
            auto recvFu = runtime.mComm.mNext.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
            sendFu.get();
            recvFu.get();
            sharedRX2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
        }

#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 3" << std::endl;
#endif

        // compute the sharedRX3.
        std::vector<i64> sharedRX3(len);
        i64Matrix prev_maskRZ(len, bit_count);
        get_random_mask(pIdx, prev_maskRZ, prevSeed);
        for (size_t i = 0; i < len; i++) {
            sharedRX3[i] = sharedRX2_matrix(i, 0) ^ prev_maskRZ(i, 0);
        }
        prev_maskRZ.resize(0, 0);
        plain_permutate(prev_inverse_permutation, sharedRX3);

        // compute the maskedRB1.
        i64Matrix maskRC(len, bit_count);
        get_random_mask(pIdx, maskRC, nextSeed);
        i64Matrix maskRB1(len, bit_count);
        for (size_t i = 0; i < len; i++) {
            maskRB1(i, 0) = sharedRX3[i] ^ maskRC(i, 0);
        }

        i64Matrix maskRB2(len, bit_count);
        // runtime.mComm.mPrev.recv(maskRB2.data(), maskRB2.size());
        round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
        last_len = len - (round - 1) * MAX_COMM_SIZE;
        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix _recv_buffer(_len, 1);
            aby3::i64Matrix _data = maskRB1.block(i * MAX_COMM_SIZE, 0, _len, 1);

            auto sendFu = runtime.mComm.mPrev.asyncSendFuture(_data.data(), _data.size());
            auto recvFu = runtime.mComm.mPrev.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
            sendFu.get();
            recvFu.get();
            maskRB2.block(i * MAX_COMM_SIZE, 0, _len, 1) = _recv_buffer;
        }
#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 4" << std::endl;
#endif
        // compute the final shares of Pi.
        for (size_t i = 0; i < len; i++) {
            Pi[i].mData[0] = maskRC(i, 0);
            Pi[i].mData[1] = maskRB1(i, 0) ^ maskRB2(i, 0);
        }


    }
    if (pIdx == 2) {

        // receive the sharedY1 from the previous party.
        i64Matrix sharedY1_matrix(len * unit_len, bit_count);

        std::vector<i64Matrix> sharedY2(len);
        std::vector<i64Matrix> next_maskZ(len);
        for (size_t i = 0; i < len; i++) {
            next_maskZ[i].resize(unit_len, bit_count);
            get_random_mask(pIdx, next_maskZ[i], nextSeed);
            sharedY2[i].resize(unit_len, bit_count);
        }

        size_t round = (size_t)ceil(sharedY1_matrix.size() / (double)MAX_COMM_SIZE);
        size_t last_len = sharedY1_matrix.size() - (round - 1) * MAX_COMM_SIZE;
        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix _recv_buffer(_len, 1);
            auto recvFu = runtime.mComm.mPrev.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
            recvFu.get();
            for(size_t j=0; j<_len; j++){
                sharedY2[(i*MAX_COMM_SIZE+j)/unit_len]((i*MAX_COMM_SIZE+j) % unit_len, 0) = _recv_buffer(j, 0) ^ next_maskZ[(i*MAX_COMM_SIZE+j)/unit_len]((i*MAX_COMM_SIZE+j) % unit_len, 0);
            }
        }
        for(size_t i=0; i<len; i++) next_maskZ[i].resize(0, 0);
        // compute the sharedY2.
        plain_permutate(next_permutation, sharedY2);

#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 1" << std::endl;
#endif

        std::vector<i64Matrix> sharedY3(len);
        for (size_t i = 0; i < len; i++) {
            sharedY3[i].resize(unit_len, bit_count);
            i64Matrix prev_maskZ_unit(unit_len, bit_count);
            get_random_mask(pIdx, prev_maskZ_unit, prevSeed);
            for (size_t j = 0; j < unit_len; j++) {
                sharedY3[i](j, 0) = sharedY2[i](j) ^ prev_maskZ_unit(j, 0);
            }
            // prev_maskZ[i].resize(0, 0);
            sharedY2[i].resize(0, 0);
        }
        plain_permutate(prev_permutation, sharedY3);


        // compute the masked C2.
        i64Matrix maskedC2_matrix(len * unit_len, bit_count);
        for (size_t i = 0; i < len; i++) {
            i64Matrix maskA_unit(unit_len, bit_count);
            get_random_mask(pIdx, maskA_unit, nextSeed);
            Tres[i].resize(unit_len, BITSIZE);
            for (size_t j = 0; j < unit_len; j++) {
                maskedC2_matrix(i * unit_len + j, 0) =
                    sharedY3[i](j, 0) ^ maskA_unit(j, 0);
                Tres[i].mShares[0](j, 0) = maskA_unit(j, 0);
            }
            sharedY3[i].resize(0, 0);
        }
        i64Matrix maskedC1_matrix(len * unit_len, bit_count);

        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            aby3::i64Matrix _send_buffer = maskedC2_matrix.block(i * MAX_COMM_SIZE, 0, _len, 1);
            aby3::i64Matrix _recv_buffer(_len, 1);
            auto sendFu = runtime.mComm.mPrev.asyncSendFuture(_send_buffer.data(), _send_buffer.size());
            auto recvFu = runtime.mComm.mPrev.asyncRecv(_recv_buffer.data(), _recv_buffer.size());
            sendFu.get();
            recvFu.get();
            for(size_t j=0; j<_len; j++) {
                Tres[(i*MAX_COMM_SIZE+j)/unit_len].mShares[1]((i*MAX_COMM_SIZE+j)%unit_len, 0) = _recv_buffer(j, 0) ^ maskedC2_matrix(i*MAX_COMM_SIZE+j, 0);
            }
        }

#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 2" << std::endl;
#endif

        // maskedC2_matrix.resize(0, 0);

        // the following is for the shared permutation.

        // compute the sharedRX1.
        std::vector<i64> sharedRX1(len);
        i64Matrix prev_maskRZ(len, bit_count);
        get_random_mask(pIdx, prev_maskRZ, prevSeed);
        for (size_t i = 0; i < len; i++) {
            sharedRX1[i] = Pi[i].mData[0] ^ Pi[i].mData[1] ^ prev_maskRZ(i, 0);
        }
        prev_maskRZ.resize(0, 0);
        plain_permutate(prev_inverse_permutation, sharedRX1);

        // compute the sharedRX2.
        std::vector<i64> sharedRX2(len);
        i64Matrix next_maskRZ(len, bit_count);
        get_random_mask(pIdx, next_maskRZ, nextSeed);
        for (size_t i = 0; i < len; i++) {
            sharedRX2[i] = sharedRX1[i] ^ next_maskRZ(i, 0);
        }
        next_maskRZ.resize(0, 0);
        plain_permutate(next_inverse_permutation, sharedRX2);

        round = (size_t)ceil(len / (double)MAX_COMM_SIZE);
        last_len = len - (round - 1) * MAX_COMM_SIZE;
        for(size_t i=0; i<round; i++){
            size_t _len = (i == round - 1) ? last_len : MAX_COMM_SIZE;
            std::vector<i64> _data(sharedRX2.begin() + i*MAX_COMM_SIZE, sharedRX2.begin() + i*MAX_COMM_SIZE + _len);
            auto sendFu = runtime.mComm.mPrev.asyncSendFuture(_data.data(), _data.size());
            sendFu.get();
        }
        sharedRX2.resize(0);

#ifdef DEBUG_SHUFFLE2
        partyfs << "permutation share 3" << std::endl;
#endif

        i64Matrix maskRC(len, bit_count);
        i64Matrix maskRA(len, bit_count);
        get_random_mask(pIdx, maskRC, prevSeed);
        get_random_mask(pIdx, maskRA, nextSeed);

        // compute the final shares of Pi.
        for (size_t i = 0; i < len; i++) {
            Pi[i].mData[0] = maskRA(i, 0);
            Pi[i].mData[1] = maskRC(i, 0);
        }

    }

#ifdef DEBUG_SHUFFLE2
    partyfs << "before return" << std::endl;
#endif

    return 0;
}