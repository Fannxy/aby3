#include "Basics.h"

static const size_t MAX_SENDING_SIZE = 1 << 25;

int large_data_sending(int pIdx, aby3::i64Matrix &sharedA,
                       aby3::Sh3Runtime &runtime, bool toNext) {
    size_t len = sharedA.rows();
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    // std::ofstream party_fs(PARTY_FILE + std::to_string(pIdx) + ".txt", std::ios::app);

    if (toNext){
        for(size_t i=0; i<round; i++){
            size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
            aby3::i64Matrix sending_data = sharedA.block(i * MAX_SENDING_SIZE, 0, sending_len, 1);
            auto sendFu = runtime.mComm.mNext.asyncSendFuture(sending_data.data(), sending_data.size());
            sendFu.get();
        }
    }
    else{
        for(size_t i=0; i<round; i++){
            size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
            aby3::i64Matrix sending_data = sharedA.block(i * MAX_SENDING_SIZE, 0, sending_len, 1);
            auto sendFu = runtime.mComm.mPrev.asyncSendFuture(sending_data.data(), sending_data.size());
            sendFu.get();
            // debug_info("success send round-" + std::to_string(i), party_fs);
        }
    }

    // party_fs.close();

    return 0;
}

int large_data_receiving(int pIdx, aby3::i64Matrix &res,
                         aby3::Sh3Runtime &runtime, bool fromPrev) {
    size_t len = res.rows();
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    if(fromPrev){
        for(size_t i=0; i<round; i++){
            size_t recv_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
            aby3::i64Matrix receiving_data(recv_len, 1);
            auto recvFu = runtime.mComm.mPrev.asyncRecv(receiving_data.data(), receiving_data.size());
            recvFu.get();
            res.block(i * MAX_SENDING_SIZE, 0, recv_len, 1) = receiving_data;
        }
    }
    else{
        for(size_t i=0; i<round; i++){
            size_t recv_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
            aby3::i64Matrix receiving_data(recv_len, 1);
            auto recvFu = runtime.mComm.mNext.asyncRecv(receiving_data.data(), receiving_data.size());
            recvFu.get();
            res.block(i * MAX_SENDING_SIZE, 0, recv_len, 1) = receiving_data;
        }
    }

    return 0;
}

std::vector<size_t> argwhere(aby3::i64Matrix &input, int target) {
    size_t len = input.size();
    std::vector<size_t> res;

    for (size_t i = 0; i < len; i++) {
        if (input(i, 0) == target) {
            res.push_back(i);
        }
    }
    return res;
}

std::vector<size_t> argwhere(std::vector<std::vector<int>> &input, int target) {
    size_t len = input.size();
    size_t unit_len = input[0].size();
    std::vector<size_t> res(len);

    for (size_t i = 0; i < len; i++) {
        for (size_t j = 0; j < unit_len; j++) {
            if (input[i][j] == target) {
                res[i] = j;  // assume only one element each row is equal to the
                             // target.
                break;
            }
        }
    }
    return res;
}

std::vector<size_t> argwhere(std::vector<size_t> &input, int target) {
    size_t len = input.size();
    std::vector<size_t> res;

    for (size_t i = 0; i < len; i++) {
        if (input[i] == target) {
            res.push_back(i);
        }
    }
    return res;
}

int large_data_encryption(int pIdx, aby3::i64Matrix &plainA,
                          aby3::sbMatrix &sharedA, aby3::Sh3Encryptor &enc,
                          aby3::Sh3Runtime &runtime) {
    size_t len = plainA.rows();
    sharedA.resize(len, BITSIZE);
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    for(size_t i=0; i<round; i++){
        size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
        aby3::i64Matrix sending_data = plainA.block(i * MAX_SENDING_SIZE, 0, sending_len, 1);
        aby3::sbMatrix encrypted_data(sending_len, 1);
        if(pIdx == 0){
            enc.localBinMatrix(runtime, sending_data, encrypted_data).get();
        }
        else{
            enc.remoteBinMatrix(runtime, encrypted_data).get();
        }
        std::memcpy(sharedA.mShares[0].data() + i * MAX_SENDING_SIZE, encrypted_data.mShares[0].data(), sending_len * sizeof(encrypted_data.mShares[0](0, 0)));
        std::memcpy(sharedA.mShares[1].data() + i * MAX_SENDING_SIZE, encrypted_data.mShares[1].data(), sending_len * sizeof(encrypted_data.mShares[1](0, 0)));
    }

    return 0;
}

int large_data_decryption(int pIdx, aby3::sbMatrix &sharedA, aby3::i64Matrix &plainA,
                          aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime){
    size_t len = sharedA.rows();
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    for(size_t i=0; i<round; i++){
        size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
        aby3::sbMatrix encrypted_data(sending_len, 1);
        std::memcpy(encrypted_data.mShares[0].data(), sharedA.mShares[0].data() + i * MAX_SENDING_SIZE, sending_len * sizeof(sharedA.mShares[0](0, 0)));
        std::memcpy(encrypted_data.mShares[1].data(), sharedA.mShares[1].data() + i * MAX_SENDING_SIZE, sending_len * sizeof(sharedA.mShares[1](0, 0)));
        aby3::i64Matrix decrypted_data(sending_len, 1);
        enc.revealAll(runtime, encrypted_data, decrypted_data).get();
        plainA.block(i * MAX_SENDING_SIZE, 0, sending_len, 1) = decrypted_data;
    }

    return 0;
}