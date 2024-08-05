#include "Basics.h"

static const size_t MAX_SENDING_SIZE = 1<<27;

int large_data_sending(int pIdx, aby3::i64Matrix &sharedA, aby3::Sh3Runtime &runtime, bool toNext){
    size_t len = sharedA.rows();
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    if(pIdx == 0) debug_info("In sending - len = " + std::to_string(len) + ", round = " + std::to_string(round) + ", last_len = " + std::to_string(last_len));

    if(round == 1){
        if(toNext) runtime.mComm.mNext.asyncSend(sharedA.data(), sharedA.size());
        else runtime.mComm.mPrev.asyncSend(sharedA.data(), sharedA.size());
        return 0;
    }

    for(size_t i=0; i<round; i++){
        size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
        aby3::i64Matrix sending_data = sharedA.block(i * MAX_SENDING_SIZE, 0, sending_len, 1);
        if(toNext) runtime.mComm.mNext.asyncSend(sending_data.data(), sending_data.size());
        else runtime.mComm.mPrev.asyncSend(sending_data.data(), sending_data.size());
    }

    return 0;
}

int large_data_receiving(int pIdx, aby3::i64Matrix &res, aby3::Sh3Runtime &runtime, bool fromPrev){
    size_t len = res.rows();
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    if(pIdx == 0) debug_info("In receving - len = " + std::to_string(len) + ", round = " + std::to_string(round) + ", last_len = " + std::to_string(last_len));

    if(round == 1){
        if(fromPrev) runtime.mComm.mPrev.recv(res.data(), res.size());
        else runtime.mComm.mNext.recv(res.data(), res.size());
        return 0;
    }

    for(size_t i=0; i<round; i++){
        size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
        aby3::i64Matrix receiving_data(sending_len, 1);
        if(fromPrev) runtime.mComm.mPrev.recv(receiving_data.data(), receiving_data.size());
        else runtime.mComm.mNext.recv(receiving_data.data(), receiving_data.size());
        res.block(i * MAX_SENDING_SIZE, 0, sending_len, 1) = receiving_data;
    }

    return 0;
}

std::vector<size_t> argwhere(aby3::i64Matrix& input, int target){
    size_t len = input.size();
    std::vector<size_t> res;

    for(size_t i=0; i<len; i++){
        if(input(i, 0) == target){
            res.push_back(i);
        }
    }
    return res;
}

std::vector<size_t> argwhere(std::vector<std::vector<int>>& input, int target){
    size_t len = input.size();
    size_t unit_len = input[0].size();
    std::vector<size_t> res(len);

    for(size_t i=0; i<len; i++){
        for(size_t j=0; j<unit_len; j++){
            if(input[i][j] == target){
                res[i] = j; // assume only one element each row is equal to the target.
                break;
            }
        }
    }
    return res;
}

std::vector<size_t> argwhere(std::vector<size_t>& input, int target){
    size_t len = input.size();
    std::vector<size_t> res;

    for(size_t i=0; i<len; i++){
        if(input[i] == target){
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

int large_data_encryption(int pIdx, aby3::i64Matrix &plainA, aby3::si64Matrix &sharedA,
                          aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime) {
    
    size_t len = plainA.rows();
    sharedA.resize(len, BITSIZE);
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    for(size_t i=0; i<round; i++){
        size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
        aby3::i64Matrix sending_data = plainA.block(i * MAX_SENDING_SIZE, 0, sending_len, 1);
        aby3::si64Matrix encrypted_data(sending_len, 1);
        if(pIdx == 0){
            enc.localIntMatrix(runtime, sending_data, encrypted_data).get();
        }
        else{
            enc.remoteIntMatrix(runtime, encrypted_data).get();
        }
        std::memcpy(sharedA.mShares[0].data() + i * MAX_SENDING_SIZE, encrypted_data.mShares[0].data(), sending_len * sizeof(encrypted_data.mShares[0](0, 0)));
        std::memcpy(sharedA.mShares[1].data() + i * MAX_SENDING_SIZE, encrypted_data.mShares[1].data(), sending_len * sizeof(encrypted_data.mShares[1](0, 0)));
    }

    return 0;

}