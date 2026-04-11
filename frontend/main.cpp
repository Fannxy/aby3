#include <map>
#include <chrono>
#include <random>
#include <thread>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <tests_cryptoTools/UnitTests.h>

#include "aby3_tests/Test.h"
#include "aby3_tests/aby3_tests.h"
#include "eric.h"

#include "../aby3-Feddb-Core/feddb.h"
#include "../aby3-Feddb-Core/utils.h"
#include "../aby3-Feddb-Core/genperm.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/debug.h"

using namespace oc;
using namespace aby3;

int main(int argc, char** argv) {
    oc::CLP cmd(argc, argv);


    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    std::string timefile = "";
    if (cmd.isSet("timefile")) {
        auto tf = cmd.getMany<std::string>("timefile");
        timefile = tf[0];
    }

    std::cout << "start aby3: " << std::to_string(role) << std::endl;

    auto t_total_start = std::chrono::high_resolution_clock::now();

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

//create
    std::vector<i64Matrix> lineitem_filtered_2(3);
    i64Matrix lineitem_filtered_2_unit_size(1, 1);
    if (role == (2 - 1)%3 ){
        std::string path = "/root/plantest/benchmarks/tpch_test/32768/db/2/data/update/lineitem_filtered_2.csv" ;
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }

        std::string lineitem_filtered_2_line;
        std::vector<std::vector<i64>> lineitem_filtered_2_data_rows;
        
        //read data from path
        while (std::getline(file, lineitem_filtered_2_line)) {
            if (lineitem_filtered_2_line.empty()) continue;
            
            std::vector<i64> row;
            std::istringstream iss(lineitem_filtered_2_line);
            std::string token;
            
            while (std::getline(iss, token, ',')) {
                if (!token.empty()) {
                    row.push_back(std::stoll(token));
                }
            }
            
            if (!row.empty()) {
                lineitem_filtered_2_data_rows.push_back(row);
            }
        }
        file.close();
        
        if (lineitem_filtered_2_data_rows.empty()) {
            throw std::runtime_error("No data found in file: " + path);
        }
        
        size_t lineitem_filtered_2_num_cols = lineitem_filtered_2_data_rows[0].size();
        size_t lineitem_filtered_2_num_rows = lineitem_filtered_2_data_rows.size();

        
        lineitem_filtered_2_unit_size(0, 0) = lineitem_filtered_2_num_rows ;
        large_data_sending(role, lineitem_filtered_2_unit_size, runtime, true);
        large_data_sending(role, lineitem_filtered_2_unit_size, runtime, false);
        
        for (size_t col = 0; col < lineitem_filtered_2_num_cols; col++) {
            lineitem_filtered_2[col].resize(lineitem_filtered_2_num_rows, 1); 
            for (size_t row = 0; row < lineitem_filtered_2_num_rows; row++) {
                lineitem_filtered_2[col](row, 0) = lineitem_filtered_2_data_rows[row][col];
            }
        }
        
    }
    else if(role == (2)%3){
        large_data_receiving(role, lineitem_filtered_2_unit_size, runtime, true);
    }
    else if(role == (2 + 1)%3){
        large_data_receiving(role, lineitem_filtered_2_unit_size, runtime, false);
    }
    
//close
    std::vector<i64Matrix> lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols;
    std::vector<si64Matrix> lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols;

    std::vector<int> lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_idx = {};
    std::vector<int> lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_idx = {0, 1, 2};

    // 合并 plain 列通信：先 pack 成一个大矩阵，一次通信，再 unpack 回各列
    if (!lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_idx.empty()) {
        size_t lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_rows = lineitem_filtered_2_unit_size(0, 0);
        size_t lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_count = lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_idx.size();
        i64Matrix lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_count *  lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_rows, 1);

        if (role == (2 - 1) % 3) {
            for (u64 i = 0; i < lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_count; ++i) {
                const auto &idx = lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_idx[i];
                std::memcpy(
                    lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged.data() + i * lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_rows,
                    lineitem_filtered_2[idx].data(),
                    lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_rows * sizeof(lineitem_filtered_2[idx](0, 0)));
            }
        }

        if (role == (2 - 1) % 3) {
            large_data_sending(role, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged, runtime, true);
            large_data_sending(role, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged, runtime, false);
        } else if (role == (2) % 3) {
            lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged.resize(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_count * lineitem_filtered_2_unit_size(0, 0), 1);
            large_data_receiving(role, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged, runtime, true);
        } else if (role == (2 + 1) % 3) {
            lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged.resize(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_count * lineitem_filtered_2_unit_size(0, 0), 1);
            large_data_receiving(role, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged, runtime, false);
        }

        for (u64 i = 0; i < lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_count; ++i) {
            i64Matrix lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_col(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_rows, 1);
            std::memcpy(
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_col.data(),
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged.data() + i * lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_rows,
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_rows * sizeof(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_merged(0, 0)));
            lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols.push_back(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_col);
        }
    }
    //只有有数据的那个role会有lineitem_filtered_2[idx]的数据，其他的只会有lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols

    /*
    //合在一起
    for(const auto &idx : lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_idx){
        if(role == (2 - 1)%3){
            large_data_sending(role, lineitem_filtered_2[idx], runtime, true);
            large_data_sending(role, lineitem_filtered_2[idx], runtime, false);
        }
        else if(role == (2)%3){
            lineitem_filtered_2[idx].resize(lineitem_filtered_2_unit_size(0, 0), 1);
            large_data_receiving(role, lineitem_filtered_2[idx], runtime, true);
        }
        else if(role == (2 + 1)%3){
            lineitem_filtered_2[idx].resize(lineitem_filtered_2_unit_size(0, 0), 1);
            large_data_receiving(role, lineitem_filtered_2[idx], runtime, false);
        }
        lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols.push_back(lineitem_filtered_2[idx]);
    }
    */


    // 合并 cipher 列加密：先 pack 成一个大矩阵，一次 local/remote，再 unpack 回各列
    if (!lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_idx.empty()) {
        size_t lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows = lineitem_filtered_2_unit_size(0, 0);
        size_t lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_count = lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_idx.size();

        i64Matrix lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_plain_merged(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_count * lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows, 1);
        si64Matrix lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_merged(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_count * lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows, 1);

        if (role == (2 - 1) % 3) {
            for (u64 i = 0; i < lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_count; ++i) {
                const auto &idx = lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_idx[i];
                std::memcpy(
                    lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_plain_merged.data() + i * lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows,
                    lineitem_filtered_2[idx].data(),
                    lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows * sizeof(lineitem_filtered_2[idx](0, 0)));
            }
            enc.localIntMatrix(runtime, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_plain_merged, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_merged).get();
        } else {
            enc.remoteIntMatrix(runtime, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_merged).get();
        }

        for (u64 i = 0; i < lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_count; ++i) {
            si64Matrix lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_i(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows, 1);
            std::memcpy(
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_i.mShares[0].data(),
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_merged.mShares[0].data() + i * lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows,
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows * sizeof(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_merged.mShares[0](0, 0)));
            std::memcpy(
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_i.mShares[1].data(),
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_merged.mShares[1].data() + i * lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows,
                lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_rows * sizeof(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_merged.mShares[1](0, 0)));
            lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols.push_back(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_i);
        }
    }

    /*
    for (const auto &idx : lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_idx){
        si64Matrix lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_i(lineitem_filtered_2_unit_size(0, 0), 1);
        if (role == (2 - 1)%3){
            enc.localIntMatrix(runtime, lineitem_filtered_2[idx], lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_i).get();
        } else {
            enc.remoteIntMatrix(runtime, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_i).get();
        }
        lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols.push_back(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_i);
    }
    */

//create
    std::vector<i64Matrix> orders_filtered_filtered_6(3);
    i64Matrix orders_filtered_filtered_6_unit_size(1, 1);
    if (role == (1 - 1)%3 ){
        std::string path = "/root/plantest/benchmarks/tpch_test/32768/db/1/data/update/orders_filtered_filtered_6.csv" ;
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }

        std::string orders_filtered_filtered_6_line;
        std::vector<std::vector<i64>> orders_filtered_filtered_6_data_rows;
        
        //read data from path
        while (std::getline(file, orders_filtered_filtered_6_line)) {
            if (orders_filtered_filtered_6_line.empty()) continue;
            
            std::vector<i64> row;
            std::istringstream iss(orders_filtered_filtered_6_line);
            std::string token;
            
            while (std::getline(iss, token, ',')) {
                if (!token.empty()) {
                    row.push_back(std::stoll(token));
                }
            }
            
            if (!row.empty()) {
                orders_filtered_filtered_6_data_rows.push_back(row);
            }
        }
        file.close();
        
        if (orders_filtered_filtered_6_data_rows.empty()) {
            throw std::runtime_error("No data found in file: " + path);
        }
        
        size_t orders_filtered_filtered_6_num_cols = orders_filtered_filtered_6_data_rows[0].size();
        size_t orders_filtered_filtered_6_num_rows = orders_filtered_filtered_6_data_rows.size();

        
        orders_filtered_filtered_6_unit_size(0, 0) = orders_filtered_filtered_6_num_rows ;
        large_data_sending(role, orders_filtered_filtered_6_unit_size, runtime, true);
        large_data_sending(role, orders_filtered_filtered_6_unit_size, runtime, false);
        
        for (size_t col = 0; col < orders_filtered_filtered_6_num_cols; col++) {
            orders_filtered_filtered_6[col].resize(orders_filtered_filtered_6_num_rows, 1); 
            for (size_t row = 0; row < orders_filtered_filtered_6_num_rows; row++) {
                orders_filtered_filtered_6[col](row, 0) = orders_filtered_filtered_6_data_rows[row][col];
            }
        }
        
    }
    else if(role == (1)%3){
        large_data_receiving(role, orders_filtered_filtered_6_unit_size, runtime, true);
    }
    else if(role == (1 + 1)%3){
        large_data_receiving(role, orders_filtered_filtered_6_unit_size, runtime, false);
    }
    
//close
    std::vector<i64Matrix> orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols;
    std::vector<si64Matrix> orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols;

    std::vector<int> orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_idx = {};
    std::vector<int> orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_idx = {0, 1, 2};

    // 合并 plain 列通信：先 pack 成一个大矩阵，一次通信，再 unpack 回各列
    if (!orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_idx.empty()) {
        size_t orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_rows = orders_filtered_filtered_6_unit_size(0, 0);
        size_t orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_count = orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_idx.size();
        i64Matrix orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_count *  orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_rows, 1);

        if (role == (1 - 1) % 3) {
            for (u64 i = 0; i < orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_count; ++i) {
                const auto &idx = orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_idx[i];
                std::memcpy(
                    orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged.data() + i * orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_rows,
                    orders_filtered_filtered_6[idx].data(),
                    orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_rows * sizeof(orders_filtered_filtered_6[idx](0, 0)));
            }
        }

        if (role == (1 - 1) % 3) {
            large_data_sending(role, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged, runtime, true);
            large_data_sending(role, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged, runtime, false);
        } else if (role == (1) % 3) {
            orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged.resize(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_count * orders_filtered_filtered_6_unit_size(0, 0), 1);
            large_data_receiving(role, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged, runtime, true);
        } else if (role == (1 + 1) % 3) {
            orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged.resize(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_count * orders_filtered_filtered_6_unit_size(0, 0), 1);
            large_data_receiving(role, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged, runtime, false);
        }

        for (u64 i = 0; i < orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_count; ++i) {
            i64Matrix orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_col(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_rows, 1);
            std::memcpy(
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_col.data(),
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged.data() + i * orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_rows,
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_rows * sizeof(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_merged(0, 0)));
            orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols.push_back(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_col);
        }
    }
    //只有有数据的那个role会有orders_filtered_filtered_6[idx]的数据，其他的只会有orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols

    /*
    //合在一起
    for(const auto &idx : orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_idx){
        if(role == (1 - 1)%3){
            large_data_sending(role, orders_filtered_filtered_6[idx], runtime, true);
            large_data_sending(role, orders_filtered_filtered_6[idx], runtime, false);
        }
        else if(role == (1)%3){
            orders_filtered_filtered_6[idx].resize(orders_filtered_filtered_6_unit_size(0, 0), 1);
            large_data_receiving(role, orders_filtered_filtered_6[idx], runtime, true);
        }
        else if(role == (1 + 1)%3){
            orders_filtered_filtered_6[idx].resize(orders_filtered_filtered_6_unit_size(0, 0), 1);
            large_data_receiving(role, orders_filtered_filtered_6[idx], runtime, false);
        }
        orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols.push_back(orders_filtered_filtered_6[idx]);
    }
    */


    // 合并 cipher 列加密：先 pack 成一个大矩阵，一次 local/remote，再 unpack 回各列
    if (!orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_idx.empty()) {
        size_t orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows = orders_filtered_filtered_6_unit_size(0, 0);
        size_t orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_count = orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_idx.size();

        i64Matrix orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_plain_merged(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_count * orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows, 1);
        si64Matrix orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_merged(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_count * orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows, 1);

        if (role == (1 - 1) % 3) {
            for (u64 i = 0; i < orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_count; ++i) {
                const auto &idx = orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_idx[i];
                std::memcpy(
                    orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_plain_merged.data() + i * orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows,
                    orders_filtered_filtered_6[idx].data(),
                    orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows * sizeof(orders_filtered_filtered_6[idx](0, 0)));
            }
            enc.localIntMatrix(runtime, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_plain_merged, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_merged).get();
        } else {
            enc.remoteIntMatrix(runtime, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_merged).get();
        }

        for (u64 i = 0; i < orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_count; ++i) {
            si64Matrix orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_i(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows, 1);
            std::memcpy(
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_i.mShares[0].data(),
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_merged.mShares[0].data() + i * orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows,
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows * sizeof(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_merged.mShares[0](0, 0)));
            std::memcpy(
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_i.mShares[1].data(),
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_merged.mShares[1].data() + i * orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows,
                orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_rows * sizeof(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_merged.mShares[1](0, 0)));
            orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols.push_back(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_i);
        }
    }

    /*
    for (const auto &idx : orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_idx){
        si64Matrix orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_i(orders_filtered_filtered_6_unit_size(0, 0), 1);
        if (role == (1 - 1)%3){
            enc.localIntMatrix(runtime, orders_filtered_filtered_6[idx], orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_i).get();
        } else {
            enc.remoteIntMatrix(runtime, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_i).get();
        }
        orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols.push_back(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_i);
    }
    */

//semi-join

   //left/right_rel to all_cipher
   std::vector<si64Matrix> orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7;
   both2cipher(role, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols_idx, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols_idx, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_plain_cols, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_cipher_cols
               , orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7, enc, eval, runtime);
   std::vector<si64Matrix> lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3;
   both2cipher(role, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols_idx, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols_idx, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_plain_cols, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_cipher_cols
               , lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3, enc, eval, runtime);


   //other_idxs
   std::vector<int> orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_join_idxs = {0} ;
   std::vector<int> lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_join_idxs = {0} ;
   std::vector<int> orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_other_idxs, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_other_idxs;

   if(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_join_idxs.size() != lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_join_idxs.size()){
      throw std::runtime_error("join key_cols is not consistent") ;
   }

   for(size_t i = 0; i < orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7.size(); i++){
      if (std::find(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_join_idxs.begin(),
                    orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_join_idxs.end(),
                    static_cast<int>(i)) == orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_join_idxs.end()) {
         orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_other_idxs.push_back(static_cast<int>(i));
      }
   }

   for(size_t i = 0; i < lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3.size(); i++){
      if (std::find(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_join_idxs.begin(),
                    lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_join_idxs.end(),
                    static_cast<int>(i)) == lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_join_idxs.end()) {
         lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_other_idxs.push_back(static_cast<int>(i));
      }
   }

   //分别得到两方的key和other列
   std::vector<si64Matrix> orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_key, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_other;
   std::vector<si64Matrix> lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_key, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_other;
   for(size_t i = 0; i < orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_join_idxs.size(); i++){
      orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_key.push_back(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7[orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_join_idxs[i]]);
   }
   for(size_t i = 0; i < orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_other_idxs.size(); i++){
      orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_other.push_back(orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7[orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_other_idxs[i]]);
   }
   for(size_t i = 0; i < lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_join_idxs.size(); i++){
      lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_key.push_back(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3[lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_join_idxs[i]]);
   }
   for(size_t i = 0; i < lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_other_idxs.size(); i++){
      lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_other.push_back(lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3[lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_other_idxs[i]]);
   }

   //semi-join
   std::vector<i64Matrix> orders_filtered_filtered_semi_join_8_plain_cols = {};
   std::vector<si64Matrix> orders_filtered_filtered_semi_join_8_cipher_cols ;

   semi_join(role, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_key, orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_other,
      lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_key, lineitem_filtered_close_to_orders_filtered_filtered_semi_join_3_other,
      orders_filtered_filtered_semi_join_8_cipher_cols, enc, eval, runtime);


   int orders_filtered_filtered_semi_join_8_size = orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_join_idxs.size() + orders_filtered_filtered_close_to_orders_filtered_filtered_semi_join_7_other_idxs.size() ;
   std::vector<int> orders_filtered_filtered_semi_join_8_plain_cols_idx = {} ;
   std::vector<int> orders_filtered_filtered_semi_join_8_cipher_cols_idx(orders_filtered_filtered_semi_join_8_size);
   for(size_t i = 0; i < orders_filtered_filtered_semi_join_8_size; i++){
      orders_filtered_filtered_semi_join_8_cipher_cols_idx[i] = i ;
   }

//project

    std::map<int, int> orders_filtered_filtered_semi_join_projected_9_select_map = {{1, 0}} ;

    std::vector<int> orders_filtered_filtered_semi_join_8_select_idx;
    for(const auto &pair : orders_filtered_filtered_semi_join_projected_9_select_map){
        orders_filtered_filtered_semi_join_8_select_idx.push_back(pair.first);
    }
    std::vector<int> orders_filtered_filtered_semi_join_8_plain_select_idx ;
    std::vector<int> orders_filtered_filtered_semi_join_8_cipher_select_idx ;

    std::vector<si64Matrix> orders_filtered_filtered_semi_join_projected_9_cipher_cols;
    std::vector<i64Matrix> orders_filtered_filtered_semi_join_projected_9_plain_cols ;

    //select_idx : plain / cipher
    for(size_t i = 0 ; i < orders_filtered_filtered_semi_join_8_plain_cols_idx.size() ; i++){
        if(std::find(orders_filtered_filtered_semi_join_8_select_idx.begin(), orders_filtered_filtered_semi_join_8_select_idx.end(), orders_filtered_filtered_semi_join_8_plain_cols_idx[i]) != orders_filtered_filtered_semi_join_8_select_idx.end()){
            orders_filtered_filtered_semi_join_8_plain_select_idx.push_back(orders_filtered_filtered_semi_join_8_plain_cols_idx[i]);
            orders_filtered_filtered_semi_join_projected_9_plain_cols.push_back(orders_filtered_filtered_semi_join_8_plain_cols[i]);
        }
    }

    for(size_t i = 0 ; i < orders_filtered_filtered_semi_join_8_cipher_cols_idx.size() ; i++){
        if(std::find(orders_filtered_filtered_semi_join_8_select_idx.begin(), orders_filtered_filtered_semi_join_8_select_idx.end(), orders_filtered_filtered_semi_join_8_cipher_cols_idx[i]) != orders_filtered_filtered_semi_join_8_select_idx.end()){
            orders_filtered_filtered_semi_join_8_cipher_select_idx.push_back(orders_filtered_filtered_semi_join_8_cipher_cols_idx[i]);
            orders_filtered_filtered_semi_join_projected_9_cipher_cols.push_back(orders_filtered_filtered_semi_join_8_cipher_cols[i]);
        }
    }
    

    std::vector<int> orders_filtered_filtered_semi_join_projected_9_plain_cols_idx;
    std::vector<int> orders_filtered_filtered_semi_join_projected_9_cipher_cols_idx;
    for(const auto &idx : orders_filtered_filtered_semi_join_8_plain_select_idx){
        orders_filtered_filtered_semi_join_projected_9_plain_cols_idx.push_back(orders_filtered_filtered_semi_join_projected_9_select_map[idx]);
    }

    for(const auto &idx : orders_filtered_filtered_semi_join_8_cipher_select_idx){
        orders_filtered_filtered_semi_join_projected_9_cipher_cols_idx.push_back(orders_filtered_filtered_semi_join_projected_9_select_map[idx]);
    }
    



//open
    if(role == (1) % 3){
        size_t orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count = orders_filtered_filtered_semi_join_projected_9_cipher_cols.size();
        if (orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count > 0) {
            size_t orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows = orders_filtered_filtered_semi_join_projected_9_cipher_cols[0].rows();
            size_t orders_filtered_filtered_semi_join_projected_open_10_cipher_merge_rows = orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count * orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows;
            si64Matrix orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge(orders_filtered_filtered_semi_join_projected_open_10_cipher_merge_rows, 1);
            for(size_t i = 0; i < orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count; i++){
                std::memcpy(
                    orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge.mShares[0].data() + i * orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows,
                    orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[0].data(),
                    orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows * sizeof(orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[0](0, 0)));
                std::memcpy(
                    orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge.mShares[1].data() + i * orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows,
                    orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[1].data(),
                    orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows * sizeof(orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[1](0, 0)));
            }
            large_data_sending(role, orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge.mShares[0], runtime, false);
        }

        /*
        for(size_t i = 0; i < orders_filtered_filtered_semi_join_projected_9_cipher_cols.size(); i++){
            large_data_sending(role, orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[0], runtime, false);
        }
        */
    }
    else if(role ==(1 - 1)%3){
        //将cipher_cols恢复成明文（合并通信后再拆分）
        size_t orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count = orders_filtered_filtered_semi_join_projected_9_cipher_cols.size();
        std::vector<i64Matrix> orders_filtered_filtered_semi_join_projected_open_10_cipher_cols(orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count);
        if (orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count > 0) {
            size_t orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows = orders_filtered_filtered_semi_join_projected_9_cipher_cols[0].rows();
            size_t orders_filtered_filtered_semi_join_projected_open_10_cipher_merge_rows = orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count * orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows;

            si64Matrix orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge(orders_filtered_filtered_semi_join_projected_open_10_cipher_merge_rows, 1);
            for(size_t i = 0; i < orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count; i++){
                std::memcpy(
                    orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge.mShares[0].data() + i * orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows,
                    orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[0].data(),
                    orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows * sizeof(orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[0](0, 0)));
                std::memcpy(
                    orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge.mShares[1].data() + i * orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows,
                    orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[1].data(),
                    orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows * sizeof(orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[1](0, 0)));
            }

            i64Matrix orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_recv_merge(orders_filtered_filtered_semi_join_projected_open_10_cipher_merge_rows, 1);
            large_data_receiving(role, orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_recv_merge, runtime, false);
            i64Matrix orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_merge = orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_recv_merge;
            orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_merge += orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge.mShares[0];
            orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_merge += orders_filtered_filtered_semi_join_projected_9_cipher_cols_merge.mShares[1];

            for(size_t i = 0; i < orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_count; i++){
                orders_filtered_filtered_semi_join_projected_open_10_cipher_cols[i].resize(orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows, 1);
                std::memcpy(
                    orders_filtered_filtered_semi_join_projected_open_10_cipher_cols[i].data(),
                    orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_merge.data() + i * orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows,
                    orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_rows * sizeof(orders_filtered_filtered_semi_join_projected_open_10_cipher_cols_merge(0, 0)));
            }
        }

        /*
        for(size_t i = 0; i < orders_filtered_filtered_semi_join_projected_open_10_cipher_cols.size(); i++){
            orders_filtered_filtered_semi_join_projected_open_10_cipher_cols[i].resize(orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].rows(), 1);
            large_data_receiving(role, orders_filtered_filtered_semi_join_projected_open_10_cipher_cols[i], runtime, false);
            orders_filtered_filtered_semi_join_projected_open_10_cipher_cols[i] += orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[0];
            orders_filtered_filtered_semi_join_projected_open_10_cipher_cols[i] += orders_filtered_filtered_semi_join_projected_9_cipher_cols[i].mShares[1];
        }
        */

        //全明文orders_filtered_filtered_semi_join_projected_open_10
        size_t orders_filtered_filtered_semi_join_projected_open_10_size = orders_filtered_filtered_semi_join_projected_9_plain_cols.size() + orders_filtered_filtered_semi_join_projected_9_cipher_cols.size();
        std::vector<i64Matrix> orders_filtered_filtered_semi_join_projected_open_10(orders_filtered_filtered_semi_join_projected_open_10_size);

        //直接加入明文列
        for(size_t i=0 ; i < orders_filtered_filtered_semi_join_projected_9_plain_cols.size() ; i++){
            orders_filtered_filtered_semi_join_projected_open_10[orders_filtered_filtered_semi_join_projected_9_plain_cols_idx[i]] = orders_filtered_filtered_semi_join_projected_9_plain_cols[i];
        }
        //加入恢复的密文列
        for(size_t i=0 ; i < orders_filtered_filtered_semi_join_projected_9_cipher_cols.size() ; i++){
            orders_filtered_filtered_semi_join_projected_open_10[orders_filtered_filtered_semi_join_projected_9_cipher_cols_idx[i]] = orders_filtered_filtered_semi_join_projected_open_10_cipher_cols[i];
        }
        

        //存储在STP下的i64Matrix
        std::ofstream file("/root/plantest/benchmarks/tpch_test/32768/db/1/data/update/orders_filtered_filtered_semi_join_projected_open_10.csv");
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: ");
        }
        for(size_t j = 0; j < orders_filtered_filtered_semi_join_projected_open_10[0].rows(); j++) {
            for(size_t i = 0; i < orders_filtered_filtered_semi_join_projected_open_10.size(); i++) {
                file << orders_filtered_filtered_semi_join_projected_open_10[i](j,0) ;
                if(i < orders_filtered_filtered_semi_join_projected_open_10.size()-1){
                    file << "," ;
                }
            }
            file << std::endl;
        }
        file.close();

    }



    auto t_total_end = std::chrono::high_resolution_clock::now();
    double run_time = std::chrono::duration<double>(t_total_end - t_total_start).count();

    std::cout << "done aby3: " << std::to_string(role) << std::endl;

    if (role == 0 && !timefile.empty()) {
        std::ofstream ofs(timefile);
        if (ofs.is_open()) {
            ofs << std::fixed << std::setprecision(3);
            ofs << "Run Time: " << run_time << " s" << std::endl;
            ofs.close();
        }
    }

    return 0;
}