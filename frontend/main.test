#include <map>
#include <chrono>
#include <random>
#include <thread>
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

    std::cout << "start aby3: " << std::to_string(role) << std::endl;

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

//create
    std::vector<i64Matrix> eq_flags_hybrid_agg_1(1);
    i64Matrix eq_flags_hybrid_agg_1_unit_size(1, 1);
    if (role == (1 - 1)%3 ){
        std::string path = "/root/plantest/benchmarks/ssn/db/1/data/update/eq_flags_hybrid_agg_1.csv" ;
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }

        std::string eq_flags_hybrid_agg_1_line;
        std::vector<std::vector<i64>> eq_flags_hybrid_agg_1_data_rows;
        
        //read data from path
        while (std::getline(file, eq_flags_hybrid_agg_1_line)) {
            if (eq_flags_hybrid_agg_1_line.empty()) continue;
            
            std::vector<i64> row;
            std::istringstream iss(eq_flags_hybrid_agg_1_line);
            std::string token;
            
            while (std::getline(iss, token, ',')) {
                if (!token.empty()) {
                    row.push_back(std::stoll(token));
                }
            }
            
            if (!row.empty()) {
                eq_flags_hybrid_agg_1_data_rows.push_back(row);
            }
        }
        file.close();
        
        if (eq_flags_hybrid_agg_1_data_rows.empty()) {
            throw std::runtime_error("No data found in file: " + path);
        }
        
        size_t eq_flags_hybrid_agg_1_num_cols = eq_flags_hybrid_agg_1_data_rows[0].size();
        size_t eq_flags_hybrid_agg_1_num_rows = eq_flags_hybrid_agg_1_data_rows.size();

        
        eq_flags_hybrid_agg_1_unit_size(0, 0) = eq_flags_hybrid_agg_1_num_rows ;
        large_data_sending(role, eq_flags_hybrid_agg_1_unit_size, runtime, true);
        large_data_sending(role, eq_flags_hybrid_agg_1_unit_size, runtime, false);
        
        for (size_t col = 0; col < eq_flags_hybrid_agg_1_num_cols; col++) {
            eq_flags_hybrid_agg_1[col].resize(eq_flags_hybrid_agg_1_num_rows, 1); 
            for (size_t row = 0; row < eq_flags_hybrid_agg_1_num_rows; row++) {
                eq_flags_hybrid_agg_1[col](row, 0) = eq_flags_hybrid_agg_1_data_rows[row][col];
            }
        }
        
    }
    else if(role == (1)%3){
        large_data_receiving(role, eq_flags_hybrid_agg_1_unit_size, runtime, true);
    }
    else if(role == (1 + 1)%3){
        large_data_receiving(role, eq_flags_hybrid_agg_1_unit_size, runtime, false);
    }
    
//close
    std::vector<si64Matrix> closed_eq_flags_hybrid_agg_1(1);

    for (size_t i = 0; i< 1; i++){
        closed_eq_flags_hybrid_agg_1[i].resize(eq_flags_hybrid_agg_1_unit_size(0, 0), 1);
        if (role == (1 - 1)%3){
            enc.localIntMatrix(runtime, eq_flags_hybrid_agg_1[i], closed_eq_flags_hybrid_agg_1[i]).get();
        } else {
            enc.remoteIntMatrix(runtime, closed_eq_flags_hybrid_agg_1[i]).get();
        }
    }

//create
    std::vector<i64Matrix> sorted_by_key_dummy_hybrid_agg_1(1);
    i64Matrix sorted_by_key_dummy_hybrid_agg_1_unit_size(1, 1);
    if (role == (1 - 1)%3 ){
        std::string path = "/root/plantest/benchmarks/ssn/db/1/data/update/sorted_by_key_dummy_hybrid_agg_1.csv" ;
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }

        std::string sorted_by_key_dummy_hybrid_agg_1_line;
        std::vector<std::vector<i64>> sorted_by_key_dummy_hybrid_agg_1_data_rows;
        
        //read data from path
        while (std::getline(file, sorted_by_key_dummy_hybrid_agg_1_line)) {
            if (sorted_by_key_dummy_hybrid_agg_1_line.empty()) continue;
            
            std::vector<i64> row;
            std::istringstream iss(sorted_by_key_dummy_hybrid_agg_1_line);
            std::string token;
            
            while (std::getline(iss, token, ',')) {
                if (!token.empty()) {
                    row.push_back(std::stoll(token));
                }
            }
            
            if (!row.empty()) {
                sorted_by_key_dummy_hybrid_agg_1_data_rows.push_back(row);
            }
        }
        file.close();
        
        if (sorted_by_key_dummy_hybrid_agg_1_data_rows.empty()) {
            throw std::runtime_error("No data found in file: " + path);
        }
        
        size_t sorted_by_key_dummy_hybrid_agg_1_num_cols = sorted_by_key_dummy_hybrid_agg_1_data_rows[0].size();
        size_t sorted_by_key_dummy_hybrid_agg_1_num_rows = sorted_by_key_dummy_hybrid_agg_1_data_rows.size();

        
        sorted_by_key_dummy_hybrid_agg_1_unit_size(0, 0) = sorted_by_key_dummy_hybrid_agg_1_num_rows ;
        large_data_sending(role, sorted_by_key_dummy_hybrid_agg_1_unit_size, runtime, true);
        large_data_sending(role, sorted_by_key_dummy_hybrid_agg_1_unit_size, runtime, false);
        
        for (size_t col = 0; col < sorted_by_key_dummy_hybrid_agg_1_num_cols; col++) {
            sorted_by_key_dummy_hybrid_agg_1[col].resize(sorted_by_key_dummy_hybrid_agg_1_num_rows, 1); 
            for (size_t row = 0; row < sorted_by_key_dummy_hybrid_agg_1_num_rows; row++) {
                sorted_by_key_dummy_hybrid_agg_1[col](row, 0) = sorted_by_key_dummy_hybrid_agg_1_data_rows[row][col];
            }
        }
        
    }
    else if(role == (1)%3){
        large_data_receiving(role, sorted_by_key_dummy_hybrid_agg_1_unit_size, runtime, true);
    }
    else if(role == (1 + 1)%3){
        large_data_receiving(role, sorted_by_key_dummy_hybrid_agg_1_unit_size, runtime, false);
    }
    
//close
    std::vector<si64Matrix> closed_sorted_by_key_hybrid_agg_1(1);

    for (size_t i = 0; i< 1; i++){
        closed_sorted_by_key_hybrid_agg_1[i].resize(sorted_by_key_dummy_hybrid_agg_1_unit_size(0, 0), 1);
        if (role == (1 - 1)%3){
            enc.localIntMatrix(runtime, sorted_by_key_dummy_hybrid_agg_1[i], closed_sorted_by_key_hybrid_agg_1[i]).get();
        } else {
            enc.remoteIntMatrix(runtime, closed_sorted_by_key_hybrid_agg_1[i]).get();
        }
    }

//index_agg
    i64Matrix closed_sorted_by_key_hybrid_agg_1_plain(closed_sorted_by_key_hybrid_agg_1[0].rows(), 1);
    enc.revealAll(runtime, closed_sorted_by_key_hybrid_agg_1[0], closed_sorted_by_key_hybrid_agg_1_plain).get();

    std::string persisted_hybrid_agg_1_name = "persisted_hybrid_agg_1";
    std::vector<si64Matrix> persisted_hybrid_agg_1;
    read_cipher(role, persisted_hybrid_agg_1_name, persisted_hybrid_agg_1);

    std::vector<si64Matrix> persisted_hybrid_agg_1_for_agg(2);
    persisted_hybrid_agg_1_for_agg[0] = persisted_hybrid_agg_1[1];
    persisted_hybrid_agg_1_for_agg[1] = persisted_hybrid_agg_1[3];

    std::vector<si64Matrix> groupedjoined(2);
    index_agg(role, closed_eq_flags_hybrid_agg_1[0], closed_sorted_by_key_hybrid_agg_1_plain, persisted_hybrid_agg_1_for_agg, groupedjoined,
                enc, eval, runtime);

//open
    if(role == (1) % 3){
        for(size_t i = 0; i < groupedjoined.size(); i++){
            large_data_sending(role, groupedjoined[i].mShares[0], runtime, false);
        }
    }
    else if(role ==(1 - 1)%3){
        std::vector<i64Matrix> groupedjoined_open(groupedjoined.size());
        for(size_t i = 0; i < groupedjoined_open.size(); i++){
            groupedjoined_open[i].resize(groupedjoined[i].rows(), 1);
            large_data_receiving(role, groupedjoined_open[i], runtime, false);
            groupedjoined_open[i] += groupedjoined[i].mShares[0];
            groupedjoined_open[i] += groupedjoined[i].mShares[1];
        }

        //存储在STP下的i64Matrix
        std::ofstream file("/root/plantest/benchmarks/ssn/db/1/data/update/groupedjoined_open.csv");
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: ");
        }
        for(size_t j = 0; j < groupedjoined_open[0].rows(); j++) {
            for(size_t i = 0; i < groupedjoined_open.size(); i++) {
                file << groupedjoined_open[i](j,0) ;
                if(i < groupedjoined_open.size()-1){
                    file << "," ;
                }
            }
            file << std::endl;
        }
        file.close();

    }



    std::cout << "done aby3: " << std::to_string(role) << std::endl;

    return 0;
}