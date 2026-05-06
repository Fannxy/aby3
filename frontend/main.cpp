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
    std::vector<i64Matrix> PostTypes_0(2);
    i64Matrix PostTypes_0_unit_size(1, 1);
    if (role == (1 - 1)%3 ){
        std::string path = "/root/plantest/benchmarks/stackoverflow/32768/db/1/data/update/PostTypes_0.csv" ;
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }

        std::string PostTypes_0_line;
        std::vector<std::vector<i64>> PostTypes_0_data_rows;
        
        //read data from path
        while (std::getline(file, PostTypes_0_line)) {
            if (PostTypes_0_line.empty()) continue;
            
            std::vector<i64> row;
            std::istringstream iss(PostTypes_0_line);
            std::string token;
            
            while (std::getline(iss, token, ',')) {
                if (!token.empty()) {
                    row.push_back(std::stoll(token));
                }
            }
            
            if (!row.empty()) {
                PostTypes_0_data_rows.push_back(row);
            }
        }
        file.close();
        
        if (PostTypes_0_data_rows.empty()) {
            throw std::runtime_error("No data found in file: " + path);
        }
        
        size_t PostTypes_0_num_cols = PostTypes_0_data_rows[0].size();
        size_t PostTypes_0_num_rows = PostTypes_0_data_rows.size();

        
        PostTypes_0_unit_size(0, 0) = PostTypes_0_num_rows ;
        large_data_sending(role, PostTypes_0_unit_size, runtime, true);
        large_data_sending(role, PostTypes_0_unit_size, runtime, false);
        
        for (size_t col = 0; col < PostTypes_0_num_cols; col++) {
            PostTypes_0[col].resize(PostTypes_0_num_rows, 1); 
            for (size_t row = 0; row < PostTypes_0_num_rows; row++) {
                PostTypes_0[col](row, 0) = PostTypes_0_data_rows[row][col];
            }
        }
        
    }
    else if(role == (1)%3){
        large_data_receiving(role, PostTypes_0_unit_size, runtime, true);
    }
    else if(role == (1 + 1)%3){
        large_data_receiving(role, PostTypes_0_unit_size, runtime, false);
    }
    //close
    std::vector<si64Matrix> PostTypes_close_to_joined_1(2);

    for (size_t i = 0; i< 2; i++){
        PostTypes_close_to_joined_1[i].resize(PostTypes_0_unit_size(0, 0), 1);
        if (role == (1 - 1)%3){
            enc.localIntMatrix(runtime, PostTypes_0[i], PostTypes_close_to_joined_1[i]).get();
        } else {
            enc.remoteIntMatrix(runtime, PostTypes_close_to_joined_1[i]).get();
        }
    }
    //create
    std::vector<i64Matrix> Posts_2(2);
    i64Matrix Posts_2_unit_size(1, 1);
    if (role == (2 - 1)%3 ){
        std::string path = "/root/plantest/benchmarks/stackoverflow/32768/db/2/data/update/Posts_2.csv" ;
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }

        std::string Posts_2_line;
        std::vector<std::vector<i64>> Posts_2_data_rows;
        
        //read data from path
        while (std::getline(file, Posts_2_line)) {
            if (Posts_2_line.empty()) continue;
            
            std::vector<i64> row;
            std::istringstream iss(Posts_2_line);
            std::string token;
            
            while (std::getline(iss, token, ',')) {
                if (!token.empty()) {
                    row.push_back(std::stoll(token));
                }
            }
            
            if (!row.empty()) {
                Posts_2_data_rows.push_back(row);
            }
        }
        file.close();
        
        if (Posts_2_data_rows.empty()) {
            throw std::runtime_error("No data found in file: " + path);
        }
        
        size_t Posts_2_num_cols = Posts_2_data_rows[0].size();
        size_t Posts_2_num_rows = Posts_2_data_rows.size();

        
        Posts_2_unit_size(0, 0) = Posts_2_num_rows ;
        large_data_sending(role, Posts_2_unit_size, runtime, true);
        large_data_sending(role, Posts_2_unit_size, runtime, false);
        
        for (size_t col = 0; col < Posts_2_num_cols; col++) {
            Posts_2[col].resize(Posts_2_num_rows, 1); 
            for (size_t row = 0; row < Posts_2_num_rows; row++) {
                Posts_2[col](row, 0) = Posts_2_data_rows[row][col];
            }
        }
        
    }
    else if(role == (2)%3){
        large_data_receiving(role, Posts_2_unit_size, runtime, true);
    }
    else if(role == (2 + 1)%3){
        large_data_receiving(role, Posts_2_unit_size, runtime, false);
    }
    //close
    std::vector<si64Matrix> Posts_close_to_joined_3(2);

    for (size_t i = 0; i< 2; i++){
        Posts_close_to_joined_3[i].resize(Posts_2_unit_size(0, 0), 1);
        if (role == (2 - 1)%3){
            enc.localIntMatrix(runtime, Posts_2[i], Posts_close_to_joined_3[i]).get();
        } else {
            enc.remoteIntMatrix(runtime, Posts_close_to_joined_3[i]).get();
        }
    }
    //join
    std::vector<si64Matrix> joined_4;

    std::vector<int> Posts_close_to_joined_3_join_idxs = {1} ;
    std::vector<int> PostTypes_close_to_joined_1_join_idxs = {0} ;
    std::vector<int> Posts_close_to_joined_3_other_idxs, PostTypes_close_to_joined_1_other_idxs;

    if(Posts_close_to_joined_3_join_idxs.size() != PostTypes_close_to_joined_1_join_idxs.size()){
        throw std::runtime_error("join key_cols is not consistent") ;
    }

    for(size_t i = 0; i < Posts_close_to_joined_3.size(); i++){
        if (std::find(Posts_close_to_joined_3_join_idxs.begin(),
                        Posts_close_to_joined_3_join_idxs.end(),
                        static_cast<int>(i)) == Posts_close_to_joined_3_join_idxs.end()) {
            Posts_close_to_joined_3_other_idxs.push_back(static_cast<int>(i));
        }
    }

    for(size_t i = 0; i < PostTypes_close_to_joined_1.size(); i++){
        if (std::find(PostTypes_close_to_joined_1_join_idxs.begin(),
                        PostTypes_close_to_joined_1_join_idxs.end(),
                        static_cast<int>(i)) == PostTypes_close_to_joined_1_join_idxs.end()) {
            PostTypes_close_to_joined_1_other_idxs.push_back(static_cast<int>(i));
        }
    }

    //分别得到两方的key和other列
    std::vector<si64Matrix> Posts_close_to_joined_3_key, Posts_close_to_joined_3_other;
    std::vector<si64Matrix> PostTypes_close_to_joined_1_key, PostTypes_close_to_joined_1_other;
    for(size_t i = 0; i < Posts_close_to_joined_3_join_idxs.size(); i++){
        Posts_close_to_joined_3_key.push_back(Posts_close_to_joined_3[Posts_close_to_joined_3_join_idxs[i]]);
    }
    for(size_t i = 0; i < Posts_close_to_joined_3_other_idxs.size(); i++){
        Posts_close_to_joined_3_other.push_back(Posts_close_to_joined_3[Posts_close_to_joined_3_other_idxs[i]]);
    }
    for(size_t i = 0; i < PostTypes_close_to_joined_1_join_idxs.size(); i++){
        PostTypes_close_to_joined_1_key.push_back(PostTypes_close_to_joined_1[PostTypes_close_to_joined_1_join_idxs[i]]);
    }
    for(size_t i = 0; i < PostTypes_close_to_joined_1_other_idxs.size(); i++){
        PostTypes_close_to_joined_1_other.push_back(PostTypes_close_to_joined_1[PostTypes_close_to_joined_1_other_idxs[i]]);
    }


    join(role, Posts_close_to_joined_3_key, Posts_close_to_joined_3_other,
    PostTypes_close_to_joined_1_key, PostTypes_close_to_joined_1_other,
    joined_4, enc, eval, runtime);

     //agg_count
   std::vector<si64Matrix> grouped_PostCount_joined_5;

   std::vector<int> joined_4_group_idxs = {2} ;
    
   std::vector<si64Matrix> joined_4_group_cols;
   std::vector<si64Matrix> grouped_PostCount_joined_5_group_cols;
   si64Matrix grouped_PostCount_joined_5_cnt_col;

   for(size_t i = 0; i < joined_4_group_idxs.size(); i++){
      joined_4_group_cols.push_back(joined_4[joined_4_group_idxs[i]]);
   }
   group_count(role, joined_4_group_cols, grouped_PostCount_joined_5_group_cols,
                  grouped_PostCount_joined_5_cnt_col, enc, eval, runtime);

   grouped_PostCount_joined_5 = grouped_PostCount_joined_5_group_cols;
   grouped_PostCount_joined_5.push_back(grouped_PostCount_joined_5_cnt_col);                

//project
    std::vector<int> grouped_PostCount_joined_projected_6_select_idx = {0, 1} ;
    std::vector<si64Matrix> grouped_PostCount_joined_projected_6(grouped_PostCount_joined_projected_6_select_idx.size()) ;

    for (size_t i = 0; i< grouped_PostCount_joined_projected_6_select_idx.size(); i++){
        grouped_PostCount_joined_projected_6[i].resize(grouped_PostCount_joined_5[grouped_PostCount_joined_projected_6_select_idx[i]].rows(), 1) ;
        grouped_PostCount_joined_projected_6[i] = grouped_PostCount_joined_5[grouped_PostCount_joined_projected_6_select_idx[i]] ; 
    }
    //open
    if(role == (1) % 3){
        for(size_t i = 0; i < grouped_PostCount_joined_projected_6.size(); i++){
            large_data_sending(role, grouped_PostCount_joined_projected_6[i].mShares[0], runtime, false);
        }
    }
    else if(role ==(1 - 1)%3){
        std::vector<i64Matrix> grouped_PostCount_joined_projected_open_7(grouped_PostCount_joined_projected_6.size());
        for(size_t i = 0; i < grouped_PostCount_joined_projected_open_7.size(); i++){
            grouped_PostCount_joined_projected_open_7[i].resize(grouped_PostCount_joined_projected_6[i].rows(), 1);
            large_data_receiving(role, grouped_PostCount_joined_projected_open_7[i], runtime, false);
            grouped_PostCount_joined_projected_open_7[i] += grouped_PostCount_joined_projected_6[i].mShares[0];
            grouped_PostCount_joined_projected_open_7[i] += grouped_PostCount_joined_projected_6[i].mShares[1];
        }

        //存储在STP下的i64Matrix
        std::ofstream file("/root/plantest/benchmarks/stackoverflow/32768/db/1/data/update/grouped_PostCount_joined_projected_open_7.csv");
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: ");
        }
        for(size_t j = 0; j < grouped_PostCount_joined_projected_open_7[0].rows(); j++) {
            for(size_t i = 0; i < grouped_PostCount_joined_projected_open_7.size(); i++) {
                file << grouped_PostCount_joined_projected_open_7[i](j,0) ;
                if(i < grouped_PostCount_joined_projected_open_7.size()-1){
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