#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <fstream>

#ifndef _DEBUG_H_
#define _DEBUG_H_

#ifndef DEBUG_FILE
#define DEBUG_FILE "/root/debug.txt"
#endif

inline std::string PARTY_FILE = "/root/aby3/party-";

#define MAIN_PARTY_DEBUG(pIdx, debug_info) \ 
if(pIdx == 0) { \
    std::cout << debug_info << std::endl; \
}

#define THROW_RUNTIME_ERROR(msg) \
    throw std::runtime_error("Error in file " + std::string(__FILE__) + " at line " + std::to_string(__LINE__) + ": " + msg)

// static std::string debugFile = "/root/aby3/debug.txt";
static std::string debugFile = DEBUG_FILE;
static std::string debugFolder= "/root/aby3/DEBUG/";

inline void printBits(int64_t num, std::ostream& stream){
    for(int i=63; i>=0; i--){
        stream << ((num >> i) & 1) << " ";
    }
    stream << std::endl;
}

extern void debug_mpi(int rank, int pIdx, std::string info);

extern void debug_info(std::string info);

extern void debug_info(std::string info, std::ofstream& ofs);

extern void debug_info(Eigen::internal::Packet4i &info);

extern void write_log(std::string log_file, std::string info);

extern void debug_output_vector(std::vector<aby3::si64>& problem_vec, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc);

extern void debug_output_matrix(aby3::si64Matrix& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc);

extern void debug_output_matrix(aby3::i64Matrix& problem_mat);

extern void debug_output_matrix(aby3::sbMatrix& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc, int pIdx, aby3::Sh3Evaluator& eval);

template <aby3::Decimal D>
extern void debug_output_matrix(aby3::sf64Matrix<D>& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc){
    aby3::u64 length = problem_mat.rows();    
    aby3::f64Matrix<D> plaininfo(length, 1);
    enc.revealAll(runtime, problem_mat, plaininfo).get();

    std::ofstream ofs(debugFile, std::ios_base::app);
    ofs << "length: " << length << std::endl;
    for(int i=0; i<length; i++) ofs << plaininfo(i, 0) << " ";
    ofs << std::endl;
    ofs.close();
}

template <aby3::Decimal D>
extern void debug_output_matrix(aby3::f64Matrix<D>& problem_mat){
    aby3::u64 length = problem_mat.rows();    
    std::ofstream ofs(debugFile, std::ios_base::app);
    ofs << "length: " << length << std::endl;
    for(int i=0; i<length; i++) ofs << problem_mat(i, 0) << " ";
    ofs << std::endl;
    ofs.close();
}

template <aby3::Decimal D>
extern void debug_output_vector(std::vector<aby3::sf64<D>>& problem_vec, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc){
    aby3::u64 length = problem_vec.size();
    aby3::sf64Matrix<D> problem_mat; problem_mat.resize(length, 1);
    for(int i=0; i<length; i++) problem_mat(i, 0, problem_vec[i]);
    return debug_output_matrix<D>(problem_mat, runtime, enc); 
}

template<typename T>
extern void debug_output_vector(std::vector<T> &problem_vec){
    std::ofstream ofs(debugFile, std::ios_base::app);
    ofs << "length: " << problem_vec.size() << std::endl;
    int print_length = (problem_vec.size() > 25) ? 25 : problem_vec.size();
    for(int i=0; i<print_length; i++) ofs << problem_vec[i] << " ";
    ofs << std::endl;
    ofs.close();
}

template<typename T>
extern void debug_output_vector(std::vector<T> &problem_vec, std::ofstream& ofs){
    ofs << "length: " << problem_vec.size() << std::endl;
    for(int i=0; i<problem_vec.size(); i++) ofs << problem_vec[i] << " ";
    ofs << std::endl;
}

template<typename T>
extern void debug_output_value(T &problem_val){
    std::ofstream ofs(debugFile, std::ios_base::app);
    ofs << "value: " << problem_val << std::endl;
    ofs.close();
}

extern void debug_secret_matrix(aby3::sbMatrix& problem_mat);

#endif