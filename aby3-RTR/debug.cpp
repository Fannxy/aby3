#include "debug.h"
#include "BuildingBlocks.h"
#include <fstream>

void debug_output_vector(std::vector<aby3::si64>& problem_vec, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc){

  aby3::u64 length = (aby3::u64) problem_vec.size();
  aby3::si64Matrix debug(length, 1);
  aby3::i64Matrix plaininfo(length, 1);

  for(int i=0; i<length; i++){
    debug(i, 0, problem_vec[i]);
  }
  enc.revealAll(runtime, debug, plaininfo).get();
  std::ofstream ofs(debugFile, std::ios_base::app);
  for(int i=0; i<length; i++) ofs << plaininfo(i, 0) << " ";
  ofs << std::endl;
  ofs.close();
}

void debug_output_matrix(aby3::si64Matrix& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc){

  aby3::u64 length = problem_mat.rows();
  aby3::i64Matrix plaininfo(length, 1);
  enc.revealAll(runtime, problem_mat, plaininfo).get();

  std::ofstream ofs(debugFile, std::ios_base::app);
  for(int i=0; i<length; i++) ofs << plaininfo(i, 0) << " ";
  ofs << std::endl;
  ofs.close();
}

void debug_output_matrix(aby3::sbMatrix& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc, int pIdx, aby3::Sh3Evaluator& eval){
  aby3::u64 length = problem_mat.rows();
  // aby3::si64Matrix ones(length, 1);
  // init_ones(pIdx, enc, runtime, ones, length);
  // eval.asyncMul(runtime, ones, problem_mat, ones).get();

  aby3::i64Matrix plaininfo(length, 1);
  enc.revealAll(runtime, problem_mat, plaininfo).get();

  std::ofstream ofs(debugFile, std::ios_base::app);
  for(int i=0; i<length; i++) ofs << plaininfo(i, 0) << " ";
  ofs << std::endl;
  ofs.close();
}



void debug_output_matrix(aby3::i64Matrix& problem_mat){
  std::ofstream ofs(debugFile, std::ios_base::app);
  int length = problem_mat.rows();
  int print_length = (length > 100) ? 100 : length;
  for(int i=0; i<print_length; i++) ofs << problem_mat(i, 0) << " ";
  ofs << std::endl;
  ofs.close();
}


void debug_mpi(int rank, int pIdx, std::string info){
  std::string debugFile_mpi = debugFolder + "DEBUG-role:" + std::to_string(pIdx) + "-rank:" + std::to_string(rank) + ".txt";

  std::ofstream ofs(debugFile_mpi, std::ios_base::app);
  ofs << info << std::endl;
  ofs.close();
}

void debug_info(std::string info){
  std::ofstream ofs(debugFile, std::ios_base::app);
  ofs << info << std::endl;
  ofs.close();
}

void debug_info(std::string info, std::ofstream& ofs){
  ofs << info << std::endl;
}

void debug_info(Eigen::internal::Packet4i &info){
  std::ofstream ofs(debugFile, std::ios_base::app);
  ofs << info << std::endl;
  ofs.close();
}

void write_log(std::string log_file, std::string info){
  std::ofstream ofs(log_file, std::ios_base::app);
  ofs << info << std::endl;
  ofs.close();
}


void debug_secret_matrix(aby3::sbMatrix& problem_mat){
  std::ofstream ofs(debugFile, std::ios_base::app);
  int length = problem_mat.rows();

  ofs << "secret 0: " << std::endl;
  for(int i=0; i<length; i++){
    printBits(problem_mat.mShares[0](i), ofs);
  }
  ofs << "secret 1: " << std::endl;
  for(int i=0; i<length; i++){
    printBits(problem_mat.mShares[1](i), ofs);
  }
  ofs << std::endl;
}