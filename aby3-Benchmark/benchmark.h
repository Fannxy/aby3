#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>

#include "../aby3-RTR/debug.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/PtATests.h"
#include "../aby3-Basic/timer.h"

// large and custom data transmission.
void create_custom_mpi_type(MPI_Datatype* mpi_type, size_t size);


// generate random data for benchmarking.
void generate_data(aby3::si64Matrix &data, size_t n, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data(std::vector<aby3::si64Matrix> &data, size_t n, size_t m, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data_parallel(std::vector<aby3::si64Matrix> &data, size_t n, size_t m, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data(aby3::sbMatrix &data, size_t n, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data(std::vector<aby3::sbMatrix> &data, size_t n, size_t m, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data(std::vector<aby3::si64> &data, size_t n, size_t m, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data(std::vector<aby3::sb64> &data, size_t n, size_t m, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data_parallel(std::vector<aby3::sbMatrix> &data, size_t n, size_t m, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data_parallel(std::vector<aby3::si64> &data, size_t n, size_t m, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);
void generate_data_parallel(std::vector<aby3::sb64> &data, size_t n, size_t m, size_t bitsize, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

// benchmark functions for search functions.
void binary_search_benchmark(oc::CLP& cmd);
void ndss22_binary_search_benchmark(oc::CLP& cmd);
void pta_binary_search_benchmark(oc::CLP& cmd);