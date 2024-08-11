
#include <cryptoTools/Common/CLP.h>
#include <tests_cryptoTools/UnitTests.h>
#include <map>
#include <mpi.h>
#include "aby3-Benchmark/benchmark.h"
#include "eric.h"

using namespace oc;
using namespace aby3;

int main(int argc, char** argv) {
  oc::CLP cmd(argc, argv);

  if(cmd.isSet("search")){
    binary_search_benchmark(cmd);
  }
  if(cmd.isSet("pta-search")){
    MPI_Init(&argc, &argv);
    pta_binary_search_benchmark(cmd);
    MPI_Finalize();
  }

  return 0;
}