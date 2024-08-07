
#include <cryptoTools/Common/CLP.h>
#include <tests_cryptoTools/UnitTests.h>
#include <map>
#include <mpi.h>
#include "aby3-RTR/PtATests.h"
#include "aby3-RTR/PtAProfile.h"
#include "eric.h"
#include "aby3_tests/Test.h"

using namespace oc;
using namespace aby3;

int main(int argc, char** argv) {
  oc::CLP cmd(argc, argv);

  MPI_Init(&argc, &argv);

  if(cmd.isSet("cipher_index")) {
	  test_cipher_index_pta(cmd);
  }

  if(cmd.isSet("max")){
    test_max_pta(cmd);
  }

  if(cmd.isSet("sort")){
    test_sort_pta(cmd);
  }

  if(cmd.isSet("sum")){
    test_sum_pta(cmd);
  }

  if(cmd.isSet("metric")){
    test_metric_pta(cmd);
  }

  if(cmd.isSet("system_profile")){
    // communication_profile(cmd);
    pta_system_profile(cmd);
  }

  if(cmd.isSet("pta_correctness")){
    correctness_cipher_index_pta(cmd);
    MPI_Barrier(MPI_COMM_WORLD);  // Wait for all processes to finish before starting the next test
    correctness_sort_pta(cmd);
    MPI_Barrier(MPI_COMM_WORLD);
    correctness_sum_pta(cmd);
    MPI_Barrier(MPI_COMM_WORLD);
    correctness_max_pta(cmd);
    MPI_Barrier(MPI_COMM_WORLD);
    correctness_metric_pta(cmd);
    MPI_Barrier(MPI_COMM_WORLD);
    correctness_search_pta(cmd);
  }

  if(cmd.isSet("task_profile")){
    std::vector<std::string> keywords_check_list = {"task", "logFolder", "startB", "gap", "endingB"};
    for (auto& keyword : keywords_check_list) {
      if (!cmd.isSet(keyword)) {
        throw std::runtime_error("Missing keyword: " + keyword);
      }
    }
    task_profile(cmd);
  }

  MPI_Finalize();

  return 0;
}