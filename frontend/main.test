
#include <cryptoTools/Common/CLP.h>
#include <tests_cryptoTools/UnitTests.h>
#include <map>
#include <mpi.h>
#include "aby3_tests/Test.h"
#include "aby3_tests/aby3_tests.h"
#include "eric.h"

using namespace oc;
using namespace aby3;

int main(int argc, char** argv) {
  oc::CLP cmd(argc, argv);
  // reinit the environment and then finalize the environment.

  // set the role for this process.
	if (cmd.isSet("Bool")){
		bool_basic_test(cmd);
		bool_basic_test2(cmd);
		get_first_zero_test(cmd);
		bool_aggregation_test(cmd);
		share_conversion_test(cmd);
	}

	if (cmd.isSet("Arith")){
		arith_basic_test(cmd);
	}

	if (cmd.isSet("Init")){
		initialization_test(cmd);
		correlation_test(cmd);
	}

	if(cmd.isSet("Comm")){
		communication_test(cmd);
	}

	if(cmd.isSet("Shuffle")){
		shuffle_test(cmd);
		// large_scale_shuffle_test(cmd); // large-scale test, costs at least 150GB+ memory.
	}

	if(cmd.isSet("ORAM")){
		pos_map_test(cmd);
		sqrt_oram_test(cmd);
	}

	if(cmd.isSet("Sort")){
		bc_sort_test(cmd);
		bc_sort_corner_test(cmd);
		bc_sort_multiple_times(cmd);
		quick_sort_test(cmd);
		// quick_sort_with_duplicate_elements_test(cmd); // too slow
		odd_even_merge_test(cmd);
	}

	if(cmd.isSet("Search")){
		// constant_dot_test(cmd);
		binary_search_test(cmd);
	}

  return 0;
}