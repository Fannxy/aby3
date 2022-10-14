
#include <cryptoTools/Common/CLP.h>
#include <map>

#include "aby3_tests/aby3_tests.h"
#include "eric.h"
#include "aby3-DB-main.h"
#include "aby3-DB_tests/UnitTests.h"
#include <tests_cryptoTools/UnitTests.h>
#include <aby3-ML/main-linear.h>
#include <aby3-ML/main-logistic.h>
#include "aby3-RTR/RTRTest.h"
#include "aby3-RTR/DistributeRTRTest.h"

#include "tests_cryptoTools/UnitTests.h"
#include "cryptoTools/Crypto/PRNG.h"

using namespace oc;
using namespace aby3;
std::vector<std::string> unitTestTag{ "u", "unitTest" };

// #define BASIC_TEST
// #define PERFORMANCE_TEST
// #define DISTRIBUTE_TEST
#define DISTRIBUTE_PERFORMANCE

void help()
{

	std::cout << "-u                        ~~ to run all tests" << std::endl;
	std::cout << "-u n1 [n2 ...]            ~~ to run test n1, n2, ..." << std::endl;
	std::cout << "-u -list                  ~~ to list all tests" << std::endl;
	std::cout << "-intersect -nn NN [-c C]  ~~ to run the intersection benchmark with 2^NN set sizes, C 32-bit data columns." << std::endl;
	std::cout << "-eric -nn NN              ~~ to run the eric benchmark with 2^NN set sizes" << std::endl;
	std::cout << "-threat -nn NN -s S       ~~ to run the threat log benchmark with 2^NN set sizes and S sets" << std::endl;
}


// int main(int argc, char** argv)
// {


// 	try {


// 		bool set = false;
// 		oc::CLP cmd(argc, argv);


// 		if (cmd.isSet(unitTestTag))
// 		{
// 			auto tests = tests_cryptoTools::Tests;
// 			tests += aby3_tests;
// 			tests += DB_tests;

// 			tests.runIf(cmd);
// 			return 0;
// 		}

// 		if (cmd.isSet("RTR"))
// 		{
// 			set = true;
// 			test_mul(cmd);
// 			std::cout << "can call functions " << std::endl;
// 		}

// 		if (cmd.isSet("linear-plain"))
// 		{
// 			set = true;
// 			linear_plain_main(cmd);
// 		}
// 		if (cmd.isSet("linear"))
// 		{
// 			set = true;
// 			linear_main_3pc_sh(cmd);
// 		}

// 		if (cmd.isSet("logistic-plain"))
// 		{
// 			set = true;
// 			logistic_plain_main(cmd);
// 		}

// 		if (cmd.isSet("logistic"))
// 		{
// 			set = true;
// 			logistic_main_3pc_sh(cmd);
// 		}

// 		if (cmd.isSet("eric"))
// 		{
// 			set = true;

// 			auto nn = cmd.getMany<int>("nn");
// 			if (nn.size() == 0)
// 				nn.push_back(16);

// 			for (auto n : nn)
// 			{
// 				eric(1 << n);
// 			}
// 		}


// 		if (cmd.isSet("intersect"))
// 		{
// 			set = true;

// 			auto nn = cmd.getMany<int>("nn");
// 			auto c = cmd.getOr("c", 0);
// 			if (nn.size() == 0)
// 				nn.push_back(1 << 16);

// 			for (auto n : nn)
// 			{
// 				auto size = 1 << n;
// 				DB_Intersect(size, c, cmd.isSet("sum"));
// 			}
// 		}


// 		if (cmd.isSet("threat"))
// 		{
// 			set = true;

// 			auto nn = cmd.getMany<int>("nn");
// 			auto c = cmd.getOr("s", 2);
// 			if (nn.size() == 0)
// 				nn.push_back(1 << 16);

// 			for (auto n : nn)
// 			{
// 				auto size = 1 << n;
// 				DB_threat(size, c);
// 			}
// 		}



// 		if (cmd.isSet("card"))
// 		{
// 			set = true;

// 			auto nn = cmd.getMany<int>("nn");
// 			if (nn.size() == 0)
// 				nn.push_back(1 << 16);

// 			for (auto n : nn)
// 			{
// 				auto size = 1 << n;
// 				DB_cardinality(size);
// 			}
// 		}
		
// 		//if (cmd.isSet("add"))
// 		//{
// 		//	set = true;

// 		//	auto nn = cmd.getMany<int>("nn");
// 		//	if (nn.size() == 0)
// 		//		nn.push_back(1 << 16);

// 		//	for (auto n : nn)
// 		//	{
// 		//		auto size = 1 << n;
// 		//		Sh3_add_test(size);
// 		//	}
// 		//}

// 		if (set == false)
// 		{
// 			help();
// 		}

// 	}
// 	catch (std::exception& e)
// 	{
// 		std::cout << e.what() << std::endl;
// 	}

// 	return 0;
// }

int main(int argc, char** argv)
{
	oc::CLP cmd(argc, argv);

	#ifdef PERFORMANCE_TEST
	// test the vectorization for basic ops (mul) and (gt).
	int repeats = int(100);

	std::vector<int> n_list = {      10,       13,       17,       23,       30,       40,
             54,       71,       95,      126,      167,      222,
            294,      390,      517,      686,      910,     1206,
           1599,     2120,     2811,     3727,     4941,     6551,
           8685,    11513,    15264,    20235,    26826,    35564,
          47148,    62505,    82864,   109854,   145634,   193069,
         255954,   339322,   449843,   596362,   790604,  1048113,
        1389495,  1842069,  2442053,  3237457,  4291934,  5689866,
        7543120, 10000000};

	std::map<int, std::map<std::string, std::vector<double>>> performance_dict;
	for(int i=0; i<n_list.size(); i++){
		int n = n_list[i];
		std::map<std::string, std::vector<double>> tmp_map;
		basic_performance(cmd, n, repeats, tmp_map);
		// dis_basic_performance(cmd, n, repeats, tmp_map);
		performance_dict[n] = tmp_map;

		// execute one evaluation and record one.
		std::cout << "\nvector size = " << n_list[i] << std::endl;
		std::cout << "mul" << std::endl;
		for(int j=0; j<3; j++){
			std::cout << performance_dict[n_list[i]]["mul"][j] << " ";
		}
		std::cout << "\ngt" << std::endl;
		for(int j=0; j<3; j++){
			std::cout << performance_dict[n_list[i]]["gt"][j] << " ";
		}
		// std::cout << std::endl;
		std::cout << "\nadd" << std::endl;
		for(int j=0; j<3; j++){
			std::cout << performance_dict[n_list[i]]["add"][j] << " ";
		}
		std::cout << std::endl;
	}

	// cout the result.
	for(int i=0; i<n_list.size(); i++){
		std::cout << "\nvector size = " << n_list[i] << std::endl;
		std::cout << "mul" << std::endl;
		for(int j=0; j<3; j++){
			std::cout << performance_dict[n_list[i]]["mul"][j] << " ";
		}
		std::cout << "\ngt" << std::endl;
		for(int j=0; j<3; j++){
			std::cout << performance_dict[n_list[i]]["gt"][j] << " ";
		}
		// std::cout << std::endl;
		std::cout << "\nadd" << std::endl;
		for(int j=0; j<3; j++){
			std::cout << performance_dict[n_list[i]]["add"][j] << " ";
		}
		std::cout << std::endl;
	}
	#endif

	#ifdef BASIC_TEST
	// test gt
	// test_gt(cmd);

	// test eq - has problems.
	// test_eq(cmd);

	// test multiplication between bits and ints.
	// test_mul(cmd);

	// test cipher_argsort
	test_argsort(cmd, 1);
	test_argsort(cmd, 0);

	// test cipher_index
	test_cipher_index(cmd, 0);
	test_cipher_index(cmd, 1);

	// test binning.
	test_cipher_binning(cmd, 0);
	test_cipher_binning(cmd, 1);
	#endif

	#ifdef DISTRIBUTE_TEST
	// basic function test
	test_mul(cmd);
	dis_test_mul(cmd);
	#endif

	#ifdef DISTRIBUTE_PERFORMANCE
	int repeats;
	std::vector<int> n_list = {      10,       13,       17,       23,       30,       40,
				54,       71,       95,      126,      167,      222,
			294,      390,      517,      686,      910,     1206,
			1599,     2120,     2811,     3727,     4941,     6551,
			8685,    11513,    15264,    20235,    26826,    35564,
			47148,    62505,    82864,   109854,   145634,   193069,
			255954,   339322,   449843,   596362,   790604,  1048113,
		1389495,  1842069,  2442053,  3237457,  4291934,  5689866,
		7543120, 10000000};
	// std::vector<int> n_list = {3237457,  4291934,  5689866, 7543120, 10000000};

	std::map<int, std::map<std::string, double>> performance_dict;
	for(int i=0; i<n_list.size(); i++){

		int n = n_list[i];
		
		// set the repeat times.
		if(i < 20) repeats = int(1e4);
		else if(i < 40) repeats = int(1e3);
		else repeats = int(100);

		std::map<std::string, double> tmp_map;
		dis_basic_performance(cmd, n, repeats, tmp_map);
		performance_dict[n] = tmp_map;

		// execute one evaluation and record one.
		std::cout << "\nvector size = " << n_list[i] << std::endl;
		std::map<std::string, double>::iterator iter;
		iter = tmp_map.begin();
		while(iter != tmp_map.end()){
			std::cout << iter->first << std::endl;
			std::cout << iter->second << std::endl;
			iter ++;
		}
	}
	#endif

	return 0;
}
