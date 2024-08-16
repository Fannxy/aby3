#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>

#include "../aby3-RTR/debug.h"

#define BASIC_TEST_INIT \ 
    int role = -1; \
    if (cmd.isSet("role")) { \
        auto keys = cmd.getMany<int>("role"); \
        role = keys[0]; \
    } \
    if (role == -1) { \
        throw std::runtime_error(LOCATION); \
    } \
    IOService ios; \
    Sh3Encryptor enc; \
    Sh3Evaluator eval; \
    Sh3Runtime runtime; \
    basic_setup((u64)role, ios, enc, eval, runtime); \

#define SPLITTED_TEST_INIT \ 
    int role = -1; \
    if (cmd.isSet("role")) { \
        auto keys = cmd.getMany<int>("role"); \
        role = keys[0]; \
    } \
    if (role == -1) { \
        throw std::runtime_error(LOCATION); \
    } \
    int rank = -1; \
    if (cmd.isSet("rank")) { \
        auto keys = cmd.getMany<int>("rank"); \
        rank = keys[0]; \
    } \
    if (rank == -1) { \
        throw std::runtime_error(LOCATION); \
    } \
    IOService ios; \
    Sh3Encryptor enc; \
    Sh3Evaluator eval; \
    Sh3Runtime runtime; \
{ \
    std::string p0_ip = "", p1_ip = ""; \
    if (cmd.isSet("p0_ip")) { \
        auto keys = cmd.getMany<std::string>("p0_ip"); \
        p0_ip = keys[0]; \
    } \
    if (p0_ip == "") { \
        throw std::runtime_error(LOCATION); \
    } \
    if (cmd.isSet("p1_ip")) { \
        auto keys = cmd.getMany<std::string>("p1_ip"); \
        p1_ip = keys[0]; \
    } \
    if (p1_ip == "") { \
        throw std::runtime_error(LOCATION); \
    } \
    splitted_setup((u64)role, rank, ios, enc, eval, runtime, p0_ip, p1_ip); \
} \


#define TEST_INIT \
    int role = -1; \
    if (cmd.isSet("role")) { \
        auto keys = cmd.getMany<int>("role"); \
        role = keys[0]; \
    } \
    if (role == -1) { \
        throw std::runtime_error(LOCATION); \
    } \
    IOService ios; \
    Sh3Encryptor enc; \
    Sh3Evaluator eval; \
    Sh3Runtime runtime; \
    basic_setup((u64)role, ios, enc, eval, runtime); \
    aby3Info party_info = aby3Info(role, enc, eval, runtime);


#define SHOW_TEST_CASE

int arith_basic_test(oc::CLP& cmd);
int arith_basic_test2(oc::CLP& cmd);
int bool_basic_test(oc::CLP& cmd);
int bool_basic_test2(oc::CLP& cmd);
int bool_aggregation_test(oc::CLP& cmd);
int get_first_zero_test(oc::CLP& cmd);
int share_conversion_test(oc::CLP& cmd);

int initialization_test(oc::CLP& cmd);
int shuffle_test(oc::CLP& cmd);
int large_scale_shuffle_test(oc::CLP& cmd);
int correlation_test(oc::CLP& cmd);
int communication_test(oc::CLP& cmd);
int splitted_shuffle_test(oc::CLP& cmd);

// oram tests
int pos_map_test(oc::CLP& cmd);
int sqrt_oram_test(oc::CLP& cmd);

// graph tests
// 1) Graph2D test.
int graph_loading_test(oc::CLP& cmd);
int graph_block_fetch_test(oc::CLP& cmd);
int basic_graph_query_test(oc::CLP& cmd);
int neighbors_find_test(oc::CLP& cmd);
int graph_sort_test(oc::CLP &cmd);

// 2) AdjGraph test.
int adj_graph_loading_test(oc::CLP& cmd);
int adj_basic_graph_query_test(oc::CLP& cmd);

// 3) node-edge list test.
int node_edge_list_basic_graph_query_test(oc::CLP& cmd);

// sort tests
int bc_sort_test(oc::CLP& cmd);
int bc_sort_corner_test(oc::CLP& cmd);
int bc_sort_multiple_times(oc::CLP& cmd);
int quick_sort_test(oc::CLP& cmd);
int quick_sort_with_duplicate_elements_test(oc::CLP& cmd);
int odd_even_merge_test(oc::CLP& cmd);
int splitted_quick_sort_with_duplicate_elements_test(oc::CLP& cmd);
int arith_sort_test(oc::CLP& cmd);

// matrix multiplication tests
int fixed_matrix_mult_test(oc::CLP& cmd);
int splitted_fixed_matrix_mult_test(oc::CLP& cmd);

int correctness_cipher_index_pta(oc::CLP& cmd);
int correctness_sort_pta(oc::CLP& cmd);
int correctness_sum_pta(oc::CLP& cmd);
int correctness_max_pta(oc::CLP& cmd);
int correctness_metric_pta(oc::CLP& cmd);

bool check_result(const std::string& func_name, aby3::i64Matrix& test,
                  aby3::i64Matrix& res);

bool check_result(const std::string& func_name, std::vector<aby3::i64Matrix> test, std::vector<aby3::i64Matrix> res);

template <aby3::Decimal D>
bool check_result(const std::string& func_name, aby3::f64Matrix<D>& test,
                  aby3::f64Matrix<D>& res) {
    return check_result(func_name, test.i64Cast(), res.i64Cast());
}

bool check_result(const std::string& func_name, aby3::i64 test, aby3::i64 res);

template <typename T>
typename std::enable_if<std::is_pod<T>::value, bool>::type check_result(
    const std::string& func_name, std::vector<T>& test, std::vector<T>& res) {
    bool check_flag = true;
    for (size_t i = 0; i < test.size(); i++) {
        if (test[i] != res[i]) {
            check_flag = false;
        }
    }
    if (!check_flag) {
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
#ifdef SHOW_TEST_CASE
        debug_info("test case: ");
        debug_output_vector(test);
        debug_info("expected result: ");
        debug_output_vector(res);
#endif
    } else {
        debug_info("\033[32m" + func_name + " SUCCESS!" + "\033[0m\n");
    }
    return check_flag;
}

template <typename T>
typename std::enable_if<std::is_pod<T>::value, bool>::type check_result(
    const std::string& func_name, T& test, T& res) {
    bool check_flag = (test == res);
    if (!check_flag) {
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
#ifdef SHOW_TEST_CASE
        debug_info("test case: ");
        debug_output_value(test);
        debug_info("expected result: ");
        debug_output_value(res);
#endif
    } else{
        debug_info("\033[32m" + func_name + " SUCCESS!" + "\033[0m\n");
    }
    return check_flag;
}
