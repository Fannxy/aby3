#pragma once
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>
#include "PtATasks.h"
#include "PtATests.h"

#define PROFILER_PREPARE \
    std::string logging_folder, probing_log, probing_res; \
    get_value("logFolder", cmd, logging_folder); \
    probing_log = logging_folder + "probing.log-" + std::to_string(role) + "-" + std::to_string(rank); \
    probing_res = logging_folder + "probing.res"; \
    bool get_unit = false; \
    if(cmd.isSet("getUnit")){ \
        get_unit = true; \
    } \
    size_t start_b, gap, ending_b; \
    get_value("startB", cmd, start_b); \
    if(get_unit){ \
        gap = -1; \
        ending_b = -1; \
    } \
    else{ \
        get_value("gap", cmd, gap); \
        get_value("endingB", cmd, ending_b); \
    } \
    size_t n=1; \
    double last_ratio = -1; \
    double ratio = -1; \
    double time_c = -1; \


#define PROFILER_RECORDER \
    if(get_unit){ \
        double unit_time = get_unit_time(start_b, evaluate_task, probing_log); \
        if(rank == 0){ \
            std::ofstream stream(probing_res, std::ios::app); \
            stream << "unit_time-" << start_b << "-" << size << ": " << unit_time << " milliseconds" << std::endl; \
            stream.close(); \
        } \
    } \
    else{ \
        size_t optimal_vector_size; \
        std::tie(optimal_vector_size, time_c) = get_optimal_vector_size(start_b, ending_b, gap, evaluate_task, probing_log); \
        if(rank == 0){ \
            std::ofstream stream(probing_res, std::ios::app); \
            stream << "optimal_B: " << optimal_vector_size << std::endl; \
            stream << "unit_time: " << time_c << " milliseconds" << std::endl; \
            stream << "ratio: " << time_c / optimal_vector_size << std::endl; \
            stream.close(); \
        } \
    } \


int pta_system_profile(oc::CLP& cmd);

int communication_profile(oc::CLP& cmd);

int task_profile(oc::CLP& cmd);

std::pair<size_t, double> get_optimal_vector_size(size_t b_start, size_t b_end, size_t gap, std::function<std::pair<double, double>(size_t)> evaluate_task, std::string logging_file);

double get_unit_time(size_t b, std::function<std::pair<double, double>(size_t)> evaluate_task, std::string logging_file);

int cipher_index_profile(oc::CLP& cmd);

int max_profile(oc::CLP& cmd);

int sort_profile(oc::CLP& cmd);

int sum_profile(oc::CLP& cmd);

int metric_profile(oc::CLP& cmd);