"""System profiling tools of PtA.
"""
import numpy as np
import pandas as pd
import os
import psutil
import re
import json
import argparse
import time
from scipy.optimize import lsq_linear
# from task_execution import *
from system_monitor import SystemMonitor, get_usage_dict, draw_usage_graph

MAX_SAMPLE_NUM=10
BITSIZE = 64
GAP = 2**10 
MAXIMUM_VECTOR_SIZE = 2**25
MAXIMUM_CORES=16
remote_compile=True
server1="aby31"
server2="aby32"
START_B = 2**15
END_B = 2**23
NETWORK_INTERFACE = "ens110"
MAXMUM_MEMORY = 256*1024*1024*1024

root_folder = "/root/aby3/"
# root_folder = "/home/tsingj_ubuntu/fanxy/PtA/aby3/"
script_folder = root_folder + "Eval/basic/"
system_profiling_folder = root_folder + "system_profiling/"
system_profiling_log = system_profiling_folder + "system_profiling.log"
system_config_file = system_profiling_folder + "system.json"
network_config_file = system_profiling_folder + "network.json"
task_profiling_foler = root_folder + "task_profiling/"

def is_sorted(a):
    return np.all(np.diff(a) >= 0)

def get_linear_coefficients(X, y):
    """Get the coefficients a, b of tt(c) = a*c + b.
    """
    if type(X) is not np.ndarray:
        X = np.array(X)
        y = np.array(y)
    
    if (not is_sorted(X)):
        sorted_args = np.argsort(X)
        X = X[sorted_args]
        y = y[sorted_args]
    X = np.vstack([X, np.ones(len(X))]).T
    
    # a, b = np.polyfit(X, y, deg=1, w=1/X)
    res = lsq_linear(X, y, bounds=(0, np.inf), method='trf', lsmr_tol='auto', verbose=0)
    a, b = res.x
    return a, b


def get_cpu_cores():
    """Get the number of CPU cores.
    """
    cpu_cores = os.cpu_count()
    if(cpu_cores > MAXIMUM_CORES):
        cpu_cores = MAXIMUM_CORES
    return cpu_cores


def get_free_memory():
    free_mem = psutil.virtual_memory().available
    if(free_mem > MAXMUM_MEMORY):
        free_mem = MAXMUM_MEMORY
    return free_mem


def system_profiling(sampling_ratio=0.1):
    """Get the system profiling information.
    """
    
    profiler = SystemMonitor(0.01)
    
    if(not os.path.exists(system_profiling_folder)):
        os.mkdir(system_profiling_folder)
    
    cpu_cores = get_cpu_cores()
    
    sample_num = int(cpu_cores * sampling_ratio)
    
    if(sample_num > MAX_SAMPLE_NUM):
        sample_num = MAX_SAMPLE_NUM
    
    test_tasks_list = 2**(np.linspace(1, np.log2(cpu_cores), sample_num)).astype("int")
    print(test_tasks_list)
    
    profiler_script = script_folder + "dis_mpi_exec.sh"
    profile_args=" -system_profile"
    log_folder_args=" -logFile " + system_profiling_log
    
    start_b = END_B
    end_b = END_B
    
    vector_args = f" -startB {start_b} -endB {end_b}"
    
    test_args = profile_args + log_folder_args + vector_args
    
    # run the profile.
    os.system(f"cd {root_folder}; python build.py") # compile the project
    if(remote_compile):
        os.system(f"scp -r {root_folder}out {server1}:{root_folder}")
        os.system(f"scp -r {root_folder}out {server2}:{root_folder}")
    
    print("System profiling...")
    start = time.time()
    for task in test_tasks_list:
        logging_file = system_profiling_folder + f"profile_task_{task}.log"
        profiler.start_cpu()
        profiler.start_memory()
        profiler.start_network(interface = NETWORK_INTERFACE)
        
        os.system(f"{profiler_script} {task} \"{test_args}\"")
        
        profiler.stop_cpu()
        profiler.stop_memory()
        profiler.stop_network()
        
        profiler.output(logging_file)
    end = time.time()
    system_profiling_time = (end - start)
    with open(f"{system_profiling_folder}cost.txt", "a") as f:
        print(f"System profiling done, using {system_profiling_time} sec", file=f)
    print("System profiling done.")
    
    return


def system_profiling2(sampling_ratio=0.1):
    if(not os.path.exists(system_profiling_folder)):
        os.mkdir(system_profiling_folder)
    
    cpu_cores = get_cpu_cores()
    
    sample_num = int(cpu_cores * sampling_ratio)
    
    if(sample_num > MAX_SAMPLE_NUM):
        sample_num = MAX_SAMPLE_NUM
    
    test_tasks_list = 2**(np.linspace(1, np.log2(cpu_cores), sample_num)).astype("int")
    print(test_tasks_list)
    
    profiler_script = script_folder + "dis_mpi_exec.sh"
    profile_args=" -system_profile"
    log_folder_args=" -logFile " + system_profiling_log
    
    test_args = profile_args + log_folder_args
    for task in test_tasks_list:
        os.system(f"{profiler_script} {task} \"{test_args}\"")
    
    return
        

def get_system_configs():
    
    task_num_list = []
    time_setup_list = []
    vector_size_list = []
    communication_time_list = []

    with open(system_profiling_log, 'r') as file:
        for line in file:
            task_num_match = re.search(r'tasknum: (\d+)', line)
            time_setup_match = re.search(r'time_setup: (\d+(\.\d+)?) milliseconds', line)
            communication_match = re.search(r'time_communication\.(\d+): (\d+(\.\d+)?) milliseconds', line)
            if task_num_match:
                task_num_list.append(int(task_num_match.group(1)))
            if time_setup_match:
                time_setup_list.append(float(time_setup_match.group(1)))
            if communication_match:
                vector_size_list.append(int(communication_match.group(1)))
                communication_time_list.append(float(communication_match.group(2)))
    
    bandwidth_utilization_list = ((np.array(vector_size_list) * 8 / 1000) * np.array(task_num_list) / np.array(communication_time_list))
    c_effect = task_num_list[np.argmax(bandwidth_utilization_list)]
    a, b = get_linear_coefficients(task_num_list, time_setup_list)
                
    return a, b, c_effect

def get_system_configs2():
    task_num_list = []
    time_setup_list = []
    
    with open(system_profiling_log, 'r') as file:
        for line in file:
            task_num_match = re.search(r'task num: (\d+)', line)
            time_setup_match = re.search(r'time_setup: (\d+(\.\d+)?) milliseconds', line)
            if task_num_match:
                task_num_list.append(int(task_num_match.group(1)))
            if time_setup_match:
                time_setup_list.append(float(time_setup_match.group(1)))
    
    print(task_num_list)
    print(time_setup_list)
    
    a, b = get_linear_coefficients(task_num_list, time_setup_list)
    return a, b, -1


def get_network_configs():
    """Get the network profiling information.
    """
    
    task_dict = {}
    
    task_num = -1
    
    with open(system_profiling_log, 'r') as file:
        lines = file.readlines()
        for line in lines:
            task_num_match = re.search(r'tasknum: (\d+)', line)
            if task_num_match:
                # task_num_list.append(int(task_num_match.group(1)))
                task_num = int(task_num_match.group(1))
                task_dict[task_num] = {
                    "latency": -1,
                    "vector_size_list": [],
                    "communication_time_list": [],
                    "bandwidth_list": [],
                }
            
            latency_match = re.search(r'time_latency: (\d+(\.\d+)?) milliseconds', line)
            if latency_match:
                task_dict[task_num]["latency"] = float(latency_match.group(1))
            
            communication_match = re.search(r'time_communication\.(\d+): (\d+(\.\d+)?) milliseconds', line)
            if communication_match:
                task_dict[task_num]["vector_size_list"].append(int(communication_match.group(1)))
                task_dict[task_num]["communication_time_list"].append(float(communication_match.group(2)))
    
    # post-process the task profiling dict. 
    for task_num in task_dict.keys():
        sort_vector_size_args = np.argsort(task_dict[task_num]["vector_size_list"])
        task_dict[task_num]["vector_size_list"] = np.array(task_dict[task_num]["vector_size_list"])[sort_vector_size_args]
        task_dict[task_num]["communication_time_list"] = np.array(task_dict[task_num]["communication_time_list"])[sort_vector_size_args]
        task_dict[task_num]["bandwidth_list"] = ((task_dict[task_num]["vector_size_list"] * 8 / 1000)  / task_dict[task_num]["communication_time_list"]) * task_num
        
        task_dict[task_num]["vector_size_list"] = task_dict[task_num]["vector_size_list"].tolist()
        task_dict[task_num]["communication_time_list"] = task_dict[task_num]["communication_time_list"].tolist()
        task_dict[task_num]["bandwidth_list"] = task_dict[task_num]["bandwidth_list"].tolist()
    
    
    with open(network_config_file, "w") as f:
        json.dump(task_dict, f)
        
    return task_dict
    

def get_deployment_configs(save_file=system_config_file):
    
    C = get_cpu_cores()
    M = get_free_memory()  
    
    p_start = time.time()
    if(not os.path.exists(system_profiling_log)):
        system_profiling2()
        
    a, b, _ = get_system_configs2()
    p_end = time.time()
    
    system_configs = {'C': C, 'M': M, 'a': a, 'b': b, 'pcost': (p_end - p_start)}
    
    with open(save_file, 'w') as file:
        json.dump(system_configs, file)
    
    return C, M, a, b


def get_c_effect():
    with open(system_config_file, 'r') as f:
        system_configs = json.load(f)
        c_effect = system_configs["c_effect"]
    return c_effect


def get_subtask_unit_time(task, c, B, spec_task_profiling_folder):
    """Get the unit time of the subtask.
    """
    profile_args=" -task_profile -getUnit"
    task_exec_args = f" -task {task} -startB {B} -logFolder {spec_task_profiling_folder} -endingB {B} -gap {GAP}"
    test_args = profile_args + task_exec_args
    
    if(not os.path.exists(spec_task_profiling_folder + "probing.res")):
        os.system(f"{script_folder}dis_mpi_exec.sh {c} \"{test_args}\"")
    
    target_file = spec_task_profiling_folder + "probing.res"
    if(not os.path.exists(target_file)):
        print(f"file {target_file} not found.")
        return -1
        
    # get the profiling results.
    with open(target_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            unit_time_match = re.search(fr'unit_time-{B}-{c}: (\d+(\.\d+)?) milliseconds', line)
            if(unit_time_match):
                unit_time = float(unit_time_match.group(1))
            else:
                unit_time = -1
    return unit_time


def get_task_maximum_c(task_dict, spec_task_profiling_folder):
    optimal_B = task_dict["optimal_B"]
    c_effect = get_c_effect()
    unit_time = task_dict["unit_time"]
    c_ratio = unit_time / c_effect
    C_MAX = get_cpu_cores()
    c = c_effect * 2
    logging_file = spec_task_profiling_folder + "probing.log"
    
    while(c < C_MAX):
        maximum_vector = int(get_free_memory() / (3 * BITSIZE * c))
        if(maximum_vector > optimal_B):
            unit_time = get_subtask_unit_time(task, c, optimal_B, spec_task_profiling_folder)
        else:
            unit_time = get_subtask_unit_time(task, c, maximum_vector, spec_task_profiling_folder)
        ratio = unit_time / c 
        
        with open(logging_file, 'a') as f:
            print(f"c = {c} | ratio = {ratio} | unit_time = {unit_time}", file=f)
            
        if(ratio <= (c_ratio * 0.9)):
            c *= 2
            c_ratio = ratio
        else:
            break

    return c

# task profiler. 
def get_task_optimal_B(task, spec_task_profiling_folder):
    c_effect = get_c_effect()
    M = get_free_memory()
    
    vector_limit = int(M / (3 * BITSIZE * c_effect))
    vector_limit = 2**np.floor(np.log2(vector_limit))
    print("vector_limit: ", vector_limit)
    
    if(vector_limit > MAXIMUM_VECTOR_SIZE):
        vector_limit = MAXIMUM_VECTOR_SIZE
    vector_start = int(2**10)

    profile_args=" -task_profile"
    task_exec_args = f" -task {task} -logFolder {spec_task_profiling_folder} -startB {vector_start} -endingB {vector_limit} -gap {GAP}"
    test_args = profile_args + task_exec_args
    
    target_file = spec_task_profiling_folder + "probing.res"
    if(not os.path.exists(target_file)):
        # run the profile.
        os.system(f"{script_folder}dis_mpi_exec.sh {c_effect} \"{test_args}\"")
    
    # get the profiling results.
    task_config_dict = {}
    with open(target_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            optimal_B_match = re.search(r'optimal_B: (\d+)', line)
            unit_time_match = re.search(r'unit_time: (\d+\.\d+) milliseconds', line)
            ratio_match = re.search(r'ratio: (\d+\.\d+)', line)
            if(optimal_B_match):
                task_config_dict["optimal_B"] = int(optimal_B_match.group(1))
            if(unit_time_match):
                task_config_dict["unit_time"] = float(unit_time_match.group(1))
            if(ratio_match):
                task_config_dict["ratio"] = float(ratio_match.group(1))
    
    # save logs into a excel. 
    logging_file = spec_task_profiling_folder + "probing.log-0-0"
    task_dict = {
        "b": [],
        "tsb": [],
        "ratio": []
    }
    with open(logging_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            b_match = re.search(r'^b: (\d+)', line)
            tsb_match = re.search(r'^time_c: (\d+\.\d+)', line)
            ratio_match = re.search(r'^ratio: (\d+\.\d+)', line)
            if(b_match):
                task_dict["b"].append(int(b_match.group(1)))
            if(tsb_match):
                task_dict["tsb"].append(float(tsb_match.group(1)))
            if(ratio_match):
                task_dict["ratio"].append(float(ratio_match.group(1)) * 10**4)
        print(task_dict)
    df = pd.DataFrame(task_dict)
    df.to_excel(spec_task_profiling_folder + "probing.xlsx")
        
    return task_config_dict


def get_task_optimal_B(task, spec_task_profiling_folder, c, b_start, b_end):
    profile_args = " -task_profile"
    task_exec_args = f" -task {task} -logFolder {spec_task_profiling_folder} -startB {b_start} -endingB {b_end} -gap {GAP}"
    print(f"shell: {task_exec_args}")
    
    target_file = spec_task_profiling_folder + "probing.res"
    # run the profile.
    os.system(f"{script_folder}dis_mpi_exec.sh {c} \"{profile_args} {task_exec_args}\"")
    
    # get the profiling results.
    task_config_dict = {}
    with open(target_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            optimal_B_match = re.search(r'optimal_B: (\d+)', line)
            unit_time_match = re.search(r'unit_time: (\d+\.\d+) milliseconds', line)
            ratio_match = re.search(r'ratio: (\d+\.\d+)', line)
            if(optimal_B_match):
                task_config_dict["optimal_B"] = int(optimal_B_match.group(1))
            if(unit_time_match):
                task_config_dict["unit_time"] = float(unit_time_match.group(1))
            if(ratio_match):
                task_config_dict["ratio"] = float(ratio_match.group(1))
    
    return task_config_dict


def get_task_configs(task, spec_task_profiling_folder):

    save_file = task_profiling_foler + task + ".json"
    if(os.path.exists(save_file)):
        with open(save_file, 'r') as file:
            task_config_dict = json.load(file)
            return task, task_config_dict
    
    start = time.time()
    task_config_dict = get_task_optimal_B(task, spec_task_profiling_folder)
    print(task_config_dict)
    
    # then we profile for the c_max. 
    c_max = get_task_maximum_c(task_config_dict, spec_task_profiling_folder)
    task_config_dict["c_max"] = c_max
    end = time.time()
    task_config_dict["pcost"] = end - start

    with open(save_file, 'w') as file:
        json.dump(task_config_dict, file)
    
    return task, task_config_dict


def is_bandwidth_limit(usage_dict, success_ratio, bandwidth_threshold_upper, bandwidth_threshold_lower, lower_ratio, cpu_threadshold):
    """bandwidth_threashold = bandwidth * 90%, which is a magic number. 
    """
    recv_list = usage_dict["network_recv"]
    send_list = usage_dict["network_send"]
    
    bandwidth_list = np.max(np.array([np.array(recv_list),np.array(send_list)]), axis=0)
    cpu_list = np.array(usage_dict["cpu"]) / get_cpu_cores()
    
    len_threshold = int(len(bandwidth_list) * success_ratio)
    len_lower_threshold = int(len_threshold * lower_ratio)
    
    sum_bandwidth = np.sum(bandwidth_list[:len_threshold])
    sum_cpu = np.sum(cpu_list[:len_threshold])
    count_bandwidth = np.sum(bandwidth_list[:len_threshold] > bandwidth_threshold_lower)
    
    
    for i in range(len_threshold, min(len(bandwidth_list), len(cpu_list))):
        avg_bandwidth = sum_bandwidth / len_threshold
        avg_cpu = sum_cpu / len_threshold
        
        if((avg_bandwidth > bandwidth_threshold_upper) and (count_bandwidth > len_lower_threshold) and (avg_cpu < cpu_threadshold)):
            return True
        sum_bandwidth = sum_bandwidth - bandwidth_list[i - len_threshold] + bandwidth_list[i]
        sum_cpu = sum_cpu - cpu_list[i - len_threshold] + cpu_list[i]
        count_bandwidth += 1 if bandwidth_list[i] > bandwidth_threshold_lower else 0
        count_bandwidth -= 1 if bandwidth_list[i - len_threshold] > bandwidth_threshold_lower else 0
    
    return False


def task_profile(task, spec_task_profiling_folder, bandwidth=25):
    """Profile the task using the new profilng mechanism!
    """
    C = get_cpu_cores()
    Mem = get_free_memory()
    vector_limit = int(Mem / (3 * BITSIZE * C))
    if(task == "metric"):
        vector_limit = int(Mem*0.8 / (3 * BITSIZE * C))
        
    vector_limit = 2**np.floor(np.log2(vector_limit))
    
    # initial profiling. 
    vector_start = vector_limit

    monitor_logging_file = spec_task_profiling_folder + "monitor-probing.log"    
    monitor = SystemMonitor(0.1)
    
    p_start = time.time()
    if(not os.path.exists(monitor_logging_file)):
        monitor.start_all(interface=NETWORK_INTERFACE)
        unit_time = get_subtask_unit_time(task, int(C), int(vector_start), spec_task_profiling_folder)
        monitor.stop_and_output(monitor_logging_file)
    else:
        unit_time = get_subtask_unit_time(task, int(C), int(vector_start), spec_task_profiling_folder)
    
    usage_dict = get_usage_dict(monitor_logging_file)
    # define whether this task is bandwidth-bounded. 
    bandwidth_upper = bandwidth * 0.9
    bandwidth_lower = bandwidth * 0.8
    lower_ratio = 0.9
    success_ratio = 0.1 # 10% of the time is successively higher.
    cpu_threadshold = 0.5 # the per-core cpu usage is lower than cpu_threshold.
    
    if(is_bandwidth_limit(usage_dict, success_ratio, bandwidth_upper, bandwidth_lower, lower_ratio, cpu_threadshold)):
        # bandwidth bounded.
        p_end = time.time()
        task_config_dict = {
            "optimal_B": vector_start,
            "unit_time": unit_time,
            "pcost": (p_end - p_start),
        }
        draw_usage_graph(usage_dict, spec_task_profiling_folder + "usage.png")
        with open(spec_task_profiling_folder + f"config.json", "w") as f:
            json.dump(task_config_dict, f)
        return task_config_dict
    
    small_vector = 1024
    task_config_dict = get_task_optimal_B(task, spec_task_profiling_folder, C, small_vector, vector_limit/2)
    p_end = time.time()
    
    draw_usage_graph(usage_dict, spec_task_profiling_folder + "usage.png")
    task_config_dict["pcost"] = (p_end - p_start)
    with open(spec_task_profiling_folder + "config.json", "w") as f:
        json.dump(task_config_dict, f)
    
    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--system_profile', action='store_true', help='profile system or not')
    parser.add_argument('--task_profile', type=str, help='which task you`d like to profile.')
    args = parser.parse_args()
    
    os.system(f"cd {root_folder}; python build.py") # compile the project
    if(remote_compile):
        os.system(f"scp -r {root_folder}out {server1}:{root_folder}")
        os.system(f"scp -r {root_folder}out {server2}:{root_folder}")

    if args.system_profile:
        get_deployment_configs()

    if(args.task_profile is not None):
        task = args.task_profile
        spec_task_profiling_folder = task_profiling_foler + task + "/"

        if(not os.path.exists(task_profiling_foler)):
            os.mkdir(task_profiling_foler)

        if(not os.path.exists(spec_task_profiling_folder)):
            os.mkdir(spec_task_profiling_folder)

        task_profile(task, spec_task_profiling_folder)