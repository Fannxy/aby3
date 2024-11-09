from task_execution import *
from system_monitor import *
from system_profile import *
import pandas as pd

tasks = ["cipher_index", "max", "sort", "sum"]
multi_step_tasks = ["sort"]
c_range = np.array([16])
b_start = 2**15
b_range = np.array([2**20])

# TEST#######
# c_range = np.array([256])
# b_start = 2**15
# b_range = np.array([2**20])
##########

target_tasks = {
    # "cipher_index": {"M_list": [1073741824], "N_list": [1]},
    "max": {"M_list": [16777216], "N_list": [1]},
    # "sort": {"M_list": [32768], "N_list": [32768]},
    # "sum": {"M_list": [1073741824], "N_list": [1]},
    # "metric": {"M_list": [1073741824], "N_list": [1]}
}

monitor = SystemMonitor(0.1)

def analyze_c_b():
    
    for task in tasks:
        df_dict = {
            "c": [],
            "b": [],
            "unit_time": []
        }
        task_log_folder = root_folder + "F-prof-Record/" + task + "/"
        if(not os.path.exists(task_log_folder)):
            os.makedirs(task_log_folder)
            
        for c in c_range:
            
            b = b_start
            
            Mem = get_free_memory()
            vector_limit = int(Mem / (3 * BITSIZE * c))
            vector_limit = 2**np.floor(np.log2(vector_limit))
            
            while b < vector_limit:
            
                logging_file = task_log_folder + f"monitor-{task}-c-{c}-b-{b}.log"
                
                monitor.start_all(interface=NETWORK_INTERFACE)
                unit_time = get_subtask_unit_time(task, c, b, task_log_folder)
                monitor.stop_and_output(logging_file)
                
                usage_dict = get_usage_dict(logging_file)
                draw_usage_graph(usage_dict, task_log_folder + f"{task}-c-{c}-b-{b}.png")
                
                df_dict["c"].append(c)
                df_dict["b"].append(b)
                df_dict["unit_time"].append(unit_time)
                df_dict["ratio"].append(unit_time / b)
        
                b *= 2
                
        df = pd.DataFrame(df_dict)
        df.to_excel(task_log_folder + f"{task}-prof.xlsx")
    
    return


def large_scale_analysis():
    
    for task in target_tasks.keys():
        
        df_dict = {
            "c": [],
            "B": [],
            "M": [],
            "N": [],
            "time_setup": [],
            "time_data_prepare": [],
            "combine": [],
            "subTask": []
        }
        
        task_log_folder = root_folder + "F-analys-Record/" + task + "/"
        if(not os.path.exists(task_log_folder)):
            os.makedirs(task_log_folder)
    
        
        for i in range(len(target_tasks[task]["M_list"])):
            M, N = target_tasks[task]["M_list"][i], target_tasks[task]["N_list"][i] 
            for c in c_range:
                Mem = get_free_memory()
                vector_limit = int(Mem / (3 * BITSIZE * c))
                vector_limit = 2**np.floor(np.log2(vector_limit))
                for b in b_range:
                    if(c * b > M * N):
                        continue
                    if(b > vector_limit):
                        continue
                    log_file = task_log_folder + f"exhv-log-config-N={N}-M={M}-c={c}-B={b}.log"
                    test_args = f" -{task} -N {N} -M {M} -c {c} -B {b} -logFile {log_file}"
                    monitor_logging_file = task_log_folder + f"monitor-exhv-log-config-N={N}-M={M}-c={c}-B={b}.log"
                    
                    monitor.start_all(interface=NETWORK_INTERFACE)
                    os.system(f"{execution_shell} {c} \"{test_args}\"")
                    monitor.stop_and_output(monitor_logging_file)
                    
                    usage_dict = get_usage_dict(monitor_logging_file)
                    draw_usage_graph(usage_dict, task_log_folder + f"{task}-c-{c}-b-{b}.png")
                    
                    with open(log_file, 'r') as f:
                        content = f.read()
                    
                    if(task in multi_step_tasks):
                        result = {
                            "time_setup": [],
                            "time_data_prepare": [],
                            "subTask": [],
                            "combine": [],
                            # "total": []
                        }
                        patterns = {
                            "time_setup": r'time_setup:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)',
                            "time_data_prepare": r'time_data_prepare:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)',
                            "subTask": r'subTask:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)',
                            "combine": r'combine:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)',
                            # "total": r'total:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)'
                        }
                        
                        for key, pattern in patterns.items():
                            matches = re.findall(pattern, content)
                            result[key].extend(float(match[0]) for match in matches)
                        
                        for key in result.keys():
                            df_dict[key].append(np.sum(result[key]))
                        
                    else:
                        performance_pattern = r'(time_setup|time_data_prepare|subTask|combine):\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)'
                        matches = re.findall(performance_pattern, content)
                        if(len(matches) != 4):
                            print(f"filename = {log_file} | matches = {matches}")
                        for match in matches:
                            df_dict[match[0]].append(float(match[1]))
                            
                    df_dict["c"].append(c)
                    df_dict["B"].append(b)
                    df_dict["M"].append(M)
                    df_dict["N"].append(N)
        
        df = pd.DataFrame(df_dict)
        df.to_excel(task_log_folder + f"{task}-analysis.xlsx")
        
    return
            
            
if __name__ == "__main__":
    
    os.system(f"python {root_folder}build.py")
    if(remote_compile):
        os.system(f"scp -r {root_folder}out {server1}:{root_folder}")
        os.system(f"scp -r {root_folder}out {server2}:{root_folder}")
    
    # analyze_c_b()
    large_scale_analysis()