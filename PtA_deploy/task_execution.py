from get_execution_parameters import *
from system_profile import *
from system_monitor import SystemMonitor

target_tasks = {
    # "cipher_index": {"M_list": [1048576], "N_list": [1024]},
    "cipher_index": {"M_list": [16777216], "N_list": [1]},
    "max": {"M_list": [16777216], "N_list": [1]},
    # "sort": {"M_list": [32768], "N_list": [32768]},
    "sum": {"M_list": [16777216], "N_list": [1]},
    "metric": {"M_list": [16777216], "N_list": [1]},
}

best_configurations = {
    "cipher_index": {
        "M_list": [1048576, 16777216, 268435456, 1073741824],
        "N_list": [1, 1, 1, 1],
        "c_list": [16, 32, 64, 256],
        "B_list": [1048576, 4194304, 4194304, 4194304]
    },
    "sum": {
        "M_list": [1048576, 16777216, 268435456, 1073741824],
        "N_list": [1, 1, 1, 1],
        "c_list": [16, 32, 128, 256],
        "B_list": [1024, 65536, 32768, 2048]
    },
    "sort": {
        "M_list": [1024, 2048, 16384, 32768],
        "N_list": [1024, 2048, 16384, 32768],
        "c_list": [16, 32, 128, 32],
        "B_list": [4194304, 1048576, 4194304, 262144]
    },
    "max": {
        "M_list": [1048576, 16777216, 268435456, 1073741824],
        "N_list": [1, 1, 1, 1],
        "c_list": [16, 32, 32, 32],
        "B_list": [1048576, 1048576, 1048576, 1048576]
    },
    "metric": {
        "M_list": [1048576, 16777216, 268435456, 1073741824],
        "N_list": [1, 1, 1, 1],
        "c_list": [16, 32, 32, 32],
        "B_list": [1048576, 262144, 4194304, 1048576]    
    }
}

# target_tasks = {
#     "cipher_index": {"M_list": [1073741824], "N_list": [1]},
# }

# target_tasks = {
#     "cipher_index": {"M_list": [1048576, ], "N_list": [1, 1, 1, 1, 1, 1]},
# }

# c_range = 2**np.linspace(1, 8, 8).astype("int")
# c_range = np.array([2**0, 2**1, 2**2, 2**3, 2**4, 2**5, 2**6, 2**7, 2**8])
# c_range = np.array([8, 16, 32, 64])
c_range = np.array([2**0, 2**1, 2**2, 2**3])
vec_range = np.array([2**10, 2**12, 2**14, 2**16, 2**18, 2**20, 2**22])
# 16 - 1048576; 256 - 324288; 32 - 65536; 

execution_shell = script_folder + "dis_mpi_exec.sh"

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--pta', action='store_true', help='run pta task')
    parser.add_argument('--exhv', action='store_true', help='run exhv task')
    parser.add_argument('--config', action='store_true', help='run task with configurations')
    args = parser.parse_args()
    print("config - ", args.config)
    PTA_TEST = args.pta
    
    monitor = SystemMonitor(0.1)
    
    # prepare the project.
    os.system(f"cp {root_folder}frontend/main.pta {root_folder}frontend/main.cpp")
    os.system(f"python {root_folder}build.py")
    if(remote_compile):
        os.system(f"scp -r {root_folder}out {server1}:{root_folder}")
        os.system(f"scp -r {root_folder}out {server2}:{root_folder}")
    
    # run the tasks. 
    # total_tasks = len(target_tasks.keys()) * len(target_tasks["cipher_index"]["M_list"])
    counter = 0
    
    if(args.config):
        for task in best_configurations:
            task_log_folder = root_folder + "Record/" + task + "/"
            if(not os.path.exists(task_log_folder)):
                os.makedirs(task_log_folder)
            M_list = best_configurations[task]["M_list"]
            N_list = best_configurations[task]["N_list"]
            c_list = best_configurations[task]["c_list"]
            B_list = best_configurations[task]["B_list"]
            
            for i in range(len(M_list)):
                M, N, c, B = M_list[i], N_list[i], c_list[i], B_list[i]
                log_file = task_log_folder + f"config-log-config-N={N}-M={M}-c={c}-B={B}.log"
                test_args = f" -{task} -N {N} -M {M} -c {c} -B {B} -logFile {log_file}"
                monitor_logging_file = task_log_folder + f"monitor-config-log-config-N={N}-M={M}-c={c}-B={B}.log"
                
                monitor.start_all(interface=NETWORK_INTERFACE)
                print(f"{execution_shell} {c} \"{test_args}\"")
                os.system(f"{execution_shell} {c} \"{test_args}\"")
                monitor.stop_and_output(monitor_logging_file)
                
        exit(0)
    
    
    for task in target_tasks.keys():
        
        # prepare the log folder.
        task_log_folder = root_folder + "Record/" + task + "/"
        if(not os.path.exists(task_log_folder)):
            os.makedirs(task_log_folder)
        
        M_list = target_tasks[task]["M_list"]
        N_list = target_tasks[task]["N_list"]
        
        # generate the execution scripts and run the tests.
        for i in range(len(M_list)):
            M, N = M_list[i], N_list[i]
            
            if(args.pta):
                c, B = get_execution_parameters(task, M, N)
                log_file = task_log_folder + f"pta-log-config-N={N}-M={M}-c={c}-B={B}.log"
                test_args = f" -{task} -N {N} -M {M} -c {c} -B {B} -logFile {log_file}"
                monitor_logging_file = task_log_folder + f"monitor-pta-log-config-N={N}-M={M}-c={c}-B={B}.log"

                monitor.start_all(interface=NETWORK_INTERFACE)
                print(f"{execution_shell} {c} \"{test_args}\"")
                os.system(f"{execution_shell} {c} \"{test_args}\"")
                monitor.stop_and_output(monitor_logging_file)
                
            elif(args.exhv): # run the exhavsive search tests. 
                for c in c_range:
                    for B in vec_range:
                        log_file = task_log_folder + f"exhv-log-config-N={N}-M={M}-c={c}-B={B}.log"
                        test_args = f" -{task} -N {N} -M {M} -c {c} -B {B} -logFile {log_file}"
                        monitor_logging_file = task_log_folder + f"monitor-exhv-log-config-N={N}-M={M}-c={c}-B={B}.log"
                        
                        rounds_per_task = N*M / (c*B)
                        
                        if(rounds_per_task < 1):
                            continue
                        
                        counter += 1
                        with open('./progress.txt', 'a') as f:
                            print(f"c={c} | B = {B} | rounds_per_task = {rounds_per_task}", file=f)
                        
                        monitor.start_all(interface=NETWORK_INTERFACE)
                        print(f"{execution_shell} {c} \"{test_args}\"")
                        os.system(f"{execution_shell} {c} \"{test_args}\"")
                        monitor.stop_and_output(monitor_logging_file)
                    