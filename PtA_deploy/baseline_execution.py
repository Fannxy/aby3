from get_execution_parameters import *
from system_profile import *
from system_monitor import SystemMonitor

target_size = [1024, 4096, 16384, 65536, 262144, 1048576]
# target_size = [1073741824]
execution_shell = script_folder + "dis_mpi_exec.sh"

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--sqrtOram', action='store_true', help='run sqrt oram baseline')
    args = parser.parse_args()
    
    monitor = SystemMonitor(0.1)
    
    if(args.sqrtOram):
    
        # prepare the project.
        os.system(f"cp {root_folder}frontend/main.soram {root_folder}frontend/main.cpp")
        os.system(f"python {root_folder}build.py")
        if(remote_compile):
            os.system(f"scp -r {root_folder}out {server1}:{root_folder}")
            os.system(f"scp -r {root_folder}out {server2}:{root_folder}")
        
        log_folder = root_folder + "Record_sqrtOram/"
        if(not os.path.exists(log_folder)):
            os.makedirs(log_folder)
        
        
        # run the tasks.
        task_args = " -BenchmarkORAM"
        access_times = [int(np.sqrt(size)) for size in target_size]
        for i in range(len(target_size)):
            M = target_size[i]
            N = access_times[i]
            
            params = f" -M {M} -N {N} "
            logFile = log_folder + f"sqrtOram-log-M={M}-N={N}.log"
            monitor_file = log_folder + f"monitor-sqrtOram-log-M={M}-N={N}.log"
            test_args = task_args + params + f" -logFile {logFile}"
            
            # print(f"{execution_shell} 1 \"{test_args}\"")
            monitor.start_all(interface=NETWORK_INTERFACE)
            os.system(f"{execution_shell} 1 \"{test_args}\"")
            monitor.stop_and_output(monitor_file)
        
            
            
            
            
        
        