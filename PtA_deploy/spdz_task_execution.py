from task_execution import target_tasks 
from system_monitor import SystemMonitor
import argparse
import os
import subprocess 
import threading

NETWORK_INTERFACE = "ens5"
LIMIT_TIME = 3

def watcher(p, remote_server1, remote_server2):
    try:
        p.wait(timeout=LIMIT_TIME * 60 * 60)
    except subprocess.TimeoutExpired:
        p.kill()
        os.system("pkill -f replicated-ring-party.x")
        os.system(f"ssh {remote_server1} 'pkill -f replicated-ring-party.x'")
        os.system(f"ssh {remote_server2} 'pkill -f replicated-ring-party.x'")

# root_folder="/home/tsingj_ubuntu/fanxy/MP-SPDZ/"

spdz_config_dict = {
    "cipher_index": {
        "M_list": [1048576],
        "N_list": [1],
        "ringsize": 256
    },
    "max": {
        "M_list": [1048576, 16777216, 268435456, 1073741824],
        "N_list": [1, 1, 1, 1],
        "ringsize": 64
    },
    "average": {
        "M_list": [1048576, 16777216, 268435456, 1073741824],
        "N_list": [1, 1, 1, 1],
        "ringsize": 64
    },
    "sort": {
        "M_list": [1024, 4096, 16384, 32768],
        "N_list": [1, 1, 1, 1],
        "ringsize": 64
    },
    "metric": {
        "M_list": [1048576, 16777216, 268435456],
        "N_list": [1, 1, 1, 1],
        "ringsize": 64
    },
}

tgl_config_dict = {
    "cipher_index":{
        "M_list": [1048576, 16777216],
        "N_list": [1, 1],
        "c_list": [64, 128, 256]
    },
    "max":{
        "M_list": [1048576, 16777216],
        "N_list": [1, 1],
        "c_list": [64, 128, 256]
    },
    "average":{
        "M_list": [1048576, 16777216],
        "N_list": [1, 1],
        "c_list": [64, 128, 256]
    },
    "sort":{
        "M_list": [1024, 4096],
        "N_list": [1, 1, 1],
        "c_list": [64, 128, 256]
    },
    "metric":{
        "M_list": [1048576, 16777216],
        "N_list": [1, 1],
        "c_list": [64, 128, 256]
    }
}

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--tgl', action='store_true', help='run pta task')
    parser.add_argument('--spdz', action='store_true', help='run spdz origional task')
    parser.add_argument('--compile_tgl', action='store_true', help='compile the spdz origional task')
    parser.add_argument('--compile_spdz', action='store_true', help='compile the spdz origional task')
    parser.add_argument('--root_folder', type=str, help="define the root folder")
    args = parser.parse_args()
    
    monitor = SystemMonitor(0.1)
    
    # root_folder="/root/MP-SPDZ/"
    # remote_root_folder="/root/MP-SPDZ/"
    root_folder="/home/ubuntu/configuration/MP-SPDZ/"
    remote_root_folder="/home/ubuntu/configuration/MP-SPDZ/"
    
    # parse the parameters.
    if args.root_folder is not None:
        root_folder = args.root_folder

    script_folder=root_folder+"Eval/"
    schedule_folder=root_folder+"Programs/Schedules/"
    record_folder=root_folder+"Record/"
    main_server = "spdz0"
    remote_server1 = "sosp1"
    remote_server2 = "sosp2"
    execution_shell = script_folder + "basic/dis_exec.sh"
    
    # prepare the project.
    if args.compile_tgl:
        os.system(f"{script_folder}tgl_benchmark_list.sh False True")
        # scp the compiled files to the remote servers.
        os.system(f"scp -r {root_folder}Programs {main_server}:{remote_root_folder}")
    
    if args.compile_spdz:
        os.system(f"{script_folder}spdz_benchmark_list.sh False True")
        # scp the compiled files to the remote servers.
        os.system(f"scp -r {root_folder}Programs {main_server}:{remote_root_folder}")
    
    # run the tasks.
    # note that you must be in the remote server. 
    protocol="replicated-ring-party.x"
    
    if args.tgl:
        source_file="pta_baseline"
        
        tgl_record_folder = record_folder + "tgl/"
        if(not os.path.exists(tgl_record_folder)):
            os.makedirs(tgl_record_folder)
            
        # # synchronize with the other servers. 
        # os.system(f"scp -r {root_folder}Programs {remote_server1}:{remote_root_folder}")
        # os.system(f"scp -r {root_folder}Programs {remote_server2}:{remote_root_folder}")
        
        for task_name in tgl_config_dict.keys():
            
            task_record_folder=tgl_record_folder+task_name+"/"
            if(not os.path.exists(task_record_folder)):
                os.makedirs(task_record_folder)
            
            for M, N in zip(tgl_config_dict[task_name]["M_list"], tgl_config_dict[task_name]["N_list"]):
                for c in tgl_config_dict[task_name]["c_list"]:
                    elog = f"{task_record_folder}tgl-{task_name}-{N}-{M}-1-{c}.log"
                    target_execution_args = f"{source_file}-{task_name}-{N}-{M}-1-{c} {protocol} {task_record_folder} {elog}"
                    
                    monitor_log = f"{task_record_folder}monitor-tgl-{task_name}-{N}-{M}-1-{c}.log"
                    
                    monitor.start_all(interface=NETWORK_INTERFACE)
                    # os.system(f"{execution_shell} {target_execution_args}")
                    
                    p = subprocess.Popen(f"{execution_shell} {target_execution_args}", shell=True)
                    t = threading.Thread(target=watcher, args=(p, remote_server1, remote_server2))
                    t.start()
                    t.join()
                    
                    monitor.stop_and_output(monitor_log)

    
    if args.spdz:
        source_file="benchmark_spdz_origional"
        
        spdz_record_folder = record_folder + "spdz/"
        if(not os.path.exists(spdz_record_folder)):
            os.makedirs(spdz_record_folder)
        
        # # synchronize with the other servers. 
        # os.system(f"scp -r {root_folder}Programs {remote_server1}:{remote_root_folder}")
        # os.system(f"scp -r {root_folder}Programs {remote_server2}:{remote_root_folder}")
        # os.system(f"scp -r {root_folder}Eval {remote_server1}:{remote_root_folder}")
        # os.system(f"scp -r {root_folder}Eval {remote_server2}:{remote_root_folder}")
        # os.system(f"scp {root_folder}HOST {remote_server1}:{remote_root_folder}")
        # os.system(f"scp {root_folder}HOST {remote_server2}:{remote_root_folder}")
        
        for task_name in spdz_config_dict.keys():
            
            task_record_folder=spdz_record_folder+task_name+"/"
            if(not os.path.exists(task_record_folder)):
                os.makedirs(task_record_folder)
            
            # recompile the virtual machine. 
            os.system(f"echo MOD = -DRING_SIZE={spdz_config_dict[task_name]['ringsize']} > {root_folder}CONFIG.mine")
            os.system(f"scp -r {root_folder}CONFIG.mine {remote_server2}:{remote_root_folder}")
            os.system(f"scp -r {root_folder}CONFIG.mine {remote_server1}:{remote_root_folder}")
            
            p1 = subprocess.Popen(f"cd {root_folder}; make -j 8 {protocol}", shell=True)
            p2 = subprocess.Popen(f"ssh -l ubuntu {remote_server1} 'cd {remote_root_folder}; make -j 8 {protocol}'", shell=True)
            p3 = subprocess.Popen(f"ssh -l ubuntu {remote_server2} 'cd {remote_root_folder}; make -j 8 {protocol}'", shell=True)
            p1.wait()
            p2.wait()
            p3.wait()
            
            for M, N in zip(spdz_config_dict[task_name]["M_list"], spdz_config_dict[task_name]["N_list"]):
                elog = f"{task_record_folder}spdz-{task_name}-{N}-{M}-1-1.log"
                target_execution_args = f"{source_file}-{task_name}-{N}-{M}-1-1 {protocol} {task_record_folder} {elog}"
                monitor_log = f"{task_record_folder}monitor-spdz-{task_name}-{N}-{M}-1-1.log"
                
                monitor.start_all(interface=NETWORK_INTERFACE)
                os.system(f"{execution_shell} {target_execution_args}")
                print(f"{execution_shell} {target_execution_args}")
                sch_file = f"{schedule_folder}{source_file}-{task_name}-{N}-{M}-1-1.sch"
                if(not os.path.exists(sch_file)):
                    print(f"{sch_file} does not exist.")
                monitor.stop_and_output(monitor_log)
                
                # exit(0)
    