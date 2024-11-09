import argparse
import re
import os
import pandas as pd
from spdz_task_execution import spdz_config_dict, tgl_config_dict
from system_monitor import get_usage_dict, draw_usage_graph


def spdz_log_total_time_extraction(logging_file):
    with open(logging_file, 'r') as file:
        content = file.read()
    total_time_match = r'Time = ([\d\.]+(?:e[+-]?\d+)?) seconds'
    matches = re.findall(total_time_match, content)
    if(len(matches) != 1):
        print(f"filename = {logging_file} | matches = {matches}")
        matches = [-1]
    total_time = float(matches[0])
    
    return total_time
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--log_folder", type=str, required=True, help="The task to be executed.")
    parser.add_argument("--spdz", action='store_true', help='analyze spdz origional task')
    parser.add_argument("--tgl", action='store_true', help='analyze spdz origional task')
    args = parser.parse_args()
    
    log_folder = args.log_folder
    if(not os.path.exists(log_folder)):
        raise FileNotFoundError("The log folder does not exist.")
    
    if(args.tgl):
        record_dict = {
            "task": [],
            "N": [],
            "M": [],
            "task_num": [],
            "total_time": []
        }
        
        monitor_graph_folder = log_folder + "monitor_graph/"
        if(not os.path.exists(monitor_graph_folder)):
            os.makedirs(monitor_graph_folder)
            
        for task in tgl_config_dict.keys():
            task_record_folder = log_folder + task + "/"
            
            if(not os.path.exists(task_record_folder)):
                raise FileNotFoundError("The task folder does not exist.")
            
            for M, N in zip(tgl_config_dict[task]["M_list"], tgl_config_dict[task]["N_list"]):
                for task_num in tgl_config_dict[task]["c_list"]:
                    elog_file = f"{task_record_folder}tgl-{task}-{N}-{M}-1-{task_num}.log"
                    monitor_log = f"{task_record_folder}monitor-tgl-{task}-{N}-{M}-1-{task_num}.log"
                    
                    if(not os.path.exists(elog_file)):
                        raise FileNotFoundError(f"The log file {elog_file} does not exist.")
                
                    total_time = spdz_log_total_time_extraction(elog_file)
                
                    record_dict["task"].append(task)
                    record_dict["N"].append(N)
                    record_dict["M"].append(M)
                    record_dict["task_num"].append(task_num)
                    record_dict["total_time"].append(total_time)
                    
                    usage_dict = get_usage_dict(monitor_log)
                    monitor_task_name = f"tgl_{task}_monitor_N={N}_M={M}_C={task_num}.png"
                    draw_usage_graph(usage_dict, os.path.join(monitor_graph_folder, monitor_task_name))
        
        tgl_performance_df = pd.DataFrame(record_dict)
        tgl_performance_df = tgl_performance_df.sort_values(by=['task', 'N', 'M', "task_num"])
        tgl_performance_df.to_excel(os.path.join(log_folder, "tgl_performance_analysis.xlsx"))
            

    if(args.spdz):
        
        record_dict = {
            "task": [],
            "N": [],
            "M": [],
            "total_time": []
        }
        
        monitor_graph_folder = log_folder + "monitor_graph/"
        if(not os.path.exists(monitor_graph_folder)):
            os.makedirs(monitor_graph_folder)
        
        for task in spdz_config_dict.keys():
            
            task_record_folder = log_folder + task + "/"
            
            if(not os.path.exists(task_record_folder)):
                raise FileNotFoundError("The task folder does not exist.")
            
            for M, N in zip(spdz_config_dict[task]["M_list"], spdz_config_dict[task]["N_list"]):
                elog_file = f"{task_record_folder}spdz-{task}-{N}-{M}-1-1.log"
                monitor_log = f"{task_record_folder}monitor-spdz-{task}-{N}-{M}-1-1.log"
                if(not os.path.exists(elog_file)):
                    raise FileNotFoundError(f"The log file {elog_file} does not exist.")

                with open(elog_file, 'r') as file:
                    content = file.read()
                total_time = spdz_log_total_time_extraction(elog_file)
                
                record_dict["task"].append(task)
                record_dict["N"].append(N)
                record_dict["M"].append(M)
                record_dict["total_time"].append(total_time)
                
                usage_dict = get_usage_dict(monitor_log)
                monitor_task_name = f"spdz_{task}_monitor_N={N}_M={M}.png"
                draw_usage_graph(usage_dict, os.path.join(monitor_graph_folder, monitor_task_name))
        
        spdz_performance_df = pd.DataFrame(record_dict)
        spdz_performance_df = spdz_performance_df.sort_values(by=['task', 'N', 'M'])
        spdz_performance_df.to_excel(os.path.join(log_folder, "spdz_performance_analysis.xlsx"))
        
                
                
                
                
    
    