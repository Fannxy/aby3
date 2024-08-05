import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob

from system_profile import *
from system_monitor import draw_usage_graph, get_usage_dict

network_analysis_folder = system_profiling_folder + "network_analysis/"
resource_analysis_folder = system_profiling_folder + "resource_analysis/"

def network_analysis():
    target_folder_prefix = os.path.join(root_folder, "system_profiling")
    folders = glob.glob(target_folder_prefix + '*')
    
    if(not os.path.exists(network_analysis_folder)):
        os.makedirs(network_analysis_folder)
        
    df_dict_list = []
    # for folder in folders:
    for k in range(len(folders)):
        folder = folders[k]
        target_file = os.path.join(folder, "network.json")
        if(os.path.exists(target_file) == False):
            continue
        with open(target_file, 'r') as file:
            _task_dict = json.load(file)
        
        df_dict = {"task_num": [], 
                   "latency": [],
                   "vector_size": [],
                   "communication_time": [],
                   "bandwidth": []}

        first_key = next(iter(_task_dict.keys()))
        unit_len = len(_task_dict[first_key]['vector_size_list'])
        for key in _task_dict.keys():
            for i in range(unit_len):
                df_dict["task_num"].append(key)
                df_dict["latency"].append(_task_dict[key]['latency'])
                df_dict["vector_size"].append(_task_dict[key]['vector_size_list'][i])
                df_dict["communication_time"].append(_task_dict[key]['communication_time_list'][i])
                df_dict["bandwidth"].append(_task_dict[key]['bandwidth_list'][i])
        
        df = pd.DataFrame.from_dict(df_dict)
        df.to_excel(os.path.join(network_analysis_folder, f"network-{k}.xlsx"))
        df_dict_list.append(df) 
    
    task_num = df_dict_list[0]["task_num"].to_numpy()
    vector_size = df_dict_list[0]["vector_size"].to_numpy()
    # analyze the mean, min, max, std of the communication time, bandwidth, latency.
    latency_mat = np.array([df_dict_list[i]["latency"].to_numpy() for i in range(len(df_dict_list))])
    std_latency = np.std(latency_mat, axis=0)
    mean_latency = np.mean(latency_mat, axis=0)
    min_latency = np.min(latency_mat, axis=0)
    max_latency = np.max(latency_mat, axis=0)
    
    communication_time_mat = np.array([df_dict_list[i]["communication_time"].to_numpy() for i in range(len(df_dict_list))])
    std_communication_time = np.std(communication_time_mat, axis=0)
    mean_communication_time = np.mean(communication_time_mat, axis=0)
    min_communication_time = np.min(communication_time_mat, axis=0)
    max_communication_time = np.max(communication_time_mat, axis=0)
    
    bandwidth_mat = np.array([df_dict_list[i]["bandwidth"].to_numpy() for i in range(len(df_dict_list))])
    std_bandwidth = np.std(bandwidth_mat, axis=0)
    mean_bandwidth = np.mean(bandwidth_mat, axis=0)
    min_bandwidth = np.min(bandwidth_mat, axis=0)
    max_bandwidth = np.max(bandwidth_mat, axis=0)
    
    # get the network analysis df.
    network_df = pd.DataFrame({"task_num": task_num, "vector_size": vector_size, 
                               "mean_latency": mean_latency, "std_latency": std_latency, 
                               "min_latency": min_latency, "max_latency": max_latency,
                               "mean_communication_time": mean_communication_time, "std_communication_time": std_communication_time, 
                               "min_communication_time": min_communication_time, "max_communication_time": max_communication_time,
                               "mean_bandwidth": mean_bandwidth, "std_bandwidth": std_bandwidth, 
                               "min_bandwidth": min_bandwidth, "max_bandwidth": max_bandwidth})
    network_df.to_excel(os.path.join(network_analysis_folder, "network_analysis.xlsx"))
    print(network_df)
        
        
def monitor_analysis():
    target_folder_prefix = os.path.join(root_folder, "system_profiling")
    folders = glob.glob(target_folder_prefix + '*')
    
    if(not os.path.exists(resource_analysis_folder)):
        os.makedirs(resource_analysis_folder)
    
    counter = 0
    
    for folder in folders:
        # traverse the folder
        saving_graph_folder = os.path.join(resource_analysis_folder, f"graph_{counter}")
        os.mkdir(saving_graph_folder)
        
        for filename in os.listdir(folder):
            match = re.match(r'profile_task_(\d+)\.log', filename)
            if match:
                task_num = int(match.group(1))
                
                usage_dict = get_usage_dict(os.path.join(folder, filename))
                graph_name = f"task_{task_num}"
                save_graph_name = os.path.join(saving_graph_folder, f"{graph_name}.png")
                
                draw_usage_graph(usage_dict, save_graph_name)
                
        counter += 1
                        
                    
if __name__ == "__main__":
    network_analysis()
    monitor_analysis()