from task_execution import *
from system_monitor import *
import pandas as pd
import glob

multi_step_tasks = ["sort"]

def task_analysis(task, save_folder, task_type="exhv"):
    
    performance_dict = {
        "task": [],
        "n": [],
        "m": [],
        "c": [],
        "B": [],
        "time_setup": [],
        "time_data_prepare": [],
        "subTask": [],
        "combine": []
    }
 
    task_log_folder = root_folder + "Record/" + task + "/"
    if(not os.path.exists((task_log_folder))):
        raise FileNotFoundError("The log folder does not exist.")

    # iterate all the files in task_log_folder
    for filename in os.listdir(task_log_folder):
        file_path = os.path.join(task_log_folder, filename)
        
        if(task_type == "exhv"):
            name_pattern = r'exhv-log-config-N=(\d+)-M=(\d+)-c=(\d+)-B=(\d+)\.log'
        elif(task_type == "pta"):
            name_pattern = r'pta-log-config-N=(\d+)-M=(\d+)-c=(\d+)-B=(\d+)\.log'
        elif(task_type == "config"):
            name_pattern = r'config-log-config-N=(\d+)-M=(\d+)-c=(\d+)-B=(\d+)\.log'
        
        task_match = re.match(name_pattern, filename)
        print(f"filename = {filename} | task_match = {task_match}")
        
        if task_match:
            N, M, c, B = task_match.groups()
            performance_dict["task"].append(task)
            performance_dict["n"].append(int(N))
            performance_dict["m"].append(int(M))
            performance_dict["c"].append(int(c))
            performance_dict["B"].append(int(B))
            
            with open(file_path, 'r') as file:
                content = file.read()
            
            if(task in multi_step_tasks):
                result = {
                    "time_setup": [],
                    "time_data_prepare": [],
                    "subTask": [],
                    "combine": []
                }
                patterns = {
                    "time_setup": r'time_setup:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)',
                    "time_data_prepare": r'time_data_prepare:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)',
                    "subTask": r'subTask:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)',
                    "combine": r'combine:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)'
                }
                
                for key, pattern in patterns.items():
                    matches = re.findall(pattern, content)
                    result[key].extend(float(match[0]) for match in matches)
                
                print(result)
                
                for key in result.keys():
                    performance_dict[key].append(np.sum(result[key]))
                
            else:
                performance_pattern = r'(time_setup|time_data_prepare|subTask|combine):\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)'
                matches = re.findall(performance_pattern, content)
                if(len(matches) != 4):
                    print(f"filename = {filename} | matches = {matches}")
                for match in matches:
                    performance_dict[match[0]].append(float(match[1]))
        
    performance_df = pd.DataFrame(performance_dict)
    performance_df = performance_df.sort_values(by=['n', 'm', 'c', 'B'])
    performance_df.to_excel(os.path.join(save_folder, f"{task}.xlsx"))
    
    return performance_df

def task_monitor_analysis(task, save_folder, task_type="exhv"):
    task_log_folder = root_folder + "Record/" + task + "/"
    if(not os.path.exists((task_log_folder))):
        raise FileNotFoundError("The log folder does not exist.")
    
    for filename in os.listdir(task_log_folder):
        file_path = os.path.join(task_log_folder, filename)

        if(task_type == "exhv"):
            monitor_name_pattern = r'monitor-exhv-log-config-N=(\d+)-M=(\d+)-c=(\d+)-B=(\d+)\.log'
        elif(task_type == "pta"):
            monitor_name_pattern = r'monitor-pta-log-config-N=(\d+)-M=(\d+)-c=(\d+)-B=(\d+)\.log'
        elif(task_type == "config"):
            monitor_name_pattern = r'monitor-config-log-config-N=(\d+)-M=(\d+)-c=(\d+)-B=(\d+)\.log'
        
            
        monitor_match = re.match(monitor_name_pattern, filename)
        
        if monitor_match:
            
            usage_dict = get_usage_dict(file_path)
            
            monitor_task_name = f"{task}_monitor_N={monitor_match.group(1)}_M={monitor_match.group(2)}_c={monitor_match.group(3)}_B={monitor_match.group(4)}"
            
            draw_usage_graph(usage_dict, os.path.join(save_folder, f"{monitor_task_name}.png"))
            
    return
       
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--pta', action='store_true', help='analyze pta task')
    parser.add_argument('--exhv', action='store_true', help='analyze exhv task')
    parser.add_argument('--config', action='store_true', help='analyze config task')
    parser.add_argument('--statistic', action='store_true', help='analyze all the tasks')
    args = parser.parse_args()
    
    if(args.pta):
        pta_prefix = "pta-log-config-"
        monitor_prefix = "monitor-pta-log-config-"
        
        analysis_folder = root_folder + "Record/Analysis/"
        if(not os.path.exists(analysis_folder)):
            os.mkdir(analysis_folder)
        
        monitor_analysis_folder = root_folder + "Record/Monitor/"
        if(not os.path.exists(monitor_analysis_folder)):
            os.mkdir(monitor_analysis_folder)
        
        task_dfs = []
        for task in target_tasks.keys():
                
            df = task_analysis(task, analysis_folder, task_type="pta")
            task_monitor_analysis(task, monitor_analysis_folder, task_type="pta")
            task_dfs.append(df)
        # concate all the task analysis xlsx into one.
        performance_df = pd.concat(task_dfs)
        performance_df.to_excel(os.path.join(analysis_folder, "performance_analysis.xlsx"))

        exit(0)

        
    if(args.exhv):
        exhv_prefix = "exhv-log-config-"
        monitor_prefix = "monitor-exhv-log-config-"
        
        analysis_folder = root_folder + "Record/Analysis/"
        if(not os.path.exists(analysis_folder)):
            os.mkdir(analysis_folder)
        
        monitor_analysis_folder = root_folder + "Record/Monitor/"
        if(not os.path.exists(monitor_analysis_folder)):
            os.mkdir(monitor_analysis_folder)
        
        task_dfs = []
        for task in target_tasks.keys():
            
            df = task_analysis(task, analysis_folder, task_type="exhv")
            task_monitor_analysis(task, monitor_analysis_folder)
            task_dfs.append(df)
        
        # concate all the task analysis xlsx into one. 
        performance_df = pd.concat(task_dfs)
        performance_df.to_excel(os.path.join(analysis_folder, "performance_analysis.xlsx"))
        
        exit(0)
    
    if(args.config):
        config_prefix = "config-log-config-"
        monitor_prefix = "monitor-config-log-config-"
        
        analysis_folder = root_folder + "Record/Analysis/"
        if(not os.path.exists(analysis_folder)):
            os.mkdir(analysis_folder)
        monitor_analysis_folder = root_folder + "Record/Monitor/"
        if(not os.path.exists(monitor_analysis_folder)):
            os.mkdir(monitor_analysis_folder)
        
        task_dfs = []
        for task in best_configurations.keys():
                
            df = task_analysis(task, analysis_folder, task_type="config")
            task_monitor_analysis(task, monitor_analysis_folder, task_type="config")
            task_dfs.append(df)
        
        # concate all the task analysis xlsx into one.
        performance_df = pd.concat(task_dfs)
        performance_df.to_excel(os.path.join(analysis_folder, "performance_analysis.xlsx"))
        exit(0)
    
    
    if(args.statistic):
        
        target_folder_prefix = os.path.join(root_folder, "Record_")
        folders = glob.glob(target_folder_prefix + '*')
        task_analysis_folder = root_folder + "Record/Analysis/"
        
        if(not os.path.exists(task_analysis_folder)):
            os.makedirs(task_analysis_folder)
        
        df_dict_list = []
        
        for k in range(len(folders)):
            folder = folders[k] 
            for filename in os.listdir(os.path.join(folder, "Analysis/")):
                if "performance_analysis.xlsx" in filename:
                    target_file = os.path.join(folder, "Analysis/", filename)
                    tmp_df = pd.read_excel(target_file)
                    df_dict_list.append(tmp_df)

        task_num = df_dict_list[0]["task"].to_numpy()
        n_list = df_dict_list[0]["n"].to_numpy()
        m_list = df_dict_list[0]["m"].to_numpy()
        c_list = df_dict_list[0]["c"].to_numpy()
        B_list = df_dict_list[0]["B"].to_numpy()

        # analyze the mean, min, max, std of the performance.
        time_setup_mat = np.array([df_dict_list[i]["time_setup"].to_numpy() for i in range(len(df_dict_list))])
        time_data_prepare_mat = np.array([df_dict_list[i]["time_data_prepare"].to_numpy() for i in range(len(df_dict_list))])
        time_subTask_mat = np.array([df_dict_list[i]["subTask"].to_numpy() for i in range(len(df_dict_list))])
        time_combine_mat = np.array([df_dict_list[i]["combine"].to_numpy() for i in range(len(df_dict_list))])

        mean_time_setup = np.mean(time_setup_mat, axis=0)
        std_time_setup = np.std(time_setup_mat, axis=0)
        mean_time_data = np.mean(time_data_prepare_mat, axis=0)
        std_time_data = np.std(time_data_prepare_mat, axis=0)
        mean_time_subTask = np.mean(time_subTask_mat, axis=0)
        std_time_subTask = np.std(time_subTask_mat, axis=0)
        mean_time_combine = np.mean(time_combine_mat, axis=0)
        std_time_combine = np.std(time_combine_mat, axis=0)
        total_time_mat = time_setup_mat + time_data_prepare_mat + time_subTask_mat + time_combine_mat
        mean_total_time = np.mean(total_time_mat, axis=0)
        std_total_time = np.std(total_time_mat, axis=0)

        # get the task analysis df.
        task_analysis_df = pd.DataFrame({"task": task_num, "n": n_list, "m": m_list, "c": c_list, "B": B_list,
                                        "mean_time_setup": mean_time_setup, "std_time_setup": std_time_setup, 
                                        "mean_time_data_prepare": mean_time_data, "std_time_data_prepare": std_time_data,
                                        "mean_time_subTask": mean_time_subTask, "std_time_subTask": std_time_subTask,
                                        "mean_time_combine": mean_time_combine, "std_time_combine": std_time_combine,
                                        "mean_total_time": mean_total_time,
                                        "std_total_time": std_total_time})
        task_analysis_df.to_excel(os.path.join(task_analysis_folder, "task_analysis.xlsx"))