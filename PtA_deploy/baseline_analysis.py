from get_execution_parameters import *
from system_profile import *
from system_monitor import SystemMonitor, get_usage_dict, draw_usage_graph
from baseline_execution import *

if __name__ == "__main__":
        
    parser = argparse.ArgumentParser()
    parser.add_argument('--sqrtOram', action='store_true', help='run sqrt oram baseline')
    args = parser.parse_args()
    
    if(args.sqrtOram):
        log_folder = root_folder + "Record_sqrtOram/"
        if(not os.path.exists(log_folder)):
            raise Exception("Please run the baseline_execution.py --sqrtOram first.")
        
        analysis_folder = log_folder + "Analysis/"
        if(not os.path.exists(analysis_folder)):
            os.makedirs(analysis_folder)
        
        monitor_folder = log_folder + "Monitor/"
        if(not os.path.exists(monitor_folder)):
            os.makedirs(monitor_folder)
            
        analysis_dict = {
            "n": [],
            "m": [],
            "time_setup": [],
            "time_data": [],
            "time_construction": [],
            "time_access": [],
            "time_access_amortized": [],
        }
        
        log_pattern = r"sqrtOram-log-M=(\d+)-N=(\d+).log"
        key_patterns = {
            "time_setup": r"time_setup:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)",
            "time_data": r"data_load:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)",
            "time_construction": r"oram_construction:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)",
            "time_access": r"oram_access:\s+([\d\.]+(?:e[+-]?\d+)?)\s+(milliseconds|ms|micros)"
        }
        
        for filename in os.listdir(log_folder):
            task_match = re.match(log_pattern, filename)
            file_path = os.path.join(log_folder, filename)
            print(f"Processing {filename}")
            print(f"logging folder = {log_folder}")
            
            # get the performance metrics. 
            if task_match:
                with open(file_path, 'r') as file:
                    content = file.read()
                m, n = int(task_match.group(1)), int(task_match.group(2))
                analysis_dict["n"].append(n)
                analysis_dict["m"].append(m)
                for key, pattern in key_patterns.items():
                    matches = re.findall(pattern, content)
                    print(matches[0][0])
                    analysis_dict[key].append(float(matches[0][0]))
            
                # draw the usage graphes. 
                monitor_name = f"monitor-sqrtOram-log-M={m}-N={n}.log"
                usage_dict = get_usage_dict(os.path.join(log_folder, monitor_name))
                draw_usage_graph(usage_dict, os.path.join(monitor_folder, f"monitor-sqrtOram-log-M={m}-N={n}.png"))
 
        
        analysis_dict["time_access_amortized"] = (np.array(analysis_dict["time_access"]) + np.array(analysis_dict["time_construction"]) + np.array(analysis_dict["time_data"]) + np.array(analysis_dict["time_setup"])) / np.array(analysis_dict["n"])

        performance_df = pd.DataFrame(analysis_dict)
        performance_df = performance_df.sort_values(by=['n', 'm'])
        performance_df.to_excel(os.path.join(analysis_folder, "sqrtOram.xlsx"))