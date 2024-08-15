from system_profile import *
from system_monitor import SystemMonitor

N_list = [1<<10, 1<<12, 1<<14, 1<<16, 1<<18, 1<<20, 1<<22, 1<<24, 1<<26, 1<<28, 1<<30]
repeat_list = [10, 10, 10, 10, 10, 5, 5, 5, 1, 1, 1]
# N_list = [1<<10, 1<<15, 1<<20, 1<<25, 1<<30]
# N_list = [1<<10, 1<<15]
tasks_list = [1, 1<<1, 1<<2, 1<<3, 1<<4, 24]

# N_list = [1<<20]
# repeat_list = [1]
# tasks_list = [24]
Func_list = ["mcomp", "comp"] # see the effects of tasks and N.
# Func_list = ["comp"]
record_folder = root_folder + "PtA/Record_search/"
gc_folder = root_folder + "PtA/Record_search/gc/"
monitor_folder = root_folder + "PtA/Record_search/monitor/"
exec_sh = root_folder + "Eval/basic/dis_mpi_exec.sh"

large_n = 1<<30
alpha_list = [1<<4, 1<<8, 1<<12, 1<<16, 1<<20]
threshold_list = [1<<2, 1<<4, 1<<8, 1<<12, 1<<14]

monitor = SystemMonitor(0.01)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--test', type=str, help="which test you want run:\n 1-normal search;\n 2-normal search with threshold;\n 3-pta search;\n 4-pta search with threshold;\n 5-ndss parallel search;\n 6-ndss parallel search with threshold.")
    args = parser.parse_args()
    
    # compile the project.
    os.system(f"cp {root_folder}frontend/main.bench {root_folder}frontend/main.cpp")
    os.system(f"python {root_folder}build.py")
    
    if(not os.path.exists(record_folder)):
        os.makedirs(record_folder)
    
    if(not os.path.exists(monitor_folder)):
        os.makedirs(monitor_folder)
    
    if(not os.path.exists(gc_folder)):  
        os.makedirs(gc_folder)
        
    if(args.test == "1" or args.test == "normal"):
    
        for N in N_list:
            for func in Func_list:
                record_file = record_folder + f"search-{func}-N={N}.log"
                exec_args = f" -search -N {N} -Func {func} -record_file {record_file}"
                os.system(f"{exec_sh} 1 \"{exec_args}\"")
                print(exec_args)
        
        # obtain the result.
        df_dict_list = {"N": [], "Func": [], "Time": []}
        for N in N_list:
            for func in Func_list:
                record_file = record_folder + f"search-{func}-N={N}.log"
                with open(record_file, "r") as f:
                    last_line = f.readlines()[-1]
                    match = re.search(r': ([\d.]+) milliseconds', last_line)
                    comp_time = float(match.group(1)) 
                df_dict_list["Time"].append(comp_time)
                df_dict_list["N"].append(N)
                df_dict_list["Func"].append(func)
        
        df = pd.DataFrame(df_dict_list).to_excel(record_folder + "normal_search.xlsx", index=False)
    
    if(args.test == "2" or args.test == "normal_threshold"):

        # test the threshold.
        target_func = "subH"
        targetN = large_n   
        for alpha in alpha_list:
            for threshold in threshold_list:
                record_file = record_folder + f"search-{target_func}-N={targetN}-threshold={threshold}-alpha={alpha}.log"
                exec_args = f" -search -N {targetN} -Func {target_func} -threshold {threshold} -alpha {alpha} -record_file {record_file}"
                print("exec_args: ", exec_args)
                os.system(f"{exec_sh} 1 \"{exec_args}\"")
                print(exec_args)
        
        # obtain the result.
        df_dict_list = {"Threshold": [], "Time": [], "alpha": [], "recursive_depth": []}
        for alpha in alpha_list:
            for threshold in threshold_list:
                record_file = record_folder + f"search-{target_func}-N={targetN}-threshold={threshold}-alpha={alpha}.log"
                with open(record_file, "r") as f:
                    last_line = f.readlines()[-1]
                    match = re.search(r': ([\d.]+) milliseconds', last_line)
                    comp_time = float(match.group(1)) 
                df_dict_list["Time"].append(comp_time)
                df_dict_list["Threshold"].append(threshold)
                df_dict_list["alpha"].append(alpha)
                df_dict_list["recursive_depth"].append(int((np.log2(large_n) - np.log2(threshold)) / np.log2(alpha)))
        
        df = pd.DataFrame(df_dict_list).to_excel(record_folder + "normal_search_threshold.xlsx", index=False)
    
    if(args.test == "3" or args.test == "pta"):
        # for N in N_list:
        for i in range(len(N_list)):
            N = N_list[i]
            repeat = repeat_list[i]
            for task_num in tasks_list:
                optB = int(N / task_num)+1 # the maximum optB, only using one round.
                # optB = 1<<25
                # optB = min(optB, 1<<22)
                for func in Func_list:
                    record_file = record_folder + f"pta-{func}-N={N}-task_num={task_num}.log"
                    monitor_file = monitor_folder + f"pta-{func}-N={N}-task_num={task_num}.log"
                    exec_args = f" -pta-search -N {N} -Func {func} -REPEAT {repeat} -record_file {record_file} -optB {optB}"
                    monitor.start_all("ens110")
                    os.system(f"{exec_sh} {task_num} \"{exec_args}\"")
                    monitor.stop_and_output(monitor_file)
                    usage_dict = get_usage_dict(monitor_file)
                    draw_usage_graph(usage_dict, monitor_file.replace(".log", ".png"))
                    print(exec_args)
        
        # obtain the result.
        df_dict_list = {"Func": [], "N": [], "task_num": [], "Time": [],  "optB": []}
        # for N in N_list:
        for i in range(len(N_list)):
            N = N_list[i]
            repeat = repeat_list[i]
            # if(N < 1<<28):
            #     continue
            for task_num in tasks_list:
                for func in Func_list:
                    record_file = record_folder + f"pta-{func}-N={N}-task_num={task_num}.log"
                    with open(record_file, "r") as f:
                        last_line = f.readlines()[-1]
                        match = re.search(r': ([\d.]+) milliseconds', last_line)
                        comp_time = float(match.group(1)) 
                    df_dict_list["Time"].append(comp_time / repeat)
                    df_dict_list["N"].append(N)
                    df_dict_list["task_num"].append(task_num)
                    df_dict_list["Func"].append(func)
                    df_dict_list["optB"].append(int(N / task_num)+1)
        df = pd.DataFrame(df_dict_list)
        df = df.sort_values(by=["Func", "N", "task_num"])
        df.to_excel(record_folder + "pta_search.xlsx", index=False)
        os.system(f"mv {record_folder}*.log {gc_folder}")
    
    if(args.test == "4" or args.test == "pta_threshold"):
        
        print("not implemented!")
    
    if(args.test == "5" or args.test == "ndss"):
        # for N in N_list:
        for i in range(len(N_list)):
            N = N_list[i]
            repeat = repeat_list[i]
            for task_num in tasks_list:
                for func in Func_list:
                    record_file = record_folder + f"ndss-{func}-N={N}-task_num={task_num}.log"
                    monitor_file = monitor_folder + f"ndss-{func}-N={N}-task_num={task_num}.log"
                    exec_args = f" -ndss-search -N {N} -Func {func} -REPEAT {repeat} -record_file {record_file}"
                    monitor.start_all("ens110")
                    os.system(f"{exec_sh} {task_num} \"{exec_args}\"")
                    monitor.stop_and_output(monitor_file)
                    usage_dict = get_usage_dict(monitor_file)
                    draw_usage_graph(usage_dict, monitor_file.replace(".log", ".png"))
                    print(exec_args)
        
        # obtain the result.
        df_dict_list = {"Func": [], "N": [],  "task_num": [], "Time": []}
        # for N in N_list:
        for i in range(len(N_list)):
            N = N_list[i]
            repeat = repeat_list[i]
            for task_num in tasks_list:
                for func in Func_list:
                    record_file = record_folder + f"ndss-{func}-N={N}-task_num={task_num}.log"
                    with open(record_file, "r") as f:
                        last_line = f.readlines()[-1]
                        match = re.search(r': ([\d.]+) milliseconds', last_line)
                        comp_time = float(match.group(1)) 
                    df_dict_list["Time"].append(comp_time / repeat)
                    df_dict_list["N"].append(N)
                    df_dict_list["task_num"].append(task_num)
                    df_dict_list["Func"].append(func)
                    
        df = pd.DataFrame(df_dict_list)
        df = df.sort_values(by=["Func", "N", "task_num"])
        df.to_excel(record_folder + "ndss_search.xlsx", index=False)   
        os.system(f"mv {record_folder}*.log {gc_folder}")