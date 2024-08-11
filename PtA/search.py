from system_profile import *
from system_monitor import SystemMonitor

N_list = [1<<30]
Pta_tasks_list = [1<<5]
# Func_list = ["mcomp", "comp", "subH"]
Func_list = ["comp", "subH"]
record_folder = root_folder + "PtA/Record_search/"
exec_sh = root_folder + "Eval/basic/dis_mpi_exec.sh"

large_n = 1<<30
alpha_list = [1<<4, 1<<8, 1<<12, 1<<16, 1<<20]
threshold_list = [1<<2, 1<<4, 1<<8, 1<<12, 1<<14]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--test', type=str, help="which test you want run:\n 1-normal search;\n 2-normal search with threshold;\n 3-pta search;\n 4-pta search with threshold")
    args = parser.parse_args()
    
    # compile the project.
    os.system(f"cp {root_folder}frontend/main.bench {root_folder}frontend/main.cpp")
    os.system(f"python {root_folder}build.py")
    
    if(not os.path.exists(record_folder)):
        os.makedirs(record_folder)
        
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
        for N in N_list:
            for task_num in Pta_tasks_list:
                for func in Func_list:
                    record_file = record_folder + f"pta-{func}-N={N}-task_num={task_num}.log"
                    exec_args = f" -pta-search -N {N} -Func {func} -record_file {record_file} -optB 65536"
                    os.system(f"{exec_sh} {task_num} \"{exec_args}\"")
                    print(exec_args)
        
        # obtain the result.
        df_dict_list = {"N": [], "task_num": [], "Time": []}
        for N in N_list:
            for task_num in Pta_tasks_list:
                for func in Func_list:
                    record_file = record_folder + f"pta-{func}-N={N}-task_num={task_num}.log"
                    with open(record_file, "r") as f:
                        last_line = f.readlines()[-1]
                        match = re.search(r': ([\d.]+) milliseconds', last_line)
                        comp_time = float(match.group(1)) 
                    df_dict_list["Time"].append(comp_time)
                    df_dict_list["N"].append(N)
                    df_dict_list["task_num"].append(task_num)
        
        df = pd.DataFrame(df_dict_list).to_excel(record_folder + "pta_search.xlsx", index=False)
    
    if(args.test == "4" or args.test == "pta_threshold"):
        
        print("not implemented!")
        