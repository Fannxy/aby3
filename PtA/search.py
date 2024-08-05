from system_profile import *
from system_monitor import SystemMonitor

N_list = [1<<10, 1<<15, 1<<20, 1<<25, 1<<30, 1<<32]
# N_list = [1<<10]
Func_list = ["mcomp", "comp", "mtag", "tag", "subH"]
record_folder = root_folder + "PtA/Record_search/"
exec_sh = root_folder + "Eval/basic/dis_mpi_exec.sh"

large_n = 1<<30
threshold_list = [1<<4, 1<<8, 1<<12, 1<<16, 1<<20]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    
    # compile the project.
    os.system(f"cp {root_folder}frontend/main.bench {root_folder}frontend/main.cpp")
    os.system(f"python {root_folder}build.py")
    
    if(not os.path.exists(record_folder)):
        os.makedirs(record_folder)
    
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
    
    df = pd.DataFrame(df_dict_list).to_excel(record_folder + "search.xlsx", index=False)

    
    # test the threshold.
    target_func = "subH"
    targetN = large_n   
    for threshold in threshold_list:
        record_file = record_folder + f"search-{target_func}-N={targetN}-threshold={threshold}.log"
        exec_args = f" -search -N {targetN} -Func {target_func} -threshold {threshold} -record_file {record_file}"
        os.system(f"{exec_sh} 1 \"{exec_args}\"")
        print(exec_args)
    
    # obtain the result.
    df_dict_list = {"Threshold": [], "Time": []}
    for threshold in threshold_list:
        record_file = record_folder + f"search-{target_func}-N={targetN}-threshold={threshold}.log"
        with open(record_file, "r") as f:
            last_line = f.readlines()[-1]
            match = re.search(r': ([\d.]+) milliseconds', last_line)
            comp_time = float(match.group(1)) 
        df_dict_list["Time"].append(comp_time)
        df_dict_list["Threshold"].append(threshold)
    
    df = pd.DataFrame(df_dict_list).to_excel(record_folder + "search_threshold.xlsx", index=False)