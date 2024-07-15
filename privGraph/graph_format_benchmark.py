"""Run basic performance benchmark tests.
"""

import time
import argparse
import json
import os
from utils import get_k

MPI = True 
MPI_TASK = 4

MAIN_FOLDER = "/root/aby3/aby3-GraphQuery"

# gtype_list = ["random", "star", "powerlaw", "bipartite", "tree"]
gtype_list = ["random"]
# gtype_list = ["geometric", "powerlaw", "bipartite", "k_regular", "random"]


format_configs = {
    "privGraph": {
        "prefix": MAIN_FOLDER + "/data/baseline/",
        "record_folder": MAIN_FOLDER + "/record/privGraph/",
        "record_pattern": "(.*?)_n-(\d+)_k-(\d+)-(\d+)",
        # "n": [1024, 2048, 4096, 8192, 16384, 32768],
        "n": [1024],
        "e": -1,
        "k": [32, 32, 32, 32, 32, 32],
        "n_stash_size": [1024, 1024, 1024, 1024, 1024, 1024],
        "n_pack_size": [16, 16, 16, 16, 16, 16],
        "e_stash_size": [1024, 1024, 1024, 1024, 1024, 1024],
        "e_pack_size": [32, 32, 32, 32, 32, 32],
        "config_keys" : ["gtype", "n", "e", "k", "se", "pe", "sn", "pn"],
        "performance_keys": ["GraphLoad", "EdgeOramInit", "NodeOramInit", "EdgeExistQuery", "OuttingEdgesCountQuery", "NeighborsGetQuery", "GraphLoad_recv", "EdgeOramInit_recv", "NodeOramInit_recv", "EdgeExistQuery_recv", "OuttingEdgesCountQuery_recv", "NeighborsGetQuery_recv", "GraphLoad_send", "EdgeOramInit_send", "NodeOramInit_send", "EdgeExistQuery_send", "OuttingEdgesCountQuery_send", "NeighborsGetQuery_send"],
    },
    "adjmat":{
        "prefix": MAIN_FOLDER + "/data/baseline/",
        "record_folder": MAIN_FOLDER + "/record/adjmat/",
        "record_pattern": "(.*?)_n-(\d+)-(\d+)",
        "n": [1024, 2048, 4096, 8192, 16384],
        "e": -1,
        "n_stash_size": [1024, 1024, 1024, 1024, 1024],
        "n_pack_size": [16, 16, 16, 16, 16],
        "e_stash_size": [1024, 1024, 1024, 1024, 1024],
        "e_pack_size": [32, 32, 32, 32, 32],
        "config_keys" : ["gtype", "n", "e", "se", "pe", "sn", "pn"],
        "performance_keys": ["GraphLoad", "EdgeOramInit", "NodeOramInit", "EdgeExistQuery", "OuttingEdgesCountQuery", "NeighborsGetQuery", "GraphLoad_recv", "EdgeOramInit_recv", "NodeOramInit_recv", "EdgeExistQuery_recv", "OuttingEdgesCountQuery_recv", "NeighborsGetQuery_recv", "GraphLoad_send", "EdgeOramInit_send", "NodeOramInit_send", "EdgeExistQuery_send", "OuttingEdgesCountQuery_send", "NeighborsGetQuery_send"],
    },
    "edgelist":{
        "prefix": MAIN_FOLDER + "/data/baseline/",
        "record_folder": MAIN_FOLDER + "/record/edgelist/",
        "record_pattern": "(.*?)_n-(\d+)-(\d+)",
        # "n": [1024, 2048, 4096, 8192, 16384, 32768],
        "n": [1024],
        "e": -1,
        "config_keys" : ["gtype", "n", "e"],
        "performance_keys": ["GraphLoad", "EdgeExistQuery", "OuttingEdgesCountQuery", "NeighborsGetQuery", "GraphLoad_recv", "EdgeExistQuery_recv", "OuttingEdgesCountQuery_recv", "NeighborsGetQuery_recv", "GraphLoad_send", "EdgeExistQuery_send", "OuttingEdgesCountQuery_send", "NeighborsGetQuery_send"],
    },
}

def data_synchronize():
    os.system("ssh aby31 'rm -rf /root/aby3/aby3-GraphQuery/data/baseline/*'")
    os.system("ssh aby32 'rm -rf /root/aby3/aby3-GraphQuery/data/baseline/*'")
    os.system("scp -r /root/aby3/aby3-GraphQuery/data/baseline/ aby31:/root/aby3/aby3-GraphQuery/data/")
    os.system("scp -r /root/aby3/aby3-GraphQuery/data/baseline/ aby32:/root/aby3/aby3-GraphQuery/data/")
    

REPEAT_TIMES = 1

if __name__ == "__main__":
    
    # get the test configs. 
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', type=str, default="privGraph", help="target graph format")
    args = parser.parse_args()
    target = args.target
    
    # prepare the file. 
    if(MPI):
        os.system("cp ./frontend/main.pgpmpi ./frontend/main.cpp; python build.py --MPI")
    else:
        os.system("cp ./frontend/main.pgp ./frontend/main.cpp; python build.py")
    
    # run the privGraph profilings.
    target_config = format_configs[target] 
    
    if(not os.path.exists(target_config["record_folder"])):
        os.makedirs(target_config["record_folder"])
    if(not os.path.exists(target_config["prefix"])):
        os.makedirs(target_config["prefix"])
    e = -1
    
    if target == "privGraph":
    
        for gtype in gtype_list:
            for i in range(len(target_config["n"])):
                n, k = target_config["n"][i], target_config["k"][i]
                n_stash_size, n_pack_size = target_config["n_stash_size"][i], target_config["n_pack_size"][i]
                e_stash_size, e_pack_size = target_config["e_stash_size"][i], target_config["e_pack_size"][i]

                # self define k.
                edge_list_meta_file = f"{target_config['prefix']}{gtype}_n-{n}_edge_list_meta.txt"
                if not os.path.exists(edge_list_meta_file):
                    print(f"File {edge_list_meta_file} does not exist. Generate the data...")
                    generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {target_config['prefix']}{gtype}_n-{n} --type {gtype} --n {n} --e {e} --saving_type edgelist"
                    os.system(generate_command)

                with open(edge_list_meta_file, "r") as f:
                    edge_list_info = f.readline()
                    numbers = edge_list_info.split(" ")
                    V, E = int(numbers[0]), int(numbers[1])
                    assert V == n
                    k = get_k(V, E)
                    
                    print(f"k = {k}")
                
                # corresponsing data files.
                data_prefix = f"{gtype}_n-{n}_k-{k}"
                file_prefix = f"{target_config['prefix']}{data_prefix}"
                meta_file = f"{file_prefix}_meta.txt"
                graph_data_file = f"{file_prefix}_2dpartition.txt"
                
                # if data do not exist, generate the data.
                if not os.path.exists(meta_file):
                    print(f"File {meta_file} does not exist. Generate the data...")
                    generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --k {k}"
                    os.system(generate_command)
                    
                # data synchronization.
                data_synchronize()
                    
                # run the profiling.
                for j in range(REPEAT_TIMES):
                    count = j+1
                    run_args = f" -privGraph -prefix {data_prefix} -rcounter {count} -noram_stash_size {n_stash_size} -noram_pack_size {n_pack_size} -eoram_stash_size {e_stash_size} -eoram_pack_size {e_pack_size}"
                    
                    print(f"Run the privGraph profiling with {run_args}")
                    
                    if(MPI):
                        os.system(f"./Eval/mpi_dis_exec.sh \"{run_args}\" {MPI_TASK}")
                    else:
                        os.system(f"./Eval/dis_exec.sh \"{run_args}\"")

                    # show some debug info. 
                    os.system(f"cat ./debug.txt; rm ./debug.txt")
                    # os.system(f"cat ./debug.txt")
        
    # run the adjmat profilings.
    if target == "adjmat":
        for gtype in gtype_list:
            for i in range(len(target_config["n"])):
                n = target_config["n"][i]
                n_stash_size, n_pack_size = target_config["n_stash_size"][i], target_config["n_pack_size"][i]
                e_stash_size, e_pack_size = target_config["e_stash_size"][i], target_config["e_pack_size"][i]
                
                data_prefix = f"{gtype}_n-{n}"
                file_prefix = f"{target_config['prefix']}{data_prefix}"
                meta_file = f"{file_prefix}_edge_list_meta.txt"
                graph_data_file = f"{file_prefix}_edge_list.txt"
                
                # if data do not exist, generate the data.
                if not os.path.exists(meta_file):
                    print(f"File {meta_file} does not exist. Generate the data...")
                    generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --saving_type edgelist"
                    os.system(generate_command)
                data_synchronize()
                
                # run the profiling.
                for j in range(REPEAT_TIMES):
                    count = j+1
                    run_args = f" -adjmat true -prefix {data_prefix} -rcounter {count} -noram_stash_size {n_stash_size} -noram_pack_size {n_pack_size} -eoram_stash_size {e_stash_size} -eoram_pack_size {e_pack_size}"
                    
                    print(f"Run the adjmat profiling with {run_args}")
                    os.system(f"./Eval/dis_exec.sh \"{run_args}\"")

                    # show some debug info. 
                    os.system(f"cat ./debug.txt; rm ./debug.txt")
        
    # run the edgelist profilings.
    if target == "edgelist":
        for gtype in gtype_list:
            for i in range(len(target_config["n"])):
                n = target_config["n"][i]
                
                data_prefix = f"{gtype}_n-{n}"
                file_prefix = f"{target_config['prefix']}{data_prefix}"
                meta_file = f"{file_prefix}_edge_list_meta.txt"
                graph_data_file = f"{file_prefix}_edge_list.txt"
                
                # if data do not exist, generate the data.
                if not os.path.exists(meta_file):
                    print(f"File {meta_file} does not exist. Generate the data...")
                    generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --saving_type edgelist"
                    os.system(generate_command)
                # data_synchronize()
                
                # run the profiling.
                for j in range(REPEAT_TIMES):
                    count = j+1
                    run_args = f" -edgelist true -prefix {data_prefix} -rcounter {count}"
                    
                    print(f"Run the edgelist profiling with {run_args}")
                    # os.system(f"./Eval/dis_exec.sh \"{run_args}\"")

                    # show some debug info. 
                    if(MPI):
                        os.system(f"./Eval/mpi_dis_exec.sh \"{run_args}\" {MPI_TASK}")
                    else:
                        os.system(f"./Eval/dis_exec.sh \"{run_args}\"")
                    os.system(f"cat ./debug.txt; rm ./debug.txt")
                

    
    