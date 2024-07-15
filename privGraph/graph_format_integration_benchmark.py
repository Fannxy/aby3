"""Run graph integration benchmark tests.
"""

import os
import argparse
from utils import get_k 

MAIN_FOLDER = "/root/aby3/aby3-GraphQuery"

# gtype_list = ["random", "geometric", "powerlaw", "bipartite"]
gtype_list = ["geometric"]
# data_provider_list = [2, 4, 8]
data_provider_list = [2]

format_configs = {
    "privGraph": {
        "prefix": MAIN_FOLDER + "/data/multiparty/",
        "record_folder": MAIN_FOLDER + "/record_offline/privGraph/",
        "record_pattern": "(.*?)_n-(\d+)_k-(\d+)_p-(\d+)-(\d+)",
        "n": [1024],
        "e": -1,
        "k": [32],
        "n_stash_size": [32],
        "n_pack_size": [16],
        "e_stash_size": [1024],
        "e_pack_size": [32],
        "config_keys" : ["gtype", "n", "e", "k", "se", "pe", "sn", "pn"],
        "performance_keys": ["GraphLoad", "EdgeOramInit", "NodeOramInit", "GraphLoad_recv", "EdgeOramInit_recv", "NodeOramInit_recv", "GraphLoad_send", "EdgeOramInit_send", "NodeOramInit_send"],
    },
    "adjmat":{
        "prefix": MAIN_FOLDER + "/data/multiparty/",
        "record_folder": MAIN_FOLDER + "/record_offline/adjmat/",
        "record_pattern": "(.*?)_n-(\d+)_p-(\d+)-(\d+)",
        "n": [1024],
        "e": -1,
        "n_stash_size": [32],
        "n_pack_size": [16],
        "e_stash_size": [1024],
        "e_pack_size": [32],
        "config_keys" : ["gtype", "n", "e", "se", "pe", "sn", "pn"],
        "performance_keys": ["GraphLoad", "EdgeOramInit", "NodeOramInit", "GraphLoad_recv", "EdgeOramInit_recv", "NodeOramInit_recv", "GraphLoad_send", "EdgeOramInit_send", "NodeOramInit_send"],
    },
    "edgelist":{
        "prefix": MAIN_FOLDER + "/data/multiparty/",
        "record_folder": MAIN_FOLDER + "/record_offline/edgelist/",
        "record_pattern": "(.*?)_n-(\d+)_p-(\d+)-(\d+)",
        "n": [1024],
        "e": -1,
        "config_keys" : ["gtype", "n", "e"],
        "performance_keys": ["GraphLoad", "GraphLoad_recv", "GraphLoad_send"],
    },
}

REPEAT_TIMES = 1

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--target', type=str, default="privGraph", help="target graph format")
    args = parser.parse_args()
    target = args.target

    os.system("cp ./frontend/main.pgp ./frontend/main.cpp; python build.py")

    target_config = format_configs[target] 

    if(not os.path.exists(target_config["record_folder"])):
        os.makedirs(target_config["record_folder"])
    if(not os.path.exists(target_config["prefix"])):
        os.makedirs(target_config["prefix"])
    e = -1

    if target == "privGraph":

        for gtype in gtype_list:
            for data_providers in data_provider_list:
                for i in range(len(target_config["n"])):
                    n, k = target_config["n"][i], target_config["k"][i]
                    n_stash_size, n_pack_size = target_config["n_stash_size"][i], target_config["n_pack_size"][i]
                    e_stash_size, e_pack_size = target_config["e_stash_size"][i], target_config["e_pack_size"][i]

                    # get k and generate the data.
                    data_exist = True
                    edge_list_meta_file = f"{target_config['prefix']}{gtype}_n-{n}_edge_list_meta_multiparty.txt"
                    if not os.path.exists(edge_list_meta_file):
                        data_exist = False

                    for p in range(data_providers):
                        party_meta_file = f"{target_config['prefix']}{gtype}_n-{n}_edge_list_meta_party-{p}.txt"
                        if not os.path.exists(party_meta_file):
                            data_exist = False

                    if not data_exist: # first generate the edge list data.
                        print(f"File {edge_list_meta_file} does not exist for {data_providers} parties. Generate the data...")
                        generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {target_config['prefix']}{gtype}_n-{n} --type {gtype} --n {n} --e {e} --p {data_providers} --saving_type edgelist"
                        os.system(generate_command)
                    
                    # get the total edge numbers e.
                    with open(edge_list_meta_file, "r") as f:
                        numbers = [int(num) for num in f.readline().strip().split(" ")]
                        v = numbers[0]
                        # n = numbers[1]
                        assert data_providers == numbers[1]
                        edges_num = 0
                        for i in range(data_providers):
                            edges_num += numbers[2 + i]
                        
                        k = get_k(v, edges_num)

                    # generate the 2d partition file.
                    data_prefix = f"{gtype}_n-{n}_k-{k}"
                    file_prefix = f"{target_config['prefix']}{data_prefix}"
                    meta_file = file_prefix + "_meta_multiparty.txt"

                    if not os.path.exists(meta_file):
                        print(f"File {meta_file} does not exist. Generate the data...")
                        generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --k {k} --p {data_providers} --saving_type 2dpartition"
                        os.system(generate_command)
                    
                    # run the profiling.
                    for j in range(REPEAT_TIMES):
                        count = j+1
                        run_args = f" -multi_privGraph -prefix {data_prefix} -rcounter {count} -noram_stash_size {n_stash_size} -noram_pack_size {n_pack_size} -eoram_stash_size {e_stash_size} -eoram_pack_size {e_pack_size}"

                        print(f"Run the privGraph profiling with {run_args}")
                        os.system(f"./Eval/dis_exec.sh \"{run_args}\"")
                        os.system(f"cat ./debug.txt; rm ./debug.txt")
                        record_file = f"{target_config['record_folder']}{data_prefix}-{count}.txt"
                        target_record_file = f"{target_config['record_folder']}{data_prefix}_p-{data_providers}-{count}.txt"
                        os.system(f"mv {record_file} {target_record_file}")
        
    # run the adjmat profilings.
    if target == "adjmat":

        for gtype in gtype_list:
            for data_providers in data_provider_list:
                for i in range(len(target_config["n"])):
                    n = target_config["n"][i]
                    n_stash_size, n_pack_size = target_config["n_stash_size"][i], target_config["n_pack_size"][i]
                    e_stash_size, e_pack_size = target_config["e_stash_size"][i], target_config["e_pack_size"][i]
                    
                    data_prefix = f"{gtype}_n-{n}"
                    file_prefix = f"{target_config['prefix']}{data_prefix}"
                    meta_file = f"{file_prefix}_edge_list_meta_multiparty.txt"
                    
                    # if data do not exist, generate the data.
                    if not os.path.exists(meta_file):
                        print(f"File {meta_file} does not exist. Generate the data...")
                        generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --p {data_providers} --saving_type edgelist"
                        os.system(generate_command)
                    
                    # run the profiling.
                    for j in range(REPEAT_TIMES):
                        count = j+1
                        run_args = f" -multi_adjmat true -prefix {data_prefix} -rcounter {count} -noram_stash_size {n_stash_size} -noram_pack_size {n_pack_size} -eoram_stash_size {e_stash_size} -eoram_pack_size {e_pack_size}"
                        
                        print(f"Run the adjmat profiling with {run_args}")
                        os.system(f"./Eval/dis_exec.sh \"{run_args}\"")
                        os.system(f"cat ./debug.txt; rm ./debug.txt")
                        record_file = f"{target_config['record_folder']}{data_prefix}-{count}.txt"
                        target_record_file = f"{target_config['record_folder']}{data_prefix}_p-{data_providers}-{count}.txt"
                        os.system(f"mv {record_file} {target_record_file}")
        
    # run the edgelist profilings.
    if target == "edgelist":
        for gtype in gtype_list:
            for data_providers in data_provider_list:
                for i in range(len(target_config["n"])):
                    n = target_config["n"][i]
                    e = -1
                    data_prefix = f"{gtype}_n-{n}"
                    file_prefix = f"{target_config['prefix']}{data_prefix}"
                    meta_file = f"{file_prefix}_edge_list_meta_multiparty.txt"
                    
                    # if data do not exist, generate the data.
                    if not os.path.exists(meta_file):
                        print(f"File {meta_file} does not exist. Generate the data...")
                        generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --p {data_providers} --saving_type edgelist"
                        os.system(generate_command)
                    
                    # run the profiling.
                    for j in range(REPEAT_TIMES):
                        count = j+1
                        run_args = f" -multi_edgelist true -prefix {data_prefix} -rcounter {count}"
                        print(f"Run the edgelist profiling with {run_args}")
                        os.system(f"./Eval/dis_exec.sh \"{run_args}\"")
                        os.system(f"cat ./debug.txt; rm ./debug.txt")
                        record_file = f"{target_config['record_folder']}{data_prefix}-{count}.txt"
                        target_record_file = f"{target_config['record_folder']}{data_prefix}_p-{data_providers}-{count}.txt"
                        os.system(f"mv {record_file} {target_record_file}")

