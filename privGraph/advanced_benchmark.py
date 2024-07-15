"""Run advanced applications.
"""
import argparse
import json
import os
import random
from utils import get_k 

MPI = True 
MPI_TASK = 4

MAIN_FOLDER = "/root/aby3/aby3-GraphQuery"
record_folder = MAIN_FOLDER + "/record/adv_application/"
data_folder = ""
data_prefix = ""

n_stash_size, n_pack_size, e_stash_size, e_pack_size = 32, 16, 1024, 32

data_type = "real_world"

party1 = "aby31"
party2 = "aby32"
party_list = [party1, party2]

data_prefix_list = []

if(data_type == "synthetic"):
    data_folder = MAIN_FOLDER + "/data/baseline/"
    data_prefix_list = ["bipartite_n-1024_k-4"]
else:
    data_folder = MAIN_FOLDER + "/data/realworld/"
    # data_prefix_list = ["slashdot", "dblp"]
    data_prefix_list = ["twitter"]
    
    
if __name__ == "__main__":
    
    # get the test configs.
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=str, default="cycle_detect", help="The target advanced applications")
    parser.add_argument('--MPI', type=bool, default=MPI, help="MPI mode")
    parser.add_argument('--MPI_TASK', type=int, default=MPI_TASK, help="MPI task numbers")
    parser.add_argument('--data_folder', type=str, default=data_folder, help="origional data folder")
    parser.add_argument('--record_folder', type=str, default=record_folder, help="record folder")
    args = parser.parse_args()
    
    data_folder = args.data_folder
    record_folder = args.record_folder
    MPI = args.MPI
    MPI_TASK = args.MPI_TASK
    
    if(MPI):
        os.system("cp ./frontend/main.pgpmpi ./frontend/main.cpp; python build.py --MPI")
    else:
        os.system("cp ./frontend/main.pgp ./frontend/main.cpp; python build.py")
    
    if(not os.path.exists(data_folder)):
        os.makedirs(data_folder)
    if(not os.path.exists(record_folder)):
        os.makedirs(record_folder)
    
    for data_prefix in data_prefix_list:
        
        if(MPI):
            data_prefix = data_prefix + f"_{MPI_TASK}"
        # if(args.target == "neighbor_stats"):
        #     # generate the property_graph.
        #     meta_file = data_folder + data_prefix + "_meta.txt"
        #     print(meta_file)
        #     with open(meta_file, "r") as f:
        #         line = f.readline()
        #         v, e, b, k, l = map(int, line.split())
        #     print(f"v: {v}, e: {e}, b: {b}, k: {k}, l: {l}")
        #     # generate the random property graph.
        #     property_file = data_folder + data_prefix + "_2dproperty.txt"
        #     with open(property_file, "w") as f:
        #         for i in range(b*b):
        #             for j in range(l):
        #                 random_number = random.randint(1, 10000)
        #                 f.write(f"{random_number}\n")
                    
        #     edgelist_meta_file = data_folder + data_prefix + "_edge_list_meta.txt"
        #     with open(edgelist_meta_file, "r") as f:
        #         line = f.readline()
        #         v, e = map(int, line.split())
        #     print(f"v: {v}, e: {e}")
            
        #     # generate the random edgelist.
        #     property_edgelist_file = data_folder + data_prefix + "_edge_list_property.txt"
        #     with open(property_edgelist_file, "w") as f:
        #         for i in range(e):
        #             random_number = random.randint(1, 10000)
        #             f.write(f"{random_number}\n")
            
            # synchronize data with other parties.
            # data synchronize
            # for party in party_list:
            #     # generate data folders on the remote party machine, also, recursive mkdir for target real_world_data_folder
            #     os.system(f"ssh {party} 'rm -rf {data_folder}*'")
            #     os.system(f"ssh {party} 'mkdir -p {data_folder}'")
            #     os.system(f"scp -r {data_folder}* {party}:{data_folder}")


        # set the graph.
        run_args = f" -{args.target} -prefix {data_prefix} -rcounter 1 -noram_stash_size {n_stash_size} -eoram_stash_size {e_stash_size} -noram_pack_size {n_pack_size} -eoram_pack_size {e_pack_size} -data_folder {data_folder} -record_folder {record_folder}"
        
        print(run_args)
        print(f"Run {args.target}...")
        if(MPI):
            os.system(f"./Eval/mpi_dis_exec.sh \"{run_args}\" {MPI_TASK}")
        else:
            os.system(f"./Eval/dis_exec.sh \"{run_args}\"")
            
        os.system(f"cat ./debug.txt; rm ./debug.txt")
        
