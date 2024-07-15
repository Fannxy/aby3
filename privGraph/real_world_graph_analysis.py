import argparse
import pandas as pd
import numpy as np
import os

MPI = False 
MPI_TASK = 4

MAIN_FOLDER = "/root/aby3/aby3-GraphQuery"
real_world_data_folder = MAIN_FOLDER + "/data/realworld/"  
main_record_folder = MAIN_FOLDER + "/record/realworld/"
mpi_data_folder = MAIN_FOLDER + "/data/realworld_mpi/"

n_stash_size, n_pack_size, e_stash_size, e_pack_size = 32, 16, 1024, 32
test_format = ["privGraph", "edgelist"]
# test_format = ["edgelist"]

def data_synchronize(data_folder):
    os.system(f"ssh aby31 'rm -rf {data_folder}*'")
    os.system(f"ssh aby32 'rm -rf {data_folder}*'")
    os.system(f"ssh aby31 'mkdir -p {data_folder}'")
    os.system(f"ssh aby32 'mkdir -p {data_folder}'")
    os.system(f"scp -r {data_folder}* aby31:{data_folder}")
    os.system(f"scp -r {data_folder}* aby32:{data_folder}")
   

REPEAT_TIMES = 1

if __name__ == "__main__":

    # prepare the file.
    os.system("cp ./frontend/main.pgp ./frontend/main.cpp; python build.py")

    # get the test configs.
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', type=str, default="", help="target real-world graph")
    parser.add_argument('--MPI', type=bool, default=MPI, help="MPI mode")
    parser.add_argument('--MPI_TASK', type=int, default=MPI_TASK, help="MPI task numbers")
    parser.add_argument('--data_folder', type=str, default=real_world_data_folder, help="origional data folder")
    args = parser.parse_args()
    
    target = args.target
    MPI_TASK = args.MPI_TASK
    real_world_data_folder = args.data_folder
    MPI = args.MPI
    if(MPI):
        target = target + f"_{MPI_TASK}"
        
    if(MPI):
        os.system("cp ./frontend/main.pgpmpi ./frontend/main.cpp; python build.py --MPI")
        # data_synchronize(real_world_data_folder)
    else:
        os.system("cp ./frontend/main.pgp ./frontend/main.cpp; python build.py")
    
    print("after data synchronization!")
    # exit(0)
    
    for gformat in test_format:
        # different graph format.
        target_record_folder = main_record_folder + f"{gformat}/"
        if(not os.path.exists(target_record_folder)):
            os.makedirs(target_record_folder)

        run_args = f" -{gformat} -prefix {target} -rcounter 1 -noram_stash_size {n_stash_size} -eoram_stash_size {e_stash_size} -noram_pack_size {n_pack_size} -eoram_pack_size {e_pack_size} -data_folder {real_world_data_folder} -record_folder {target_record_folder}"

        print(run_args)
        if(MPI):
            os.system(f"./Eval/mpi_dis_exec.sh \"{run_args}\" {MPI_TASK}")
        else:
            print(">>>>>>>>> in this branch!!!!!!")
            os.system(f"./Eval/dis_exec.sh \"{run_args}\"")
            
        os.system(f"cat ./debug.txt; rm ./debug.txt")
