import argparse
import pandas as pd
import numpy as np
import os

MPI_TASK = 4

MAIN_FOLDER = "/root/aby3/aby3-GraphQuery"
origional_data_folder = MAIN_FOLDER + "/data/realworld/" 
mpi_data_folder = MAIN_FOLDER + "/data/realworld_mpi/"


if __name__ == "__main__":
    
    # get the test configs.
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', type=str, default="", help="target real-world graph")
    parser.add_argument('--MPI_TASK', type=int, default=2, help="MPI task numbers")
    parser.add_argument('--origional_folder', type=str, default=origional_data_folder, help="origional data folder")
    parser.add_argument('--mpi_folder', type=str, default=mpi_data_folder, help="mpi data folder")
    args = parser.parse_args()
    target = args.target
    MPI_TASK = args.MPI_TASK
    origional_data_folder = args.origional_folder
    mpi_data_folder = args.mpi_folder
    
    if(not os.path.exists(mpi_data_folder)):
        os.makedirs(mpi_data_folder)
        
    # read the edgelist.
    edge_list = np.loadtxt(origional_data_folder + target + "_edge_list.txt", dtype=int)
    edge_list_property = np.loadtxt(origional_data_folder + target + "_edge_list_property.txt", dtype=int)
    print("edge_list size = ", edge_list.shape)
    target_size = np.ceil(edge_list.shape[0] / MPI_TASK).astype(int)
    print("target_size = ", target_size)
    target_edges = edge_list[:target_size]
    target_property = edge_list_property[:target_size]
    print("target_edges size = ", target_edges.shape)
    
    target_file = mpi_data_folder + target + f"_{MPI_TASK}_edge_list.txt"
    target_meta_file = mpi_data_folder + target + f"_{MPI_TASK}_edge_list_meta.txt"
    target_property_file = mpi_data_folder + target + f"_{MPI_TASK}_edge_list_property.txt"
    # save the edge file
    np.savetxt(target_file, target_edges, fmt="%d", delimiter=' ')
    np.savetxt(target_property_file, target_property, fmt="%d", delimiter=' ')
    # save the meta file.
    with open(target_meta_file, "w") as f:
        f.write(f"{target_edges.shape[0]} {target_edges.shape[1]}")
    
    # read the 2d-partition data.
    partition_data = np.loadtxt(origional_data_folder + target + "_2dpartition.txt", dtype=int)
    partition_meta_file = origional_data_folder + target + "_meta.txt"
    with open(partition_meta_file, "r") as f:
        partition_meta = f.readline().strip().split(" ")
        v, e, b, k, l = map(int, partition_meta)
    
    print("partition_data size = ", partition_data.shape)
    print("expected size = ", (b**2 * l))
    target_l = np.ceil(l / MPI_TASK).astype(int)
    print("target_l = ", target_l)
    target_expected_size = b**2 * target_l
    target_2d_partition = partition_data[:target_expected_size]
    print("target_2d_partition size = ", target_2d_partition.shape)
    target_partition_meta_file = mpi_data_folder + target + f"_{MPI_TASK}_meta.txt"
    partition_property_data = np.loadtxt(origional_data_folder + target + "_2dproperty.txt", dtype=int)
    target_partition_property_file = mpi_data_folder + target + f"_{MPI_TASK}_2dproperty.txt"
    target_partition_property = partition_property_data[:target_expected_size]
    
    np.savetxt(mpi_data_folder + target + f"_{MPI_TASK}_2dpartition.txt", target_2d_partition, fmt="%d", delimiter=' ')
    np.savetxt(target_partition_property_file, target_partition_property, fmt="%d", delimiter=' ')
    meta_dict = {"v": v, "e": e, "b": b, "k": k, "l": target_l}
    with open(target_partition_meta_file, "w") as f:
        f.write(f"{v} {e} {b} {k} {target_l}")
    
    
    