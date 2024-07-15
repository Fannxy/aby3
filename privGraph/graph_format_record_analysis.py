import argparse
import pandas as pd
import numpy as np
import json
import os
import re
# from graph_format_benchmark import format_configs, gtype_list
import graph_format_benchmark as bb
import graph_format_integration_benchmark as ib

basic_querys = ["EdgeExistQuery", "OuttingEdgesCountQuery", "NeighborsGetQuery"]
basic_query_communications = ["EdgeExistQuery_recv", "OuttingEdgesCountQuery_recv", "NeighborsGetQuery_recv", "EdgeExistQuery_send", "OuttingEdgesCountQuery_send", "NeighborsGetQuery_send"]
result_online_folder = "./privGraph/results/"
result_offline_folder = "./privGraph/results_offline/"

if(not os.path.exists(result_online_folder)):
    os.makedirs(result_offline_folder)

if(not os.path.exists(result_offline_folder)):
    os.makedirs(result_offline_folder)

analysis_dict = {
    "privGraph": {
        "gtype": [],
        "n": [],
        "e": [],
        "k": [],
        "l": [],
        "se": [],
        "pe": [],
        "sn": [],
        "pn": [],
        "GraphLoad": [],
        "EdgeOramInit": [],
        "NodeOramInit": [],
        "GraphLoad_recv": [],
        "EdgeOramInit_recv": [],
        "NodeOramInit_recv": [],
        "GraphLoad_send": [],
        "EdgeOramInit_send": [],
        "NodeOramInit_send": [],
    },
    "adjmat": {
        "gtype": [],
        "n": [],
        "e": [],
        "se": [],
        "pe": [],
        "sn": [],
        "pn": [],
        "GraphLoad": [],    
        "EdgeOramInit": [],
        "NodeOramInit": [],
        "GraphLoad_recv": [],
        "EdgeOramInit_recv": [],
        "NodeOramInit_recv": [],
        "GraphLoad_send": [],
        "EdgeOramInit_send": [],
        "NodeOramInit_send": [],
    },
    "edgelist": {
        "gtype": [],
        "n": [],
        "e": [],
        "GraphLoad": [],
        "GraphLoad_recv": [],
        "GraphLoad_send": [],
    }
}


def logging2dict(log_file, result_dict):
    with open(log_file, "r") as f:
        for line in f:
            if "Time taken" in line:
                parts = line.strip().split(" ")
                key = parts[3][:-1]
                value = float(parts[4])
                result_dict[key] = value
            elif " : " in line:
                parts = line.strip().split(" : ")
                key = parts[0]
                value = int(parts[1])
                result_dict[key] = value
            elif "Communicaitions" in line:
                parts = line.strip().split(" ")
                key = parts[3][:-1]
                # print("key = ", key)
                value = float(parts[4])
                result_dict[key] = value

    return result_dict


def get_privGraph_record(target_folder, gtype, n, k, count):
    """fetch a result from the target file.
    """
    target_file = f"{target_folder}{gtype}_n-{n}_k-{k}-{count}.txt"
    result_dict =  {"gtype": gtype, "n": n, "k": k}
    return logging2dict(target_file, result_dict)


def get_privGraph_record_offline(target_folder, gtype, n, k, party, count):
    """fetch a result from the target file.
    """
    target_file = f"{target_folder}{gtype}_n-{n}_k-{k}_p-{party}-{count}.txt"
    result_dict =  {"gtype": gtype, "n": n, "k": k, "party": party}
    return logging2dict(target_file, result_dict)


def get_baseline_record(target_folder, gtype, n, count):
    """fetch a result from the target file.
    """
    target_file = f"{target_folder}{gtype}_n-{n}-{count}.txt"
    result_dict = {"gtype": gtype, "n": n}
    
    return logging2dict(target_file, result_dict)


def get_baseline_record_offline(target_folder, gtype, n, party, count):
    """fetch a result from the target file.
    """
    target_file = f"{target_folder}{gtype}_n-{n}_p-{party}-{count}.txt"
    result_dict = {"gtype": gtype, "n": n, "party": party}
    
    return logging2dict(target_file, result_dict)


def record_dict_construct(graph_format, target="online"):
    
    record_dict = analysis_dict[graph_format]
    
    target_config = None
    if(target == "online"):
        target_config = bb.format_configs[graph_format]
    elif(target == "offline"):
        target_config = ib.format_configs[graph_format]
    
    target_folder = target_config["record_folder"]
    
    if(target == "online"):
        for filename in os.listdir(target_folder):
            match = re.match(target_config["record_pattern"], filename)
            if(match):
                if(graph_format == "privGraph"):
                    gtype = match.group(1)
                    n = int(match.group(2))
                    k = int(match.group(3))
                    count = int(match.group(4))
                    result_dict = get_privGraph_record(target_folder, gtype, n, k, count)
                else:
                    gtype = match.group(1)
                    n = int(match.group(2))
                    count = int(match.group(3))
                    result_dict = get_baseline_record(target_folder, gtype, n, count)

    elif(target == "offline"):
        for filename in os.listdir(target_folder):
            match = re.match(target_config["record_pattern"], filename)
            if(match):
                if(graph_format == "privGraph"):
                    gtype = match.group(1)
                    n = int(match.group(2))
                    k = int(match.group(3))
                    party = int(match.group(4))
                    count = int(match.group(5))
                    result_dict = get_privGraph_record_offline(target_folder, gtype, n, k, party, count)
                else:
                    gtype = match.group(1)
                    n = int(match.group(2))
                    party = int(match.group(3))
                    count = int(match.group(4))
                    result_dict = get_baseline_record_offline(target_folder, gtype, n, party, count)
        
    for key in record_dict:
        if key in result_dict:
            record_dict[key].append(result_dict[key])
        else:
            print(f"key {key} not found in {result_dict}")
        
    record_df = pd.DataFrame(record_dict)
    grouped = record_df.groupby(target_config["config_keys"])
    record_df = grouped.agg({key: [np.mean, np.std] for key in target_config["performance_keys"]})
    record_df.columns = ['_'.join(col).strip() for col in record_df.columns.values]  
    # record_df.to_excel(result_folder + f"{graph_format}_record.xlsx")

    return record_df
        
        
if __name__ == "__main__":
        
        parser = argparse.ArgumentParser()
        parser.add_argument('--target', type=str, default="privGraph", help="target graph format")
        # parser.add_argument('--online', type=bool, default=True, help="online or offline")
        parser.add_argument('--stage', type=str, default="online", help="online analysis")
        args = parser.parse_args()

        target = args.target
        online = True
        if(args.stage == "offline"):
            online = False
        
        if(online):
            for key in analysis_dict:
                analysis_dict[key].update({query: [] for query in basic_querys})
                analysis_dict[key].update({query: [] for query in basic_query_communications})
            result_folder = result_online_folder
            record_df = record_dict_construct(target, "online")
            record_df.to_excel(result_folder + f"{target}_record.xlsx")
        else:
            result_folder = result_offline_folder
            record_df = record_dict_construct(target, "offline")
            record_df.to_excel(result_folder + f"{target}_record.xlsx")
        
        print("Record analysis done.")