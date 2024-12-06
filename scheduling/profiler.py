import sys
import os
import argparse
import threading
import re
import math
import pandas as pd
import json
import subprocess
import shlex
from task_assigning import assign_task
import copy
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PtA_deploy')))
from system_monitor import *

root_folder = "/root/aby3/"

def calculate_expression(n, expression):
    allowed_functions = {name: getattr(math, name) for name in dir(math) if not name.startswith("__")}
    allowed_functions['n'] = n
    result = eval(expression, allowed_functions)
    return result

def get_bandwidth(i, time, server_host, ip_address, test_server, parallel=1):
    os.system(f"ssh {server_host} 'iperf3 -s -D'")
    result = os.popen(f"ssh {test_server} 'iperf3 -c {ip_address} -t {time} -P {parallel} -J'").read()
    data = json.loads(result)
    bandwidth = data["end"]["sum_received"]["bits_per_second"] / (2**30)
    os.system(f"ssh {server_host} 'pkill iperf3'")
    return bandwidth

def run_command(command):
    os.system(command)

def analysis(server_host, keyword, command, record_folder, interface):
    os.system(f"ssh {server_host} python {root_folder}/scheduling/monitor_dis_run_profiler.py --keyword {keyword} --command ' {command.replace('-', '+')} ' --record_folder {record_folder} --interface {interface}")

def collect_network_usage(data_size, args, record_folder, keyword, server_host, ip_address, network_interface):
    threads = []
    for i in range(3):
        command = f"{root_folder}out/build/linux/frontend/frontend -dataSize {data_size} -role {i} {args} -p0_ip {ip_address[0]} -p1_ip {ip_address[1]} -rank 0"
        print(command)
        thread = threading.Thread(target=analysis, args=(server_host[i], f"{keyword}-{data_size}", command, record_folder, network_interface[i]))
        threads.append(thread)
    
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    
    os.system(f"mv {record_folder}/monitor-{keyword}-{data_size}.log {record_folder}/monitor-{keyword}-{data_size}-0.log")
    for i in range(1, n):
        os.system(f"scp -r {server_host[i]}:{record_folder}/monitor-{keyword}-{data_size}.log {record_folder}/monitor-{keyword}-{data_size}-{i}.log")
    return

def get_profile_usage_dict(data_size, role, args, record_folder, keyword, server_host, ip_address, network_interface):
    if not os.path.exists(f"{record_folder}/monitor-{keyword}-{data_size}-{role}.log"):
        collect_network_usage(data_size, args, record_folder, keyword, server_host, ip_address, network_interface)
    return get_usage_dict(f"{record_folder}/monitor-{keyword}-{data_size}-{role}.log")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--record_folder', type=str, default=root_folder+"scheduling/Record_test/", help='record folder')
    parser.add_argument('--config_folder', type=str, default=root_folder+"scheduling/config/", help='config folder')

    parser.add_argument('--args', type=str, help='run args')
    parser.add_argument('--args_agg', type=str, help='run args for aggregation')
    parser.add_argument('--keyword', type=str, help='keyword')
    parser.add_argument('--task', type=str, help='task name')

    parser.add_argument('--num_parties', type=int, default=3, help='number of parties')
    parser.add_argument('--server_host', type=str, nargs='+', default=["aby30", "aby31", "aby32"], help='server host')
    parser.add_argument('--ip_address', type=str, nargs='+', default=["10.3.0.13", "10.3.0.16", "10.3.0.17"], help='ip address')
    parser.add_argument('--network_interface', type=str, nargs='+', default=["ibs110", "ibs110", "ibs110"], help='network interface')
    parser.add_argument('--profile_ip_address', type=str, nargs='+', default=["10.3.0.13", "10.3.0.16", "10.3.0.17"], help='ip address')
    parser.add_argument('--profile_network_interface', type=str, nargs='+', default=["ibs110", "ibs110", "ibs110"], help='network interface')

    parser.add_argument('--data_size', type=int, help='data size')
    parser.add_argument('--get_bandwidth_time', type=int, default=5, help='get bandwidth time')
    parser.add_argument('--fitting_length', type=int, default=10, help='fitting length')
    parser.add_argument('--fitting_step', type=int, default=100, help='fitting step')
    parser.add_argument('--complexity', type=str, nargs='+', default=['1', 'n'], help='communication complexity of each stage')
    parser.add_argument('--parallelism_limit', type=int, default=64, help='parallelism limit')
    parser.add_argument('--fix_parallelism', action='store_true', help='fix parallelism')
    parser.add_argument('--run_tasks', action='store_true', help='run tasks')
    parser.add_argument('--MPI', action='store_true', help='run in MPI')
    parser.add_argument('--assignment_strategy', type=str, choices=['baseline', 'roundrole'], default='roundrole', help='baseline or roundrole')

    args = parser.parse_args()
    n = args.num_parties
    data_size = args.data_size
    
    # mkdir the record folder.
    if not os.path.exists(args.record_folder):
        os.makedirs(args.record_folder)
    
    if not os.path.exists(args.config_folder):
        os.makedirs(args.config_folder)
    
    # get bandwidth
    print("Getting bandwidth")
    bandwidth = []
    # if the bandwidth file exists, load the bandwidth from the file.
    if os.path.exists(f"{args.config_folder}/bandwidth.json"):
        with open(f"{args.config_folder}/bandwidth.json", "r") as f:
            network_dict = json.load(f)
        for i in range(n):
            bandwidth.append(network_dict[args.server_host[i]])
    else: # otherwise get the bandwidth from the servers.
        for i in range(n):
            test_server = args.server_host[i]
            for j in range(n):
                if args.server_host[j] != args.server_host[i]:
                    test_server = args.server_host[j]
                    break
            bandwidth_i = []
            for _ in range(5):
                bandwidth_i.append(get_bandwidth(i, args.get_bandwidth_time, args.server_host[i], args.ip_address[i], test_server))
            bandwidth.append(max(bandwidth_i))
        for i in range(n):
            print(f"{args.server_host[i]}: {bandwidth[i]}Gb/s")
        network_dict = {args.server_host[i]: bandwidth[i] for i in range(n)}
        with open(f"{args.config_folder}/bandwidth.json", "w") as f:
            json.dump(network_dict, f)

    # collect network usage
    print("Collecting network usage")
    length = args.fitting_length
    step = args.fitting_step
    usage_dict = [{} for _ in range(3)]
    usage_dict_agg = [{} for _ in range(3)]
    for size in range(step, step * length + 1, step):
        for role in range(3):
            print(f"Collecting network usage: {size} {role}")
            usage_dict[role][size] = get_profile_usage_dict(
                data_size=size,
                role=role,
                args=args.args,
                record_folder=args.record_folder,
                keyword=args.task,
                server_host=args.server_host,
                ip_address=args.profile_ip_address,
                network_interface=args.profile_network_interface
            )
    
    recv_mean_usage = [{} for _ in range(3)]
    send_mean_usage = [{} for _ in range(3)]
    agg_time = {}
    size = step
    while size < data_size:
        max_mean_usage = 0
        for role in range(n):
            usage_dict[role][size] = get_profile_usage_dict(
                data_size=size,
                role=role,
                args=args.args,
                record_folder=args.record_folder,
                keyword=args.task,
                server_host=args.server_host,
                ip_address=args.profile_ip_address,
                network_interface=args.profile_network_interface
            )
            usage_dict_agg[role][size] = get_profile_usage_dict(
                data_size=size,
                role=role,
                args=args.args_agg,
                record_folder=args.record_folder,
                keyword=f"{args.task}-agg",
                server_host=args.server_host,
                ip_address=args.profile_ip_address,
                network_interface=args.profile_network_interface
            )
            recv_mean_usage[role][size] = np.mean(usage_dict[role][size]["network_recv"])
            send_mean_usage[role][size] = np.mean(usage_dict[role][size]["network_send"])
            max_mean_usage = max(max_mean_usage, recv_mean_usage[role][size], send_mean_usage[role][size])
        agg_time[size] = np.mean([len(usage_dict_agg[role][size]["network_recv"]) for role in range(n)])
        if max_mean_usage > max(bandwidth):
            break
        size *= 2
    
    # for size in agg_time.keys():
    #     comp_time = np.mean([len(usage_dict[role][size]["network_recv"]) for role in range(n)])
    #     print(f"size: {size}, agg_time: {agg_time[size]}, comp_time: {comp_time}")
    
    # fit network usage with given expression
    print("Fitting network usage")
    coef_recv = []
    coef_send = []
    for role in range(3):
        sizes = range(step, step * length + 1, step)
        sum_recv = []
        sum_send = []
        for size in sizes:
            network_recv, network_send = usage_dict[role][size]["network_recv"], usage_dict[role][size]["network_send"]
            sum_recv.append(np.sum(network_recv))
            sum_send.append(np.sum(network_send))
        complexity_matrix = np.array([
            [calculate_expression(size, complexity_expr) for complexity_expr in args.complexity]
            for size in sizes
        ])
        coef_recv.append(np.linalg.lstsq(complexity_matrix, sum_recv, rcond=None)[0])
        coef_send.append(np.linalg.lstsq(complexity_matrix, sum_send, rcond=None)[0])
    expr_recv = ["0"] * n
    expr_send = ["0"] * n
    for i in range(3):
        expr_recv[i] = " + ".join(f"({coef}) * ({complexity_expr})" for coef, complexity_expr in zip(coef_recv[i], args.complexity))
        expr_send[i] = " + ".join(f"({coef}) * ({complexity_expr})" for coef, complexity_expr in zip(coef_send[i], args.complexity))
    print(expr_recv)
    print(expr_send)

    # assign tasks
    time_stamp_file_prefix = f"{args.record_folder}/stamp-{args.keyword}-{data_size}"

    print("Assigning tasks")
    comp_assignment, aggr_assignment = assign_task(strategy=args.assignment_strategy, bandwidth=bandwidth, expr_recv=expr_recv, expr_send=expr_send, recv_mean_usage=recv_mean_usage, send_mean_usage=send_mean_usage, data_size=data_size, parallelism_limit=args.parallelism_limit, agg_time=agg_time, parallelism_min=(args.parallelism_limit if args.fix_parallelism else 1))
    print("Task assignment:", comp_assignment)
    print("Aggregation assignment:", aggr_assignment)
    coef_recv = [calculate_expression(data_size + 1, expr) - calculate_expression(data_size, expr) for expr in expr_recv]
    coef_send = [calculate_expression(data_size + 1, expr) - calculate_expression(data_size, expr) for expr in expr_send]
    max_coef = max(max(coef_recv), max(coef_send))
    with open(f"{args.config_folder}/task_assignment-{args.assignment_strategy}.txt", "a") as f:
        for i in range(3):
            f.write(f"role {i} recv: %.2f\n" % (coef_recv[i] / max_coef))
            f.write(f"role {i} send: %.2f\n" % (coef_send[i] / max_coef))
        f.write(f"parallelism: {len(comp_assignment)}\n")
        json.dump(comp_assignment, f)
        json.dump(aggr_assignment, f)
        f.write(f"\ntask: {args.keyword}\n")
    
    time_stamp_file_prefix = f"{args.record_folder}/stamp-{args.keyword}-{data_size}-{args.assignment_strategy}"

    # run the tasks
    if args.run_tasks:
        print("Running tasks")

        commands = [[] for i in range(n)]
        for rank, param in enumerate(comp_assignment):
            role_assignment = param[0]
            subtask_size = param[1]
            party_id = []
            for role in range(3):
                for i in range(n):
                    if role_assignment[i] == role:
                        party_id.append(i)
                        break
            for role in range(3):
                p0_ip = args.ip_address[party_id[0]]
                p1_ip = args.ip_address[party_id[1]]
                command = f"{root_folder}out/build/linux/frontend/frontend -dataSize {subtask_size} -role {role} {args.args} -rank {rank} -p0_ip {p0_ip} -p1_ip {p1_ip} -stampFile {time_stamp_file_prefix}-{rank}.txt"
                commands[party_id[role]].append(command)
        
        commands_agg = []
        for agg_layer in aggr_assignment:
            commands_agg_layer = [[] for i in range(n)]
            for rank, param in enumerate(agg_layer):
                role_assignment = param[0]
                subtask_size = param[1]
                party_id = []
                for role in range(3):
                    for i in range(n):
                        if role_assignment[i] == role:
                            party_id.append(i)
                            break
                for role in range(3):
                    p0_ip = args.ip_address[party_id[0]]
                    p1_ip = args.ip_address[party_id[1]]
                    command = f"{root_folder}out/build/linux/frontend/frontend -dataSize {subtask_size} -role {role} {args.args_agg} -rank {rank} -p0_ip {p0_ip} -p1_ip {p1_ip}"
                    commands_agg_layer[party_id[role]].append(command)
            commands_agg.append(commands_agg_layer)
        
        total_command = []
        for node in range(n):
            if not args.MPI:
                command = " \& ".join(commands[node]) + " \& wait"
                for agg_layer in commands_agg:
                    command += " \; "
                    command += " \& ".join(agg_layer[node]) + " \& wait"
                total_command.append(command)
            else:
                mpi_command = " mpirun"
                for i, command in enumerate(commands[node]):
                    mpi_command += " -np 1 " + command
                    if i != len(commands[node]) - 1:
                        mpi_command += " :"
                for agg_layer in commands_agg:
                    mpi_command += " \; mpirun"
                    for i, command in enumerate(agg_layer[node]):
                        mpi_command += " -np 1 " + command
                        if i != len(agg_layer[node]) - 1:
                            mpi_command += " :"
                total_command.append(mpi_command)

        threads = []
        for i in range(n):
            thread = threading.Thread(target=analysis, args=(args.server_host[i], f"{args.keyword}-{data_size}", total_command[i], args.record_folder, args.network_interface[i]))
            threads.append(thread)

        start_time = time.time()
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        end_time = time.time()
        
        # concate all the time stamp files.
        os.system(f"cat {time_stamp_file_prefix}-*.txt > {time_stamp_file_prefix}.txt")
        os.system(f"rm {time_stamp_file_prefix}-*.txt")
        
        os.system(f"mv {time_stamp_file_prefix}.txt {time_stamp_file_prefix}-0.txt")
        os.system(f"mv {args.record_folder}/monitor-{args.keyword}-{data_size}.log {args.record_folder}/monitor-{args.keyword}-{data_size}-0.log")
        for i in range(1, n):
            os.system(f"ssh {args.server_host[i]} \"cat {time_stamp_file_prefix}-*.txt > {time_stamp_file_prefix}.txt\"")
            os.system(f"ssh {args.server_host[i]} \"rm {time_stamp_file_prefix}-*.txt\"")
            os.system(f"scp -r {args.server_host[i]}:{args.record_folder}/monitor-{args.keyword}-{data_size}.log {args.record_folder}/monitor-{args.keyword}-{data_size}-{i}.log")
            os.system(f"scp -r {args.server_host[i]}:{time_stamp_file_prefix}.txt {time_stamp_file_prefix}-{i}.txt")

        keyword = f"{args.keyword}"
        recv_utilization = []
        send_utilization = []
        total_utilization = 0
        total_bandwidth = 0
        for i in range(n):
            usage_dict = get_usage_dict(f"{args.record_folder}/monitor-{args.keyword}-{data_size}-{i}.log")
            time_stamp_dict = get_time_stamp(f"{time_stamp_file_prefix}-{i}.txt")
            if not time_stamp_dict:
                draw_usage_graph(usage_dict, f"{args.record_folder}/{keyword}-{data_size}-{i}.png")
            else:
                draw_usage_graph(usage_dict, f"{args.record_folder}/{keyword}-{data_size}-{i}.png", time_stamp_dict)
            recv_u = np.mean(usage_dict["network_recv"]) / bandwidth[i]
            send_u = np.mean(usage_dict["network_send"]) / bandwidth[i]
            total_utilization += np.mean(usage_dict["network_recv"]) + np.mean(usage_dict["network_send"])
            total_bandwidth += bandwidth[i]
            recv_utilization.append(recv_u)
            send_utilization.append(send_u)
        total_utilization /= total_bandwidth*2

        file_path = f"{args.record_folder}/record.xlsx"
        if not os.path.exists(f"{args.record_folder}/record.xlsx"):
            df = pd.DataFrame()
            df.to_excel(file_path, index=False)
        
        df = pd.read_excel(file_path, engine='openpyxl')
        new_row = {}
        new_row['keyword'] = [f"{keyword}"]
        new_row['data size'] = [data_size]
        new_row['time'] = [end_time - start_time]
        new_row["parallelism"] = [len(comp_assignment)]
        for i in range(n):
            new_row[f"recv utilization-{i}"] = [recv_utilization[i]]
            new_row[f"send utilization-{i}"] = [send_utilization[i]]
        new_row[f"total utilization-{i}"] = [total_utilization]
        df = pd.concat([df, pd.DataFrame(new_row)], ignore_index=True)
        df.to_excel(file_path, index=False)