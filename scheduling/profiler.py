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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PtA_deploy')))
from system_monitor import *

root_folder = "/root/aby3/"

def calculate_expression(n, expression):
    allowed_functions = {name: getattr(math, name) for name in dir(math) if not name.startswith("__")}
    allowed_functions['n'] = n
    result = eval(expression, allowed_functions)
    return result

def get_bandwidth(i, time, server_host, ip_address, test_server):
    os.system(f"ssh {server_host} 'iperf3 -s -D'")
    result = os.popen(f"ssh {test_server} 'iperf3 -c {ip_address} -t {time} -J'").read()
    data = json.loads(result)
    bandwidth = data["end"]["sum_received"]["bits_per_second"] / (2**30)
    os.system(f"ssh {server_host} 'pkill iperf3'")
    return bandwidth

def run_command(command):
    os.system(command)

def collect_network_usage(data_size, role, args, record_folder, keyword, server_host, ip_address, network_interface):
    monitor = SystemMonitor(0.01)
    monitor.start_all(interface=network_interface[0])

    threads = []
    role_assignment = [role, (role + 1) % 3, (role + 2) % 3]
    index = list(range(3))
    for i in range(3):
        index[role_assignment[i]] = i
    for i in range(3):
        command = f"{root_folder}out/build/linux/frontend/frontend -dataSize {data_size} -role {role_assignment[i]} {args} -p0_ip {ip_address[index[0]]} -p1_ip {ip_address[index[1]]} -rank 0"
        if i != 0:
            command = f"ssh {server_host[i]} " + command
        print(command)
        thread = threading.Thread(target=run_command, args=(command,))
        threads.append(thread)
        thread.start()
    
    for thread in threads:
        thread.join()

    monitor.stop_and_output(f"{record_folder}/monitor-{keyword}.log")
    usage_dict = get_usage_dict(f"{record_folder}/monitor-{keyword}.log")
    return usage_dict

def analysis(server_host, keyword, command, record_folder, interface):
    os.system(f"ssh {server_host} python {root_folder}/scheduling/monitor_dis_run_profiler.py --keyword {keyword} --command ' {command.replace('-', '+')} ' --record_folder {record_folder} --interface {interface}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--args', type=str, help='run args')
    parser.add_argument('--record_folder', type=str, help='record folder')
    parser.add_argument('--keyword', type=str, help='keyword')
    parser.add_argument('--num_parties', type=int, default=3, help='number of parties')
    parser.add_argument('--server_host', type=str, nargs='+', default=["aby30", "aby31", "aby32"], help='server host')
    parser.add_argument('--ip_address', type=str, nargs='+', default=["10.3.0.15", "10.3.0.16", "10.3.0.17"], help='ip address')
    parser.add_argument('--network_interface', type=str, nargs='+', default=["ibs110", "ibs110", "ibs110"], help='network interface')
    parser.add_argument('--data_size', type=int, help='data size')
    parser.add_argument('--get_bandwidth_time', type=int, default=5, help='get bandwidth time')
    parser.add_argument('--skip_monitor', action='store_true', help='skip monitor')
    parser.add_argument('--fitting_length', type=int, default=10, help='fitting length')
    parser.add_argument('--fitting_step', type=int, default=100, help='fitting step')
    parser.add_argument('--complexity', type=str, nargs='+', default=['1', 'n'], help='communication complexity of each stage')
    parser.add_argument('--parallelism_limit', type=int, default=64, help='parallelism limit')
    parser.add_argument('--run_tasks', action='store_true', help='run tasks')
    parser.add_argument('--MPI', action='store_true', help='run in MPI')
    parser.add_argument('--baseline', action='store_true', help='run baseline')

    args = parser.parse_args()
    n = args.num_parties
    data_size = args.data_size
    
    # mkdir the record folder.
    if not os.path.exists(args.record_folder):
        os.makedirs(args.record_folder)
    
    # get bandwidth
    print("Getting bandwidth")
    bandwidth = []
    for i in range(n):
        test_server = args.server_host[i]
        for j in range(n):
            if args.server_host[j] != args.server_host[i]:
                test_server = args.server_host[j]
                break
        bandwidth.append(get_bandwidth(i, args.get_bandwidth_time, args.server_host[i], args.ip_address[i], test_server))
    for i in range(n):
        print(f"{args.server_host[i]}: {bandwidth[i]}Gb/s")
    
    # collect network usage
    print("Collecting network usage")
    length = args.fitting_length
    step = args.fitting_step
    usage_dict = [{} for _ in range(3)]
    if not args.skip_monitor:
        for role in range(3):
            sizes = range(step, step * length + 1, step)
            for size in sizes:
                usage_dict[role][size] = collect_network_usage(size, role, args.args, args.record_folder, f"{args.keyword}-{role}-{size}", args.server_host, args.ip_address, args.network_interface)
                time.sleep(0.5)
    else:
        for role in range(3):
            for size in range(step, step * length + 1, step):
                usage_dict[role][size] = get_usage_dict(f"{args.record_folder}/monitor-{args.keyword}-{role}-{size}.log")
    
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

    # decide parallelism
    print("Deciding parallelism")
    parallelism = args.parallelism_limit // 2
    min_bandwidth = min(bandwidth)
    min_bandwidth_index = bandwidth.index(min_bandwidth)
    while parallelism > 2:
        if not args.skip_monitor:
            usage_dict = collect_network_usage(data_size // parallelism, min_bandwidth_index, args.args, args.record_folder, f"{args.keyword}-{min_bandwidth_index}-{data_size}", args.server_host, args.ip_address, args.network_interface)
        else:
            usage_dict = get_usage_dict(f"{args.record_folder}/monitor-{args.keyword}-{min_bandwidth_index}-{data_size}.log")
        network_recv, network_send = usage_dict["network_recv"], usage_dict["network_send"]
        network_usage = [max(recv, send) for recv, send in zip(network_recv, network_send)]
        tmp = [np.mean(network_usage[i:i+64]) for i in range(0, len(network_usage), 64)]
        tmp = [x for x in tmp if x > min_bandwidth / 64]
        mean_usage = np.mean(tmp)
        if mean_usage * parallelism < min_bandwidth:
            break
        parallelism //= 2
    parallelism *= 2
    print(f"parallelism: {parallelism}")


    # assign tasks
    
    if args.baseline:
        res = [[[i for i in range(n)], data_size // parallelism] for i in range(parallelism)]
    else:
        res = assign_task(bandwidth=bandwidth, expr_recv=expr_recv, expr_send=expr_send, parallelism=parallelism, data_size=data_size)
    print("Task assignment:", res)

    # run the tasks
    if args.run_tasks:
        print("Running tasks")

        commands = [[] for i in range(n)]
        for rank, param in enumerate(res):
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
                command = f"{root_folder}out/build/linux/frontend/frontend -dataSize {subtask_size} -role {role} {args.args} -rank {rank} -p0_ip {p0_ip} -p1_ip {p1_ip}"
                commands[party_id[role]].append(command)
        
        total_command = []
        for node in range(n):
            if not args.MPI:
                command = " \& ".join(commands[node]) + " \& wait"
                total_command.append(command)
            else:
                mpi_command = " mpirun"
                for i, command in enumerate(commands[node]):
                    mpi_command += " -np 1 " + command
                    if i != len(commands[node]) - 1:
                        mpi_command += " :"
                total_command.append(mpi_command)

        threads = []
        for i in range(n):
            thread = threading.Thread(target=analysis, args=(args.server_host[i], f"{args.keyword}-{data_size}", total_command[i], args.record_folder, args.network_interface[i]))
            threads.append(thread)

        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        
        os.system(f"mv {args.record_folder}/monitor-{args.keyword}-{data_size}.log {args.record_folder}/monitor-{args.keyword}-{data_size}-0.log")
        for i in range(1, n):
            os.system(f"scp -r {args.server_host[i]}:{args.record_folder}/monitor-{args.keyword}-{data_size}.log {args.record_folder}/monitor-{args.keyword}-{data_size}-{i}.log")

        for i in range(n):
            usage_dict = get_usage_dict(f"{args.record_folder}/monitor-{args.keyword}-{data_size}-{i}.log")
            keyword = f"{args.keyword}"
            if args.baseline:
                keyword += "-baseline"
            keyword += f"-{i}"
            draw_usage_graph(usage_dict, f"{args.record_folder}/{keyword}-{data_size}.png")
            file_path = f"{args.record_folder}/record.xlsx"
            if not os.path.exists(f"{args.record_folder}/record.xlsx"):
                df = pd.DataFrame(columns=['keyword', 'data size', 'time stamps', 'utilization'])
                df.to_excel(file_path, index=False)
            df = pd.read_excel(file_path, engine='openpyxl')
            utilization = (np.mean(usage_dict["network_recv"]) + np.mean(usage_dict["network_send"])) / (2 * bandwidth[i])
            new_row = pd.DataFrame({'keyword': [f"{keyword}"], 'data size': [data_size], 'time stamps': [len(usage_dict["network_recv"])], 'utilization': [utilization]})
            df = pd.concat([df, new_row], ignore_index=True)
            df.to_excel(file_path, index=False)