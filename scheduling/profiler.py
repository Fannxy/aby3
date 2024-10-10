import sys
import os
import argparse
import threading
import re
import math
from task_assigning import assign_task
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PtA_deploy')))
from system_monitor import *

root_folder = "/root/aby3/"
NETWORK_INTERFACE = "ibs110"

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

def collect_network_usage(data_size, role, args, record_folder, keyword, server_host, ip_address):
    monitor = SystemMonitor(0.01)
    monitor.start_all(interface=NETWORK_INTERFACE)

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

def get_network_usage(record_folder, keyword):
    usage_dict = get_usage_dict(f"{record_folder}/monitor-{keyword}.log")
    network_recv = usage_dict["network_recv"]
    network_send = usage_dict["network_send"]
    return network_recv, network_send

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--args', type=str, help='run args')
    parser.add_argument('--record_folder', type=str, help='record folder')
    parser.add_argument('--keyword', type=str, help='keyword')
    parser.add_argument('--num_parties', type=int, default=3, help='number of parties')
    parser.add_argument('--server_host', type=str, nargs='+', default=["aby30", "aby31", "aby32"], help='server host')
    parser.add_argument('--ip_address', type=str, nargs='+', default=["10.3.0.15", "10.3.0.16", "10.3.0.17"], help='ip address')
    parser.add_argument('--data_size', type=int, help='data size')
    parser.add_argument('--get_bandwidth_time', type=int, default=5, help='get bandwidth time')
    parser.add_argument('--skip_monitor', action='store_true', help='skip monitor')
    parser.add_argument('--fitting_length', type=int, default=10, help='fitting length')
    parser.add_argument('--fitting_step', type=int, default=100, help='fitting step')
    parser.add_argument('--complexity', type=str, nargs='+', default=['1', 'n'], help='communication complexity of each stage')
    parser.add_argument('--parallelism_limit', type=int, default=64, help='parallelism limit')
    parser.add_argument('--run_tasks', action='store_true', help='run tasks')

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
        bandwidth.append(get_bandwidth(i, args.get_bandwidth_time, args.server_host[i], args.ip_address[i], args.server_host[(i+1)%n]))
    for i in range(n):
        print(f"{args.server_host[i]}: {bandwidth[i]}Gb/s")
    
    # collect network usage
    print("Collecting network usage")
    length = args.fitting_length
    step = args.fitting_step
    if not args.skip_monitor:
        for role in range(3):
            sizes = range(step, step * length + 1, step)
            for size in sizes:
                collect_network_usage(size, role, args.args, args.record_folder, f"{args.keyword}-{role}-{size}", args.server_host, args.ip_address)
                time.sleep(0.5)
    
    # fit network usage with given expression
    print("Fitting network usage")
    coef_recv = []
    coef_send = []
    for role in range(3):
        sizes = range(step, step * length + 1, step)
        mean_recv = []
        mean_send = []
        for size in sizes:
            network_recv, network_send = get_network_usage(args.record_folder, f"{args.keyword}-{role}-{size}")
            mean_recv.append(np.mean(network_recv))
            mean_send.append(np.mean(network_send))

        complexity_matrix = np.array([
            [calculate_expression(size, complexity_expr) for complexity_expr in args.complexity]
            for size in sizes
        ])

        coef_recv.append(np.linalg.lstsq(complexity_matrix, mean_recv, rcond=None)[0])
        coef_send.append(np.linalg.lstsq(complexity_matrix, mean_send, rcond=None)[0])

    expr_recv = ["0"] * n
    expr_send = ["0"] * n

    for i in range(3):
        expr_recv[i] = " + ".join(f"({coef}) * ({complexity_expr})" for coef, complexity_expr in zip(coef_recv[i], args.complexity))
        expr_send[i] = " + ".join(f"({coef}) * ({complexity_expr})" for coef, complexity_expr in zip(coef_send[i], args.complexity))
    
    print(expr_recv)
    print(expr_send)

    # assign tasks
    
    res = assign_task(bandwidth, expr_recv, expr_send, data_size)
    print("Task assignment:", res)

    # run the tasks
    if args.run_tasks:
        print("Running tasks")
        monitor = SystemMonitor(0.01)
        monitor.start_all(interface=NETWORK_INTERFACE)

        threads = []
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
        
        threads = []
        for node in range(n):
            command = " & ".join(commands[node]) + " & wait"
            if node != 0:
                command = f"ssh {args.server_host[node]} '{command}'"
            print(command)
            thread = threading.Thread(target=run_command, args=(command,))
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join()
        
        monitor.stop_and_output(f"{args.record_folder}/monitor-{args.keyword}-{data_size}.log")
        usage_dict = get_usage_dict(f"{args.record_folder}/monitor-{args.keyword}-{data_size}.log")
        draw_usage_graph(usage_dict, f"{args.record_folder}/{args.keyword}-{data_size}.png")