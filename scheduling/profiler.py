import sys
import os
import argparse
import threading
import re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PtA_deploy')))
from system_monitor import *

root_folder = "/root/aby3/"
SERVER_HOST = ["aby30", "aby31", "aby32"]
HOSTNAME = ["10.1.0.15", "10.1.0.16", "10.1.0.17"]
NETWORK_INTERFACE = ["ibs110", "ibs110", "ibs110"]
IP_ADDRESS = ["10.3.0.15", "10.3.0.16", "10.3.0.17"]

# def get_bandwidth(server_host, ip_address, time):
#     os.system(f"ssh {server_host} 'iperf3 -s -D'")
#     result = os.popen(f"ssh {server_host} 'iperf3 -c {ip_address} -t {time} -J'").read()
#     data = json.loads(result)
    
#     print("test data: ", data)
    
#     bandwidth = data["end"]["sum_received"]["bits_per_second"] / (2**30)
#     os.system(f"ssh {server_host} 'pkill iperf3'")
#     return bandwidth

def get_bandwidth(i, time):
    server_host = SERVER_HOST[i]
    ip_address = IP_ADDRESS[i]
    test_server = SERVER_HOST[(i+1)%3] # same-node iperf will cause error, as all the same-node communication are handled with different network interfaces. the original code lead to 30+Gbps bandwidith for all cases.
    os.system(f"ssh {server_host} 'iperf3 -s -D'")
    result = os.popen(f"ssh {test_server} 'iperf3 -c {ip_address} -t {time} -J'").read()
    data = json.loads(result)
    bandwidth = data["end"]["sum_received"]["bits_per_second"] / (2**30)
    os.system(f"ssh {server_host} 'pkill iperf3'")
    return bandwidth
    

def run_command(command):
    os.system(command)

def collect_network_usage(data_size, role, args, record_folder, keyword):
    monitor = SystemMonitor(0.01)
    monitor.start_all(interface=NETWORK_INTERFACE[0])

    threads = []
    role_assignment = [role, (role + 1) % 3, (role + 2) % 3]
    index = list(range(3))
    for i in range(3):
        index[role_assignment[i]] = i
    for i in range(3):
        command = f"{root_folder}out/build/linux/frontend/frontend -dataSize {data_size} -role {role_assignment[i]} {args} -p0_ip {IP_ADDRESS[index[0]]} -p1_ip {IP_ADDRESS[index[1]]} -rank 0"
        if i != 0:
            command = f"ssh {SERVER_HOST[i]} " + command
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
    parser.add_argument('--data_size', type=int, help='data size')
    parser.add_argument('--get_bandwidth_time', type=int, default=5, help='get bandwidth time')
    parser.add_argument('--skip_monitor', action='store_true', help='skip monitor')
    parser.add_argument('--fitting_length', type=int, default=10, help='fitting length')
    parser.add_argument('--fitting_step', type=int, default=100, help='fitting step')
    parser.add_argument('--fitting_degree', type=int, default=1, help='fitting degree')
    parser.add_argument('--parallelism_limit', type=int, default=64, help='parallelism limit')
    parser.add_argument('--run_tasks', action='store_true', help='run tasks')

    args = parser.parse_args()
    
    # mkdir the record folder.
    if not os.path.exists(args.record_folder):
        os.makedirs(args.record_folder)

    # get bandwidth
    print("Getting bandwidth")
    bandwidth = []
    for i in range(3):
        bandwidth.append(get_bandwidth(i, args.get_bandwidth_time))
    for i in range(3):
        print(f"{SERVER_HOST[i]}: {bandwidth[i]}Gb/s")
    
    
    # collect network usage
    print("Collecting network usage")
    length = args.fitting_length
    step = args.fitting_step
    if not args.skip_monitor:
        for role in range(3):
            sizes = range(step, step * length + 1, step)
            for size in sizes:
                collect_network_usage(size, role, args.args, args.record_folder, f"{args.keyword}-{role}-{size}")
                time.sleep(0.5)
    
    # fit network usage with polynomial
    print("Fitting network usage")
    degree = args.fitting_degree
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
            # print(f"Data size: {size} Mean recv: {mean_recv[-1]} Mean send: {mean_send[-1]}")
        coef_recv.append(np.polyfit(np.log(sizes), mean_recv, degree))
        coef_send.append(np.polyfit(np.log(sizes), mean_send, degree))
        # print(f"Data size: {np.array(sizes)}")
        # print(f"Mean recv: {np.array(mean_recv)}")
        # print(f"Mean send: {np.array(mean_send)}")
        # print(f"MSE recv: {np.mean((np.polyval(coef_recv[-1], np.log(sizes)) - mean_recv)**2)}")
        # print(f"MSE send: {np.mean((np.polyval(coef_send[-1], np.log(sizes)) - mean_send)**2)}")
    # print("Coef recv:")
    # print(coef_recv)
    # print("Coef send:")
    # print(coef_send)

    # decide the degree of parallelism
    sum_coef_recv = np.zeros(degree + 1)
    sum_coef_send = np.zeros(degree + 1)
    for i in range(3):
        sum_coef_recv += coef_recv[i]
        sum_coef_send += coef_send[i]
    
    data_size = args.data_size
    left = 1
    right = min(data_size // (3 * step), args.parallelism_limit)
    while left < right:
        mid = (left + right) // 2
        sum_recv = np.polyval(sum_coef_recv, np.log(data_size / (mid * 3)))
        sum_send = np.polyval(sum_coef_send, np.log(data_size / (mid * 3)))
        bandwidth_usage = max(sum_recv, sum_send) * mid
        if bandwidth_usage <= bandwidth[0]:
            left = mid + 1
        else:
            right = mid
    task_num = left
    subtask_data_size = []
    for i in range(task_num * 3):
        subtask_data_size.append((data_size + i) // (task_num * 3))
    print(f"Degree of parallelism: {task_num}")

    # run the tasks
    if args.run_tasks:
        print("Running tasks")
        monitor = SystemMonitor(0.01)
        monitor.start_all(interface=NETWORK_INTERFACE[0])

        threads = []
        commands = [[] for i in range(3)]
        for node in range(3):
            for rank in range(task_num):
                for role in range(3):
                    index = list(range(3))
                    for i in range(3):
                        index[(role + i) % 3] = i
                    p0_ip = IP_ADDRESS[index[0]]
                    p1_ip = IP_ADDRESS[index[1]]
                    subtask_size = subtask_data_size[rank * 3 + role]
                    command = f"{root_folder}out/build/linux/frontend/frontend -dataSize {subtask_size} -role {(role + node) % 3} {args.args} -rank {rank * 3 + role} -p0_ip {p0_ip} -p1_ip {p1_ip}"
                    commands[node].append(command)
        
        threads = []
        for node in range(3):
            command = " & ".join(commands[node]) + " & wait"
            if node != 0:
                command = f"ssh {SERVER_HOST[node]} '{command}'"
            print(command)
            thread = threading.Thread(target=run_command, args=(command,))
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join()
        
        monitor.stop_and_output(f"{args.record_folder}/monitor-{args.keyword}-{data_size}.log")
        usage_dict = get_usage_dict(f"{args.record_folder}/monitor-{args.keyword}-{data_size}.log")
        draw_usage_graph(usage_dict, f"{args.record_folder}/{args.keyword}-{data_size}.png")