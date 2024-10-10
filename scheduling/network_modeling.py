import os
import argparse
import threading
import subprocess
import queue
import json

root_folder = "/root/aby3/"
SERVER_HOST = ["aby30", "aby31", "aby32"]
HOSTNAME = ["10.1.0.12", "10.1.0.14", "10.1.0.4"]
NETWORK_INTERFACE = ["ibs110", "ibs110", "ibs110"] # do not use
IP_ADDRESS = ["10.3.0.12", "10.3.0.14", "10.3.0.4"]
START_PORT = 5201
TIME_INTERVAL = 60
network_res_file = root_folder + "scheduling/network.json"


def server_start(hostname, target_links):
    iperf_start_args_prefix = f"iperf3 -s -D -p"
    for i in range(target_links):
        port = START_PORT + i
        os.system(f"ssh {hostname} '{iperf_start_args_prefix} {port}'")
    return


def run_iperf3_client(hostname, target_ip_address, time, port, result_queue):
    command = f"ssh {hostname} 'iperf3 -c {target_ip_address} -t {time} -p {port} -J'"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
    result_queue.put((hostname, target_ip_address, result))
    return


def get_bandwidth(hostname_list, target_ip_addresses_list, target_ports_list, time):
    # ssh into the target host (string), and run iperf3 client command for the target ip addresses (list of strings) for the given time (int).
    threads = []
    result_queue = queue.Queue()
        
    for i in range(len(hostname_list)):
        hostname = hostname_list[i]
        target_ip_addresses = target_ip_addresses_list[i]
        target_ports = target_ports_list[i]
        for j in range(len(target_ip_addresses)):
            ip_address = target_ip_addresses[j]
            port = target_ports[j]
            thread = threading.Thread(target=run_iperf3_client, args=(hostname, ip_address, time, port, result_queue))
            threads.append(thread)
            thread.start()
    
    for thread in threads:
        thread.join()

    # analyze the result.
    return result_queue


def kill_iperf3_process(hostname):
    os.system(f"ssh {hostname} 'pkill iperf3'")
    return

    
if __name__ == "__main__":
    
    # get parameters.
    N = len(SERVER_HOST)
    
    # start the iperf3 server on each host.
    for host in SERVER_HOST:
        server_start(host, N)
    
    # run the iperf3 client on each host.
    ip_address_list = []
    ports_list = []
    for i in range(N):
        target_ip_addresses = [IP_ADDRESS[(i+j)%N] for j in range(1, N)]
        target_ports = [START_PORT + i for _ in range(N-1)]
        ip_address_list.append(target_ip_addresses)
        ports_list.append(target_ports)
    
    result_queue = get_bandwidth(SERVER_HOST, ip_address_list, ports_list, TIME_INTERVAL)
    
    # kill the iperf3 server on each host.
    for host in SERVER_HOST:
        kill_iperf3_process(host)
    
    # analyze the result.
    bandwidth_dict = {}
    while not result_queue.empty():
        hostname, target_ip_address, result = result_queue.get()
        res_dict = json.loads(result.stdout)
        bandwidth = res_dict["end"]["sum_received"]["bits_per_second"] / (2**30)
        if hostname not in bandwidth_dict:
            bandwidth_dict[hostname] = {target_ip_address: bandwidth}
        else:
            bandwidth_dict[hostname][target_ip_address] = bandwidth
    
    with open(network_res_file, "w") as f:
        json.dump(bandwidth_dict, f)
    
    