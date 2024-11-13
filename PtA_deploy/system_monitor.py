import psutil
import threading
import time
import json
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import os
import numpy as np

def draw_usage_graph(usage_dict, graph_file_name):
    """
    Draw the usage graph.
    """
    plt.figure(figsize=(16, 4))
                
    plt.subplot(1, 4, 1)
    plt.plot(usage_dict["cpu"], label="cpu")
    plt.plot(usage_dict["user_cpu"], label="user")
    plt.plot(usage_dict["sys_cpu"], label="sys")
    print(f"stamps of cpu = {len(usage_dict['cpu'])}")
    plt.title("cpu usage")
    plt.xlabel("time")
    plt.ylabel("usage (%)")
    plt.legend()
    
    plt.subplot(1, 4, 2)
    plt.plot(usage_dict["memory"])
    print(f"stamps of memory = {len(usage_dict['memory'])}")
    plt.title("memory usage")
    plt.xlabel("time")
    plt.ylabel("usage (GB)")
    
    plt.subplot(1, 4, 3)
    plt.plot(usage_dict["network_send"], label="network send")
    plt.plot(usage_dict["network_recv"], label="network recv")
    print(f"stamps of network = {len(usage_dict['network_send'])}")
    plt.title("network usage")
    plt.xlabel("time")
    plt.ylabel("usage (Gb/s)")
    plt.legend()
    
    plt.subplot(1, 4, 4)
    plt.plot(usage_dict["recv_drop"], label="drop")
    plt.title("network recvdrop")
    plt.xlabel("time")
    plt.ylabel("drop rate")
    plt.legend()
    
    plt.tight_layout()
    print(graph_file_name)
    plt.savefig(graph_file_name, dpi=300)
    
    return

def get_usage_dict(monitor_log):
    usage_dict = {
        "cpu": [],
        "user_cpu": [],
        "sys_cpu": [],
        "memory": [],
        "network_send": [],
        "network_recv": [],
        "recv_drop": [],
    }
    
    with open(monitor_log, 'r') as file:
        tmp_dict = json.load(file)
    
    for key in tmp_dict.keys():
        if(key != "network" and key != "user_cpu" and key != "sys_cpu"):
            usage_dict[key] = [item[1] for item in tmp_dict[key]]
        else:
            send_list = [item[1][0] for item in tmp_dict["network"]]
            recv_list = [item[1][1] for item in tmp_dict["network"]]
            recv_drop_list = [item[1][6] for item in tmp_dict["network"]]
            time_stamp_list = [item[0] for item in tmp_dict["network"]]
            
            diff_send = np.diff(send_list)
            diff_recv = np.diff(recv_list)
            diff_drop = np.diff(recv_drop_list)
            diff_time = np.diff(time_stamp_list)
            
            send_speed = (diff_send / (diff_time * (2**30))) * 8
            recv_speed = (diff_recv / (diff_time * (2**30))) * 8
            drop_rate = diff_drop / diff_time
            
            usage_dict["network_recv"] = recv_speed.tolist()
            usage_dict["network_send"] = send_speed.tolist()
            usage_dict["recv_drop"] = drop_rate.tolist()
            
            sys_cpu_list = [item[1] for item in tmp_dict["sys_cpu"]]
            user_cpu_list = [item[1] for item in tmp_dict["user_cpu"]]
            time_stamp_list_cpu = [item[0] for item in tmp_dict["sys_cpu"]]
            
            diff_time_cpu = np.diff(time_stamp_list_cpu)
            sys_cpu_list = np.diff(sys_cpu_list)
            user_cpu_list = np.diff(user_cpu_list)
            sys_cpu_rate = np.array(sys_cpu_list) / diff_time_cpu
            user_cpu_rate = np.array(user_cpu_list) / diff_time_cpu
            
            usage_dict["sys_cpu"] = sys_cpu_rate.tolist()
            usage_dict["user_cpu"] = user_cpu_rate.tolist()
    
    return usage_dict


class SystemMonitor:
    
    def __init__(self, interval):
        self.interval = interval
        self.cpu_list = []
        self.network_list = []
        self.memory_list = []
        self.sys_cpu_list = []
        self.user_cpu_list = []
        self.cpu_running = False
        self.network_running = False
        self.memory_running = False
    
    def start_cpu(self):
        self.cpu_running = True
        self.thread_cpu = threading.Thread(target=self.record_cpu)
        self.thread_cpu.start()
    
    def record_cpu(self):
        while self.cpu_running:
            cpu_percent = psutil.cpu_percent()
            time_stamp = time.time()
            self.cpu_list.append((time_stamp, cpu_percent))
            self.sys_cpu_list.append((time_stamp, psutil.cpu_times().system))
            self.user_cpu_list.append((time_stamp, psutil.cpu_times().user))
            
            time.sleep(self.interval)
            
    def stop_cpu(self):
        self.cpu_running = False
        self.thread_cpu.join()
    
    def start_network(self, interface):
        self.network_running = True
        self.thread_network = threading.Thread(target=self.record_network, args=(interface,))
        self.thread_network.start()
    
    def record_network(self, interface):
        while self.network_running:
            network_percent = psutil.net_io_counters(pernic=True).get(interface)
            time_stamp = time.time()
            self.network_list.append((time_stamp, network_percent))
            time.sleep(self.interval)
    
    def stop_network(self):
        self.network_running = False
        self.thread_network.join()
    
    def start_memory(self):
        self.memory_running = True
        self.thread_memory = threading.Thread(target=self.record_memory)    
        self.thread_memory.start()
    
    def record_memory(self):
        while self.memory_running:
            memory_percent = (psutil.virtual_memory().used) / (2**30)
            time_stamp = time.time()
            self.memory_list.append((time_stamp, memory_percent))
            time.sleep(self.interval)
    
    def stop_memory(self):
        self.memory_running = False
        self.thread_memory.join()
        
    def start_all(self, interface):
        self.start_cpu()
        self.start_memory()
        self.start_network(interface)
    
    def output(self, save_file):
        with open(save_file, 'w') as file:
            json.dump({"cpu": self.cpu_list, "network": self.network_list, "memory": self.memory_list, "user_cpu": self.user_cpu_list, "sys_cpu": self.sys_cpu_list}, file)
        # clear the memory.
        self.cpu_list = []
        self.network_list = []
        self.memory_list = []
        self.sys_cpu_list = []
        self.user_cpu_list = []
    
    def stop_and_output(self, save_file):
        self.stop_cpu()
        self.stop_memory()
        self.stop_network()
        self.output(save_file)