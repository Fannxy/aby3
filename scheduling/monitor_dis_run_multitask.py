import sys
import os
import argparse
import threading
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PtA_deploy')))
from system_monitor import *

root_folder = "/root/aby3/"
server1="aby31"
server2="aby32"
NETWORK_INTERFACE0 = "ens121f0"
NETWORK_INTERFACE1 = "ens121f0"
NETWORK_INTERFACE2 = "ens121f0"
IP_ADDRESS0 = "10.1.0.15"
IP_ADDRESS1 = "10.1.0.16"
IP_ADDRESS2 = "10.1.0.17"
# NETWORK_INTERFACE0 = "ens110"
# NETWORK_INTERFACE1 = "ens11"
# NETWORK_INTERFACE2 = "ens11"
# IP_ADDRESS0 = "10.5.0.15"
# IP_ADDRESS1 = "10.5.0.16"
# IP_ADDRESS2 = "10.5.0.17"

def run_command(command):
    os.system(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--keyword', type=str, help='keyword')
    parser.add_argument('--args', type=str, help='run args')
    parser.add_argument('--record_folder', type=str, help='record folder')
    parser.add_argument('--role', type=int, help='role')
    parser.add_argument('--analysis', action='store_true', help='analysis')
    parser.add_argument('--symmetric', action='store_true', help='symmetric')
    parser.add_argument('--order', type=str, action='append', help='order')
    parser.add_argument('--repeat', type=int, help='repeat')
    
    args = parser.parse_args()
        
    # if(os.path.exi args.record_folder)
    if(not os.path.exists(args.record_folder)):
        os.makedirs(args.record_folder)
        
    NETWORK_INTERFACE = eval(f"NETWORK_INTERFACE{args.role}")
        
    
    if(not args.analysis):
        monitor = SystemMonitor(0.01)
        monitor.start_all(interface=NETWORK_INTERFACE)

        threads = []
        role_assignments = []
        for role_assignment in args.order:
            role_assignments.append(list(map(int, role_assignment.split(','))))
        role_assignments = role_assignments * args.repeat
        for rank, role_assignment in enumerate(role_assignments):
            role = role_assignment[args.role]
            index = list(range(3))
            for i in range(3):
                index[role_assignment[i]] = i
            p0_ip = eval(f"IP_ADDRESS{index[0]}")
            p1_ip = eval(f"IP_ADDRESS{index[1]}")
            command = f"{root_folder}out/build/linux/frontend/frontend -role {role} {args.args} -rank {rank} -p0_ip {p0_ip} -p1_ip {p1_ip}"
            print(command)
            thread = threading.Thread(target=run_command, args=(command,))
            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()

        monitor.stop_and_output(args.record_folder + f"/monitor-{args.keyword}.log")
        exit(0)

    if(args.analysis):  
        if(args.role == 0):
            os.system(f"mv {args.record_folder}/monitor-{args.keyword}.log {args.record_folder}/monitor-{args.keyword}-0.log")
            os.system(f"scp -r aby31:{args.record_folder}/monitor-{args.keyword}.log {args.record_folder}/monitor-{args.keyword}-1.log")
            os.system(f"scp -r aby32:{args.record_folder}/monitor-{args.keyword}.log {args.record_folder}/monitor-{args.keyword}-2.log")
            
            print("scp done!!!!!!")
            # analyze the monitors.
            for i in range(3):
                usage_dict = get_usage_dict(args.record_folder + f"/monitor-{args.keyword}-{i}.log")
                # print(usage_dict)
                draw_usage_graph(usage_dict, args.record_folder + f"/{args.keyword}-{i}.png")
    
    