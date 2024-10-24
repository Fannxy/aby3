import sys
import os
import argparse
import threading
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PtA_deploy')))
from system_monitor import *

root_folder = "/root/aby3/"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--keyword', type=str, help='keyword')
    parser.add_argument('--command', type=str, nargs='+', help='run command')
    parser.add_argument('--record_folder', type=str, help='record folder')
    parser.add_argument('--interface', type=str, help='network interface')
    
    args = parser.parse_args()

    if(not os.path.exists(args.record_folder)):
        os.makedirs(args.record_folder)
    
    monitor = SystemMonitor(0.01)
    monitor.start_all(interface=args.interface)
    command = " ".join(args.command).replace('+', '-')
    os.system(command)
    monitor.stop_and_output(args.record_folder + f"/monitor-{args.keyword}.log")
    exit(0)
