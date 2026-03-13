#!/usr/bin/env python3
"""
监测脚本：使用 SystemMonitor 类监测系统资源使用情况
用法: python monitor_frontend.py <output_file> <interface> <interval>
"""
import sys
import os
import signal
import time
import json

# 添加当前脚本所在目录到 Python 路径，以便导入 system_monitor
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from system_monitor import SystemMonitor, get_usage_dict, draw_usage_graph

def signal_handler(sig, frame):
    """处理停止信号"""
    global monitor
    print(f"\n收到停止信号，正在保存数据...")
    monitor.stop_and_output(output_file)
    print(f"数据已保存到: {output_file}")
    
    # 生成图表
    try:
        # 检查文件是否存在
        if not os.path.exists(output_file):
            print(f"警告: 数据文件不存在: {output_file}", file=sys.stderr)
            return
        
        # 检查文件是否为空
        if os.path.getsize(output_file) == 0:
            print(f"警告: 数据文件为空: {output_file}", file=sys.stderr)
            return
        
        print("正在生成监测图表...")
        usage_dict = get_usage_dict(output_file)
        
        # 检查数据是否有效
        if not usage_dict or len(usage_dict.get("cpu", [])) == 0:
            print("警告: 监测数据为空，跳过图表生成", file=sys.stderr)
            return
        
        # 生成图表文件名（将 .json 替换为 .png）
        graph_file = output_file.rsplit('.', 1)[0] + '.png'
        draw_usage_graph(usage_dict, graph_file)
        print(f"图表已保存到: {graph_file}")
    except FileNotFoundError:
        print(f"错误: 找不到数据文件: {output_file}", file=sys.stderr)
    except json.JSONDecodeError as e:
        print(f"错误: JSON 数据格式错误: {e}", file=sys.stderr)
    except Exception as e:
        print(f"生成图表时出错: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    
    sys.exit(0)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("用法: python monitor_frontend.py <output_file> <interface> <interval>")
        print("示例: python monitor_frontend.py /tmp/monitor.json eth0 0.1")
        sys.exit(1)
    
    output_file = sys.argv[1]
    interface = sys.argv[2]
    interval = float(sys.argv[3])
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # 创建监测器
    monitor = SystemMonitor(interval=interval)
    
    # 注册信号处理器
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    # 启动所有监测
    print(f"开始监测 - 输出文件: {output_file}, 接口: {interface}, 间隔: {interval}秒")
    monitor.start_all(interface)
    
    try:
        # 持续运行直到收到停止信号
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        signal_handler(None, None)

        

