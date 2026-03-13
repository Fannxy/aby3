#!/bin/bash

# 清理可能存在的旧 frontend 进程，避免端口冲突
pkill -f "frontend.*-role" || true
sleep 1  # 等待进程完全退出

# 测量 build 时间
BUILD_START=$(date +%s%N)
cp ./frontend/main.test ./frontend/main.cpp
python build.py 
BUILD_END=$(date +%s%N)
BUILD_DURATION=$(awk "BEGIN {printf \"%.3f\", ($BUILD_END - $BUILD_START) / 1000000000}")  # 转换为秒，保留3位小数


# ./out/build/linux/frontend/frontend  -role 0 &
# ./out/build/linux/frontend/frontend  -role 1 &
# ./out/build/linux/frontend/frontend  -role 2 &
# wait
# 设置监测参数
MONITOR_INTERVAL=0.1  # 采样间隔（秒）
# Mininet 中每个 host 的接口名称不同：h1-eth0, h2-eth0, h3-eth0
NETWORK_INTERFACE_H1=h1-eth0  # h1 的网络接口名称
NETWORK_INTERFACE_H2=h2-eth0  # h2 的网络接口名称
NETWORK_INTERFACE_H3=h3-eth0  # h3 的网络接口名称
MONITOR_DIR=./Monitor  # 监测数据保存目录

# 创建监测目录
mkdir -p ${MONITOR_DIR}

# 在每个 host 上启动监测脚本（后台运行）
echo "启动监测脚本..."
ssh h1 "cd /root/plantest/include/external-lib/aby3/; mkdir -p ${MONITOR_DIR}; nohup python3 Eval/monitor_frontend.py ${MONITOR_DIR}/monitor_h1_role0.json ${NETWORK_INTERFACE_H1} ${MONITOR_INTERVAL} > ${MONITOR_DIR}/monitor_h1.log 2>&1 & echo \$! > ${MONITOR_DIR}/monitor_h1.pid"
ssh h2 "cd /root/plantest/include/external-lib/aby3/; mkdir -p ${MONITOR_DIR}; nohup python3 Eval/monitor_frontend.py ${MONITOR_DIR}/monitor_h2_role1.json ${NETWORK_INTERFACE_H2} ${MONITOR_INTERVAL} > ${MONITOR_DIR}/monitor_h2.log 2>&1 & echo \$! > ${MONITOR_DIR}/monitor_h2.pid"
ssh h3 "cd /root/plantest/include/external-lib/aby3/; mkdir -p ${MONITOR_DIR}; nohup python3 Eval/monitor_frontend.py ${MONITOR_DIR}/monitor_h3_role2.json ${NETWORK_INTERFACE_H3} ${MONITOR_INTERVAL} > ${MONITOR_DIR}/monitor_h3.log 2>&1 & echo \$! > ${MONITOR_DIR}/monitor_h3.pid"

# 等待一下确保监测脚本启动
sleep 2

# 测量运行时间
RUN_START=$(date +%s%N)
# 运行 frontend 程序
echo "启动 frontend 程序..."
#./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} &
ssh h1 "cd /root/plantest/include/external-lib/aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list}" &
frontend_h1_pid=$!

ssh h2 "cd /root/plantest/include/external-lib/aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list}" &
frontend_h2_pid=$!

ssh h3 "cd /root/plantest/include/external-lib/aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list}" &
frontend_h3_pid=$!

# 等待所有 frontend 程序完成
wait;

# 停止监测脚本
echo "停止监测脚本..."
ssh h1 "cd /root/plantest/include/external-lib/aby3/; if [ -f ${MONITOR_DIR}/monitor_h1.pid ]; then kill -TERM \$(cat ${MONITOR_DIR}/monitor_h1.pid) 2>/dev/null; rm -f ${MONITOR_DIR}/monitor_h1.pid; fi || pkill -f 'monitor_frontend.py.*monitor_h1_role0'" 2>/dev/null
ssh h2 "cd /root/plantest/include/external-lib/aby3/; if [ -f ${MONITOR_DIR}/monitor_h2.pid ]; then kill -TERM \$(cat ${MONITOR_DIR}/monitor_h2.pid) 2>/dev/null; rm -f ${MONITOR_DIR}/monitor_h2.pid; fi || pkill -f 'monitor_frontend.py.*monitor_h2_role1'" 2>/dev/null
ssh h3 "cd /root/plantest/include/external-lib/aby3/; if [ -f ${MONITOR_DIR}/monitor_h3.pid ]; then kill -TERM \$(cat ${MONITOR_DIR}/monitor_h3.pid) 2>/dev/null; rm -f ${MONITOR_DIR}/monitor_h3.pid; fi || pkill -f 'monitor_frontend.py.*monitor_h3_role2'" 2>/dev/null

# 等待监测脚本完成数据保存
sleep 3

echo "监测完成，数据保存在 ${MONITOR_DIR}/ 目录下"


RUN_END=$(date +%s%N)
RUN_DURATION=$(awk "BEGIN {printf \"%.3f\", ($RUN_END - $RUN_START) / 1000000000}")  # 转换为秒，保留3位小数

# 清理 frontend 进程
pkill -f "frontend.*-role" || true

# 写入时间结果到文件（如果环境变量设置了）
if [ ! -z "$FEDDB_ABY3_TIME_FILE" ]; then
    cat > "$FEDDB_ABY3_TIME_FILE" << EOF
Aby3 Dispatch Time Measurement
===============================
Build Time: $BUILD_DURATION s
Run Time: $RUN_DURATION s
EOF
fi