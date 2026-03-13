args_list=$1
echo ${args_list}

# ./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list} &
# wait;


# THE FOLLOWING IS FOR DISTRIBUTED TEST
# scp ./out/build/linux/frontend/frontend h2:/root/plantest/include/external-lib/aby3/out/build/linux/frontend &
# scp ./out/build/linux/frontend/frontend h3:/root/plantest/include/external-lib/aby3/out/build/linux/frontend &
# wait;

#设置监测参数
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

# 运行 frontend 程序
echo "启动 frontend 程序..."
#./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} &
ssh h1 "cd /root/plantest/include/external-lib/aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} "&
frontend_h1_pid=$!

ssh h2 "cd /root/plantest/include/external-lib/aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list} "&
frontend_h2_pid=$!

ssh h3 "cd /root/plantest/include/external-lib/aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list} "&
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
