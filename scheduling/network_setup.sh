# start network and set the bandwidth
bw_list="$1 $2 $3"
mn -c
root_folder=/root/aby3
log_file=${root_folder}/network_setup.log

MININET_SESSION="mininet"

# 检查 tmux 会话是否存在
if tmux has-session -t $MININET_SESSION 2>/dev/null; then
    echo "tmux session $MININET_SESSION already exists"
else
    # 创建一个新的 tmux 会话
    tmux new-session -d -s $MININET_SESSION
    echo "tmux session $MININET_SESSION created"
fi

cd ${root_folder}/Net/mininet;
# tmux new-session -d -s $MININET_SESSION "python ./examples/p2p_3pc_net.py --ip 10.1.0.12 --bw $bw_list; bash"
tmux send-keys -t $MININET_SESSION "python ./examples/p2p_3pc_net.py --ip 10.1.0.12 --bw $bw_list; bash" C-m

# wait for the network to be set up
# 等待网络设置完成
hosts=("aby30" "aby31" "aby32")
for host in "${hosts[@]}"; do
    while ! ssh -o ConnectTimeout=2 -o StrictHostKeyChecking=no $host "exit" 2>/dev/null; do
        echo "Waiting for $host to be ready..." >> $log_file
        sleep 5
    done
    echo "$host is ready" >> $log_file
done

echo "All hosts are ready" >> $log_file