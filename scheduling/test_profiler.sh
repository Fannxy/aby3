ip_address="$1 $2 $3"
profile_ip_address="$4 $5 $6"
keyword=$7
task=$8
data_size=$9
assignment_strategy=${10}
net_config=${11}
parallelism_limit=${12}
fix_strategy=${13}
micro_benchmark="${14:-""}"

root_folder=/root/aby3
port=22
# port=1022
num_parties=3
server_host="aby30 aby31 aby32"
get_bandwidth_time=2

fitting_length=16
fitting_step=2048576
complexity="n"

network_interface=""
for ip in $ip_address; do
    interface=$(ssh $ip -o StrictHostKeyChecking=no -p $port "ip -o -4 addr show | grep $ip" | awk '{print $2}')
    network_interface="$network_interface $interface"
done

profile_network_interface=""
for ip in $profile_ip_address; do
    interface=$(ssh $ip -o StrictHostKeyChecking=no -p $port "ip -o -4 addr show | grep $ip" | awk '{print $2}')
    profile_network_interface="$profile_network_interface $interface"
done

echo "ip_address: $ip_address"
echo "keyword: $keyword"
echo "task: $task"
echo "data_size: $data_size"
echo "network_interface: $network_interface"
echo "assignment_strategy: $assignment_strategy"
echo "net_config: $net_config"

python ${root_folder}/scheduling/profiler.py \
  --args " ${micro_benchmark} -${task}" \
  --args_agg " -${task}-agg" \
  --record_folder ${root_folder}/scheduling/Record_test \
  --config_folder ${root_folder}/scheduling/Result/${net_config} \
  --keyword ${keyword} \
  --task ${task} \
  --num_parties ${num_parties} \
  --server_host ${server_host} \
  --ip_address ${ip_address} \
  --network_interface ${network_interface} \
  --profile_ip_address ${profile_ip_address} \
  --profile_network_interface ${profile_network_interface} \
  --data_size ${data_size} \
  --fitting_length ${fitting_length} \
  --fitting_step ${fitting_step} \
  --get_bandwidth_time ${get_bandwidth_time} \
  --parallelism_limit ${parallelism_limit} \
  --complexity ${complexity} \
  --assignment_strategy ${assignment_strategy} \
  --balance_fix ${fix_strategy} \
  --run_tasks