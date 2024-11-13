ip_address="$1 $2 $3"
profile_ip_address="$4 $5 $6"
keyword=$7
task=$8
data_size=$9
assignment_strategy=${10}
net_config=${11}
micro_benchmark="${12:-""}"

root_folder=/root/aby3
port=7897
# port=1022
num_parties=3
server_host="aby30 aby31 aby32"
get_bandwidth_time=2
parallelism_limit=48

declare -A fitting_length
declare -A fitting_step
declare -A complexity
fitting_length["Matrix"]=8
fitting_step["Matrix"]=256
complexity["Matrix"]="1 n"
fitting_length["Sort"]=8
fitting_step["Sort"]=4096
complexity["Sort"]="1 n n*log(n)"

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
  --fitting_length ${fitting_length[${task}]} \
  --fitting_step ${fitting_step[${task}]} \
  --get_bandwidth_time ${get_bandwidth_time} \
  --parallelism_limit ${parallelism_limit} \
  --complexity ${complexity[${task}]} \
  --assignment_strategy ${assignment_strategy} \
  --run_tasks