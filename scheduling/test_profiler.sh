ip_address="$1 $2 $3"
keyword=$4
task=$5
data_size=$6

root_folder=/root/aby3
port=7897
num_parties=3
server_host="aby30 aby31 aby32"
fitting_length=16
fitting_step=128
complexity="1 n"
get_bandwidth_time=5
parallelism_limit=128

network_interface=""
for ip in $ip_address; do
    interface=$(ssh $ip -o StrictHostKeyChecking=no -p $port "ip -o -4 addr show | grep $ip" | awk '{print $2}')
    network_interface="$network_interface $interface"
done

echo "ip_address: $ip_address"
echo "keyword: $keyword"
echo "task: $task"
echo "data_size: $data_size"
echo "network_interface: $network_interface"

python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${keyword} --task ${task} \
--num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
--data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
--run_tasks --skip_monitor

python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${keyword}-baseline --task ${task} \
--num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
--data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
--run_tasks --skip_monitor --baseline