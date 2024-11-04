root_folder=/root/aby3
result_folder=${root_folder}/scheduling/Result

# rm -r ${root_folder}/scheduling/Record_test/record.xlsx
port=7897
num_parties=3
server_host="aby30 aby31 aby32"
fitting_length=8
fitting_step=256
complexity="1 n"
get_bandwidth_time=1
parallelism_limit=96

if [ ! -d ${result_folder} ]; then
    mkdir ${result_folder}
fi

# # Sync the schedule
# scp -r ./scheduling aby31:${root_folder}/ &
# scp -r ./scheduling aby32:${root_folder}/ &
# wait;

# # prepare the test cpp.
# cp ${root_folder}/frontend/main.test ${root_folder}/frontend/main.cpp
# python ${root_folder}/build.py
# scp -r ${root_folder}/out/build/linux/frontend/frontend aby31:${root_folder}/out/build/linux/frontend/ &
# scp -r ${root_folder}/out/build/linux/frontend/frontend aby32:${root_folder}/out/build/linux/frontend/ &
# wait;

declare -A ip_addr_dict
ip_addr_dict["Homo-35G"]="10.3.0.12 10.3.0.14 10.3.0.4"
ip_addr_dict["Homo-10G"]="10.5.0.12 10.5.0.14 10.5.0.4"
ip_addr_dict["Homo-1G"]="10.1.0.12 10.1.0.14 10.1.0.4"
ip_addr_dict["Hetero-10G-35G-35G"]="10.5.0.12 10.3.0.14 10.3.0.4"
ip_addr_dict["Hetero-10G-10G-35G"]="10.5.0.12 10.5.0.14 10.3.0.4"
ip_addr_dict["Hetero-1G-35G-35G"]="10.1.0.12 10.3.0.14 10.3.0.4"
ip_addr_dict["Hetero-1G-1G-35G"]="10.1.0.12 10.1.0.14 10.3.0.4"
ip_addr_dict["Hetero-1G-10G-10G"]="10.1.0.12 10.5.0.14 10.5.0.4"
ip_addr_dict["Hetero-1G-1G-10G"]="10.1.0.12 10.1.0.14 10.5.0.4"

net_config_list=("Homo-35G" "Homo-10G" "Homo-1G" 
"Hetero-10G-35G-35G" "Hetero-10G-10G-35G" "Hetero-1G-35G-35G" "Hetero-1G-10G-10G" "Hetero-1G-1G-10G")

net_config_list=("Homo-35G")

network_config_list=("Homo-10G")
task_list=("Matrix" "Sort")
# data_size=33554432
data_size=4194304

for task in ${task_list[@]}; do
    # ip_address=${ip_addr_dict["Homo-35G"]}
    # network_interface=""
    # for ip in $ip_address; do
    #     interface=$(ssh $ip -o StrictHostKeyChecking=no -p $port "ip -o -4 addr show | grep $ip" | awk '{print $2}')
    #     network_interface="$network_interface $interface"
    # done

    python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${task} --task ${task} \
    --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
    --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} 

    for net_config in ${net_config_list[@]}; do
        ip_address=${ip_addr_dict[${net_config}]}
        keyword="${task}-${net_config}"
        bash ${root_folder}/scheduling/test_profiler.sh ${ip_address} ${keyword} ${task} ${data_size} ${net_config}
    done
done
