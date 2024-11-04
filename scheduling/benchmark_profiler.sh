root_folder=/root/aby3
result_folder=${root_folder}/scheduling/Result

# rm -r ${root_folder}/scheduling/Record_test/record.xlsx

if [ ! -d ${result_folder} ]; then
    mkdir ${result_folder}
fi

# Sync the schedule
scp -r ./scheduling aby31:${root_folder}/ &
scp -r ./scheduling aby32:${root_folder}/ &
wait;

# prepare the test cpp.
cp ${root_folder}/frontend/main.test ${root_folder}/frontend/main.cpp
python ${root_folder}/build.py
scp -r ${root_folder}/out/build/linux/frontend/frontend aby31:${root_folder}/out/build/linux/frontend/ &
scp -r ${root_folder}/out/build/linux/frontend/frontend aby32:${root_folder}/out/build/linux/frontend/ &
wait;

num_parties=3
server_host="aby30 aby31 aby32"
ip_address="10.5.0.13 10.3.0.16 10.5.0.17"
network_interface="ens11 ibs110 ens110"

declare -A network_ip_dict
network_ip_dict["Homo-1G"]="10.1.0.12 10.1.0.14 10.1.0.4"
network_ip_dict["Homo-10G"]="10.5.0.12 10.5.0.14 10.5.0.4"
network_ip_dict["Homo-35G"]="10.3.0.12 10.3.0.14 10.3.0.4"
network_ip_dict["Heter-1G-10G"]="10.1.0.12 10.5.0.14 10.1.0.4"
network_ip_dict["Heter-1G-35G"]="10.1.0.12 10.3.0.14 10.1.0.4"
network_ip_dict["Heter-10G-35G"]="10.5.0.12 10.3.0.14 10.5.0.4"

declare -A network_interface_dict
network_interface_dict["Homo-1G"]="ens121f0 ens121f0 ens121f0"
network_interface_dict["Homo-10G"]="ens110 ens11 ens11"
network_interface_dict["Homo-35G"]="ibs110 ibs110 ibs110"
network_interface_dict["Heter-1G-10G"]="ens121f0 ens11 ens121f0"
network_interface_dict["Heter-1G-35G"]="ens121f0 ibs110 ens121f0"
network_interface_dict["Heter-10G-35G"]="ens110 ibs110 ens11"


fitting_length=16
fitting_step=128
complexity="1 n"
get_bandwidth_time=10
parallelism_limit=48

# task_list=("Index" "Max" "Metric")
# data_size=268435456

# for task in ${task_list[@]}; do
#     python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${task} \
#     --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
#     --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
#     --run_tasks --MPI

#     python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${task} \
#     --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
#     --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
#     --run_tasks --MPI --skip_monitor --baseline
# done

network_config_list=("Homo-1G" "Homo-10G" "Homo-35G" "Heter-1G-10G" "Heter-1G-35G" "Heter-10G-35G")
# network_config_list=("Homo-35G")
task_list=("Matrix" "Sort")
# data_size=33554432
data_size=4194304

for task in ${task_list[@]}; do
    # for i in 0 1 2 3; do
    for network_config in ${network_config_list[@]}; do
        ip_address=${network_ip_dict[$network_config]}
        network_interface=${network_interface_dict[$network_config]}

        skip_monitor="--skip_monitor"
        if [ $i -gt 0 ]; then
            skip_monitor="--skip_monitor"
        fi

        python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test_${network_config} --keyword ${task}-${i} \
        --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
        --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
        --run_tasks ${skip_monitor} --config_folder ${root_folder}/scheduling/config_${network_config}

        python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test_${network_config} --keyword ${task}-${i} \
        --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
        --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
        --run_tasks --skip_monitor --baseline --config_folder ${root_folder}/scheduling/config_${network_config}

        cp ${root_folder}/scheduling/Record_test_${network_config}/record.xlsx ${result_folder}/${task}-${network_config}.xlsx
    done
done
