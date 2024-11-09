root_folder=/root/aby3

num_parties=3
# node_id=(12 14 4)
node_id=(13 16 17)
server_host="aby30 aby31 aby32"

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
ip_addr_dict["Homo-35G"]="10.3.0.${node_id[0]} 10.3.0.${node_id[1]} 10.3.0.${node_id[2]}"
ip_addr_dict["Homo-10G"]="10.5.0.${node_id[0]} 10.5.0.${node_id[1]} 10.5.0.${node_id[2]}"
ip_addr_dict["Homo-1G"]="10.1.0.${node_id[0]} 10.1.0.${node_id[1]} 10.1.0.${node_id[2]}"
ip_addr_dict["Hetero-10G-35G-35G"]="10.5.0.${node_id[0]} 10.3.0.${node_id[1]} 10.3.0.${node_id[2]}"
ip_addr_dict["Hetero-10G-10G-35G"]="10.5.0.${node_id[0]} 10.5.0.${node_id[1]} 10.3.0.${node_id[2]}"
ip_addr_dict["Hetero-1G-35G-35G"]="10.1.0.${node_id[0]} 10.3.0.${node_id[1]} 10.3.0.${node_id[2]}"
ip_addr_dict["Hetero-1G-1G-35G"]="10.1.0.${node_id[0]} 10.1.0.${node_id[1]} 10.3.0.${node_id[2]}"
ip_addr_dict["Hetero-1G-10G-10G"]="10.1.0.${node_id[0]} 10.5.0.${node_id[1]} 10.5.0.${node_id[2]}"
ip_addr_dict["Hetero-1G-1G-10G"]="10.1.0.${node_id[0]} 10.1.0.${node_id[1]} 10.5.0.${node_id[2]}"
ip_addr_dict["Hetero-35G-35G-10G"]="10.3.0.${node_id[0]} 10.3.0.${node_id[1]} 10.5.0.${node_id[2]}"

net_config_list=("Homo-35G" "Homo-10G" "Homo-1G" 
"Hetero-10G-35G-35G" "Hetero-10G-10G-35G" "Hetero-1G-35G-35G" "Hetero-1G-1G-35G" "Hetero-1G-10G-10G" "Hetero-1G-1G-10G")
task_list=("Matrix" "Sort")
assignment_strategy_list=("roundrole" "baseline")
data_size=33554432

net_config_list=("Homo-10G" "Hetero-10G-35G-35G" "Hetero-10G-10G-35G")
task_list=("Matrix" "Sort")
assignment_strategy_list=("roundrole" "baseline")
data_size=4194304

for task in ${task_list[@]}; do
    for net_config in ${net_config_list[@]}; do
        for assignment_strategy in ${assignment_strategy_list[@]}; do
            ip_address=${ip_addr_dict[${net_config}]}
            keyword="${task}-${net_config}-${assignment_strategy}"
            bash ${root_folder}/scheduling/test_profiler.sh ${ip_address} ${ip_addr_dict["Homo-35G"]} ${keyword} ${task} ${data_size} ${assignment_strategy} ${net_config} $((parallelism_limit * 3))
        done
    done
done
