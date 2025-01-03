root_folder=/root/aby3

num_parties=3
# node_id=(12 14 4)
node_id=(11 12 13)
server_host="aby30 aby31 aby32"
parallelism_limit=64

# prepare the test cpp.
cp ${root_folder}/frontend/main.test ${root_folder}/frontend/main.cpp
python ${root_folder}/build.py
wait;

NETNAME="Homo-100M"

declare -A ip_addr_dict
ip_addr_dict[${NETNAME}]="10.0.0.${node_id[0]} 10.0.0.${node_id[1]} 10.0.0.${node_id[2]}"


# net_config_list=(${NETNAME})
# task_list=("Matrix" "Sort" "ORAM")
assignment_strategy_list=("roundrole" "baseline")
# declare -A data_size
# data_size["Matrix"]=1048576
# data_size["Sort"]=4194304
# data_size["ORAM"]=33554432

# for task in ${task_list[@]}; do
#     for net_config in ${net_config_list[@]}; do
#         for assignment_strategy in ${assignment_strategy_list[@]}; do
#             ip_address=${ip_addr_dict[${net_config}]}
#             keyword="${task}-${net_config}-${assignment_strategy}"
#             bash ${root_folder}/scheduling/test_profiler.sh ${ip_address} ${ip_address} ${keyword} ${task} ${data_size[${task}]} ${assignment_strategy} ${net_config} ${parallelism_limit}
#         done
#     done
# done

mkdir ${root_folder}/scheduling/Record_test

net_config_list=(${NETNAME})
data_size_micro=1073741824
micro_benchmarks=("shuffle" "ff-mul" "ib-mul" "a2b" "b2a-single")
for task in ${micro_benchmarks[@]}; do
    for net_config in ${net_config_list[@]}; do
        ip_address=${ip_addr_dict[${net_config}]}
        for assignment_strategy in  ${assignment_strategy_list[@]}; do
            keyword="${task}-${net_config}-${assignment_strategy}"
            bash ${root_folder}/scheduling/test_profiler.sh ${ip_address} ${ip_address} ${keyword} ${task} ${data_size_micro} ${assignment_strategy} ${net_config} ${parallelism_limit} -Micro
        done
    done
done

cp -r ${root_folder}/scheduling/Record_test ${root_folder}/scheduling/Result/
mv ${root_folder}/scheduling/Result/Record_test ${root_folder}/scheduling/Result/Record_${NETNAME}
# cp ${root_folder}/scheduling/Record_test/record.xlsx ${root_folder}/scheduling/Result/
# rm -rf ${root_folder}/scheduling/Record_test