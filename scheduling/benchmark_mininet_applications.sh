NETNAME=$1

root_folder=/root/aby3

num_parties=3
# node_id=(12 14 4)
node_id=(11 12 13)
server_host="aby30 aby31 aby32"
parallelism_limit=96

# prepare the test cpp.
cp ${root_folder}/frontend/main.test ${root_folder}/frontend/main.cpp
python ${root_folder}/build.py
wait;

declare -A ip_addr_dict
ip_addr_dict[${NETNAME}]="10.0.0.${node_id[0]} 10.0.0.${node_id[1]} 10.0.0.${node_id[2]}"


net_config_list=(${NETNAME})
task_list=("LogReg-0" "LogReg-1" "LogReg-2" "Sort" "ORAM")
# task_list=("Shuffle")
assignment_strategy_list=("roundrole" "baseline")

declare -A data_size
data_size["Matrix"]=1048576
data_size["Sort"]=4194304
data_size["ORAM"]=33554432
data_size["LogReg-0"]=67108864
data_size["LogReg-1"]=67108864
data_size["LogReg-2"]=67108864
data_size["Shuffle"]=4194304

if [ $NETNAME == "Hetero-"* ]; then
    echo "Network bandwidth is 10G, scaling up data sizes..."
    data_size["Matrix"]=$((data_size["Matrix"] * 2))
    data_size["Sort"]=$((data_size["Sort"] * 4))
    data_size["Shuffle"]=$((data_size["Shuffle"] * 4))
    data_size["ORAM"]=$((data_size["ORAM"] * 2))
    data_size["LogReg-0"]=$((data_size["LogReg-0"] * 4))
    data_size["LogReg-1"]=$((data_size["LogReg-1"] * 4))
    data_size["LogReg-2"]=$((data_size["LogReg-2"] * 4))
fi


for task in ${task_list[@]}; do
    for net_config in ${net_config_list[@]}; do
        for assignment_strategy in ${assignment_strategy_list[@]}; do
            ip_address=${ip_addr_dict[${net_config}]}
            keyword="${task}-${net_config}-${assignment_strategy}"
            bash ${root_folder}/scheduling/test_profiler.sh ${ip_address} ${ip_address} ${keyword} ${task} ${data_size[${task}]} ${assignment_strategy} ${net_config} ${parallelism_limit} true "True"
        done
    done
done

cp -r ${root_folder}/scheduling/Record_test ${root_folder}/scheduling/Result/
mv ${root_folder}/scheduling/Result/Record_test ${root_folder}/scheduling/Result/Record_${NETNAME}
rm -rf ${root_folder}/scheduling/Record_test