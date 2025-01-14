NETNAME=$1
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

declare -A ip_addr_dict
ip_addr_dict[${NETNAME}]="10.0.0.${node_id[0]} 10.0.0.${node_id[1]} 10.0.0.${node_id[2]}"


# net_config_list=(${NETNAME})
# task_list=("Matrix" "Sort" "ORAM")
# assignment_strategy_list=("roundrole" "baseline")
assignment_strategy_list=("roundrole" "baseline")
# fix_balance_list=("False" "True")

net_config_list=(${NETNAME})
data_size_micro=1073741824
if [ $NETNAME = "Homo-10G" ]; then
    data_size_micro=2147483648
fi
if [[ $NETNAME == "Hetero-"* ]]; then
    data_size_micro=2147483648
fi
# micro_benchmarks=("ff-mul" "ib-mul" "a2b" "b2a-single" "shuffle")
# micro_benchmarks=("fake_test" "fake_test2")
micro_benchmarks=("fake_test" "fake_test2")
for task in ${micro_benchmarks[@]}; do
    for net_config in ${net_config_list[@]}; do
        ip_address=${ip_addr_dict[${net_config}]}
        for assignment_strategy in  ${assignment_strategy_list[@]}; do
            # for fix_strategy in ${fix_balance_list[@]}; do
            keyword="${task}-${net_config}-${assignment_strategy}"
            bash ${root_folder}/scheduling/test_profiler.sh ${ip_address} ${ip_address} ${keyword} ${task} ${data_size_micro} ${assignment_strategy} ${net_config} ${parallelism_limit} True -Micro
            # done
        done
    done
done
cp -r ${root_folder}/scheduling/Record_test ${root_folder}/scheduling/Result/
mv ${root_folder}/scheduling/Result/Record_test ${root_folder}/scheduling/Result/Record_${NETNAME}
rm -rf ${root_folder}/scheduling/Record_test