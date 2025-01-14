NETNAME="Net-${1}-${2}-${3}M"

root_folder=/root/aby3

# test the roundrole assignment strategy
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

assignment_strategy_list=("roundrole" "baseline")

net_config_list=(${NETNAME})

data_size_micro=33554432
micro_benchmarks=("fake_test" "fake_test2")

for task in ${micro_benchmarks[@]}; do
    for net_config in ${net_config_list[@]}; do
        ip_address=${ip_addr_dict[${net_config}]}
        for assignment_strategy in  ${assignment_strategy_list[@]}; do
            keyword="${task}-${net_config}-${assignment_strategy}"
            bash ${root_folder}/scheduling/test_profiler.sh ${ip_address} ${ip_address} ${keyword} ${task} ${data_size_micro} ${assignment_strategy} ${net_config} ${parallelism_limit} -Micro
        done
    done
done

if [ ! -d "${root_folder}/scheduling/test" ]; then
    mkdir ${root_folder}/scheduling/test
fi

mv ${root_folder}/scheduling/Record_test ${root_folder}/scheduling/test/Record_${NETNAME}