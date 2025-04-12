#!/bin/bash

NETNAME=$1
echo "NETNAME >>> " $NETNAME
root_folder=/root/aby3

num_parties=3
# node_id=(12 14 4)
node_id=(11 12 13)
server_host="aby30 aby31 aby32"
parallelism_limit=96
min_parallelism=16

# prepare the test cpp.
cp ${root_folder}/frontend/main.test ${root_folder}/frontend/main.cpp
python ${root_folder}/build.py
wait;

declare -A ip_addr_dict
ip_addr_dict[${NETNAME}]="10.0.0.${node_id[0]} 10.0.0.${node_id[1]} 10.0.0.${node_id[2]}"


assignment_strategy_list=("roundrole" "baseline")

net_config_list=(${NETNAME})
data_size_micro=1416810830
if [ $NETNAME = "Homo-10G" ]; then
    data_size_micro=1416810830
fi
if [[ $NETNAME == "Hetero-"* ]]; then
    data_size_micro=1416810830
fi

# micro_benchmarks=("ff-mul" "ib-mul" "a2b" "b2a-single" "shuffle")
micro_benchmarks=("ff-mul" "ib-mul")
# micro_benchmarks=("shuffle")
for task in ${micro_benchmarks[@]}; do
    for net_config in ${net_config_list[@]}; do
        ip_address=${ip_addr_dict[${net_config}]}
        for assignment_strategy in  ${assignment_strategy_list[@]}; do
            # for fix_strategy in ${fix_balance_list[@]}; do
            keyword="${task}-${net_config}-${assignment_strategy}"
            bash ${root_folder}/scheduling/test_profiler.sh ${ip_address} ${ip_address} ${keyword} ${task} ${data_size_micro} ${assignment_strategy} ${net_config} ${parallelism_limit} True ${min_parallelism} -Micro
            # done
        done
    done
done

DIR="/ssdshare/fanxy/roundrole1/Result"

if [ ! -d "$DIR" ]; then
  mkdir -p "$DIR"
  echo "目录 $DIR 已创建。"
else
  echo "目录 $DIR 已存在。"
fi

cp -r ${root_folder}/scheduling/Record_test /ssdshare/fanxy/roundrole1/Result
mv /ssdshare/fanxy/roundrole1/Result/Record_test /ssdshare/fanxy/roundrole1/Result/Record_${NETNAME}
rm -rf ${root_folder}/scheduling/Record_test