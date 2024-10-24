root_folder=/root/aby3

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
ip_address="10.3.0.15 10.3.0.16 10.3.0.17"
network_interface="ibs110 ibs110 ibs110"
fitting_length=16
fitting_step=128
complexity="1 n"
get_bandwidth_time=1
parallelism_limit=64

task_list=("Index" "Max" "Metric")
data_size=268435456

for task in ${task_list[@]}; do
    python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${task} \
    --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
    --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
    --run_tasks --MPI

    python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${task} \
    --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
    --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
    --run_tasks --MPI --skip_monitor --baseline
done

task_list=("Sort" "Matrix")
data_size=33554432

for task in ${task_list[@]}; do
    python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${task} \
    --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
    --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
    --run_tasks

    python ${root_folder}/scheduling/profiler.py --args " -${task}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${task} \
    --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} --network_interface ${network_interface} \
    --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
    --run_tasks --skip_monitor --baseline
done