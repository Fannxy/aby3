root_folder=/root/aby3
keyword="matrix"
args=" -Matrix"
# num_parties=4
# server_host="aby30 aby32 aby31 aby32"
# ip_address="10.3.0.15 10.5.0.17 10.3.0.16 10.3.0.17"
num_parties=3
server_host="aby30 aby31 aby32"
ip_address="10.3.0.15 10.3.0.16 10.3.0.17"
data_size=2097152
fitting_length=16
fitting_step=64
complexity="1 n"
get_bandwidth_time=1
parallelism_limit=64

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

python ${root_folder}/scheduling/profiler.py --args "${args}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${keyword} \
  --num_parties ${num_parties} --server_host ${server_host} --ip_address ${ip_address} \
  --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --complexity ${complexity} \
  --run_tasks \
#   --skip_monitor