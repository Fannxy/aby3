root_folder=/root/aby3
keyword="sort"
args=" -Sort"
data_size=1048576
fitting_length=16
fitting_step=64
fitting_degree=3
get_bandwidth_time=1
parallelism_limit=64

Sync the schedule
scp -r ./scheduling aby31:${root_folder}/ &
scp -r ./scheduling aby32:${root_folder}/ &
wait;

prepare the test cpp.
cp ${root_folder}/frontend/main.test ${root_folder}/frontend/main.cpp
python ${root_folder}/build.py
scp -r ${root_folder}/out/build/linux/frontend/frontend aby31:${root_folder}/out/build/linux/frontend/ &
scp -r ${root_folder}/out/build/linux/frontend/frontend aby32:${root_folder}/out/build/linux/frontend/ &
wait;

python ${root_folder}/scheduling/profiler.py --args "${args}" --record_folder ${root_folder}/scheduling/Record_test --keyword ${keyword} --data_size ${data_size} --fitting_length ${fitting_length} --fitting_step ${fitting_step} --fitting_degree ${fitting_degree} --get_bandwidth_time ${get_bandwidth_time} --parallelism_limit ${parallelism_limit} --run_tasks