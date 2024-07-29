args=" -Matrix -sizeX 200 -sizeY 200 -sizeZ 200 -numTasks 3"
keyword="matrix"
# options=" --symmetric"
options=""
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
scp -r ${root_folder}/scheduling/monitor_dis_run_multitask.py aby31:${root_folder}/scheduling/monitor_dis_run_multitask.py &
scp -r ${root_folder}/scheduling/monitor_dis_run_multitask.py aby32:${root_folder}/scheduling/monitor_dis_run_multitask.py &
wait;

# run the tests with monitor.
python ${root_folder}/scheduling/monitor_dis_run_multitask.py --keyword ${keyword} --record_folder ${root_folder}/scheduling/Record_test --role 0 --args "${args}" ${options} &
ssh aby31 "cd ./aby3/; python ${root_folder}/scheduling/monitor_dis_run_multitask.py --keyword ${keyword} --record_folder ./scheduling/Record_test --role 1 --args \"${args}\" ${options}" &
ssh aby32 "cd ./aby3/; python ${root_folder}/scheduling/monitor_dis_run_multitask.py --keyword ${keyword} --record_folder ./scheduling/Record_test --role 2 --args \"${args}\" ${options}" &
wait;

# analyze the records.
python ${root_folder}/scheduling/monitor_dis_run.py --record_folder ${root_folder}/scheduling/Record_test --keyword ${keyword} --role 0 --analysis;