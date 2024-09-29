dataSize=$1
symmetric=$2
repeat=$3
task=$4 # possible options: "Index, Max, Metric"

args=" -"${task}" -dataSize "${dataSize}
keyword=${task}"-"${symmetric}"-"${dataSize}
if [ ${symmetric} -eq 0 ]; then
    options=" --order 0,1,2 --order 0,1,2 --order 0,1,2 --repeat "${repeat}
else 
    options=" --order 0,1,2 --order 1,2,0 --order 2,0,1 --repeat "${repeat}
fi
root_folder=/root/aby3

# Sync the schedule
scp -r ./scheduling aby31:${root_folder}/ &
scp -r ./scheduling aby32:${root_folder}/ &
wait;

rm -r ./debug.txt

# prepare the test cpp.
cp ${root_folder}/frontend/main.test ${root_folder}/frontend/main.cpp
python ${root_folder}/build.py
scp -r ${root_folder}/out/build/linux/frontend/frontend aby31:${root_folder}/out/build/linux/frontend/ &
scp -r ${root_folder}/out/build/linux/frontend/frontend aby32:${root_folder}/out/build/linux/frontend/ &
wait;

# run the tests with monitor.
python ${root_folder}/scheduling/monitor_dis_run_multitask.py --keyword ${keyword} --record_folder ${root_folder}/scheduling/Record_test --role 0 --args "${args}" ${options} --MPI &
ssh aby31 "cd ./aby3/; python ${root_folder}/scheduling/monitor_dis_run_multitask.py --keyword ${keyword} --record_folder ./scheduling/Record_test --role 1 --args \"${args}\" ${options} --MPI" &
ssh aby32 "cd ./aby3/; python ${root_folder}/scheduling/monitor_dis_run_multitask.py --keyword ${keyword} --record_folder ./scheduling/Record_test --role 2 --args \"${args}\" ${options} --MPI" &
wait;

# analyze the records.
python ${root_folder}/scheduling/monitor_dis_run.py --record_folder ${root_folder}/scheduling/Record_test --keyword ${keyword} --role 0 --analysis;