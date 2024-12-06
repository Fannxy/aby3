args=" -LogReg -dataSize 100000 -p0_ip 10.3.0.13 -p1_ip 10.3.0.16 -rank 0"
keyword="LogReg"
root_folder=/root/aby3

# Sync the schedule
scp -r ./scheduling/*.py aby31:${root_folder}/scheduling/ &
scp -r ./scheduling/*.py aby32:${root_folder}/scheduling/ &
scp -r ./scheduling/*.sh aby31:${root_folder}/scheduling/ &
scp -r ./scheduling/*.sh aby32:${root_folder}/scheduling/ &
wait;

# # prepare the test cpp.
cp ${root_folder}/frontend/main.test ${root_folder}/frontend/main.cpp
python ${root_folder}/build.py;
scp -r ${root_folder}/out/build/linux/frontend/frontend aby31:${root_folder}/out/build/linux/frontend/ &
scp -r ${root_folder}/out/build/linux/frontend/frontend aby32:${root_folder}/out/build/linux/frontend/ &
wait;

# run the tests with monitor.
python ${root_folder}/scheduling/monitor_dis_run.py --keyword ${keyword} --record_folder ${root_folder}/scheduling/Record_test --role 0 --args "${args}" &
ssh aby31 "cd ./aby3/; python ${root_folder}/scheduling/monitor_dis_run.py --keyword ${keyword} --record_folder ./scheduling/Record_test --role 1 --args \"${args}\"" &
ssh aby32 "cd ./aby3/; python ${root_folder}/scheduling/monitor_dis_run.py --keyword ${keyword} --record_folder ./scheduling/Record_test --role 2 --args \"${args}\"" &
wait;

# analyze the records.
python ${root_folder}/scheduling/monitor_dis_run.py --record_folder ${root_folder}/scheduling/Record_test --keyword ${keyword} --role 0 --analysis;