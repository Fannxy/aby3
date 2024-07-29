USER_FOLDER=/root/aby3

times=5

for i in $(seq 1 $times); do
    python ${USER_FOLDER}/PtA_deploy/task_execution.py --config
    python ${USER_FOLDER}/PtA_deploy/task_analysis.py --config
    mv ${USER_FOLDER}/Record ${USER_FOLDER}/Record_${i}
done

python ${USER_FOLDER}/PtA_deploy/task_analysis.py --statistic
mkdir ${USER_FOLDER}/Record/Log
mv ${USER_FOLDER}/Record_* ${USER_FOLDER}/Record/Log/