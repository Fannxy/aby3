USER_FOLDER=/home/tsingj_ubuntu/fanxy/PtA/aby3

times=3

# for i in $(seq 1 $times); do
#     python ${USER_FOLDER}/PtA_deploy/task_execution.py --exhv
#     python ${USER_FOLDER}/PtA_deploy/task_analysis.py --exhv
#     mv ${USER_FOLDER}/Record ${USER_FOLDER}/Record_${i}
# done

python ${USER_FOLDER}/PtA_deploy/task_analysis.py --statistic
mkdir ${USER_FOLDER}/Record/Log
mv ${USER_FOLDER}/Record_* ${USER_FOLDER}/Record/Log/