USER_FOLDER=/home/tsingj_ubuntu/fanxy/PtA/aby3  
task_list=("cipher_index" "max" "sort")
# task_list=("cipher_index")

for task in ${task_list[@]}; do
    python ${USER_FOLDER}/PtA_deploy/system_profile.py --task_profile ${task}
done