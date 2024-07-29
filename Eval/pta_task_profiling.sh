USER_FOLDER=/root/aby3  
# task_list=("cipher_index" "max" "sort" "metric" "sum")
# task_list=("metric")
task_list=("cipher_index")
# task_list=("cipher_index" "metric")
# task_list=("sort" "metric" "sum")
# task_list=("cipher_index")

for task in ${task_list[@]}; do
    python ${USER_FOLDER}/PtA_deploy/system_profile.py --task_profile ${task}
done