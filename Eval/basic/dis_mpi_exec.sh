task_num=$1
args_list=$2
echo ${args_list}

USER_FOLDER=/root/aby3
# USER_FOLDER=/home/ubuntu/configuration/aby3

if [ $task_num -lt 256 ]; then
    max_core=$task_num
else
    max_core=256
fi

ulimit -n 65536;
taskset -c 0-255 mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 0 ${args_list} >> ./log 2>&1 &
ssh sosp1 "ulimit -n 65536; taskset -c 0-255 mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 1 ${args_list} &" &
ssh sosp2 "ulimit -n 65536; taskset -c 0-255 mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 2 ${args_list} &" &

wait;