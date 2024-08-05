task_num=$1
args_list=$2
echo ${args_list}

USER_FOLDER=/root/aby3

if [ $task_num -lt 24 ]; then
    max_core=$task_num
else
    max_core=24
fi

# sync with others.
scp $USER_FOLDER/out/build/linux/frontend/frontend aby31:$USER_FOLDER/out/build/linux/frontend/
scp $USER_FOLDER/out/build/linux/frontend/frontend aby32:$USER_FOLDER/out/build/linux/frontend/
wait;

ulimit -n 65536;
# taskset -c 0-255 mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 0 ${args_list} >> ./log 2>&1 &
# ssh aby31 "ulimit -n 65536; taskset -c 0-255 mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 1 ${args_list} &" &
# ssh aby32 "ulimit -n 65536; taskset -c 0-255 mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 2 ${args_list} &" &

mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 0 ${args_list} >> ./log 2>&1 &
ssh aby31 "ulimit -n 65536; mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 1 ${args_list} &" &
ssh aby32 "ulimit -n 65536; mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 2 ${args_list} &" &

wait;