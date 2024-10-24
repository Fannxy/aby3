# mpi tasks
dataSize=33554432
task_list=("Index" "Max" "Metric")
symmetric=0
repeat=1

for task in "${task_list[@]}"; do
    ./scheduling/multitask_mpi_monitor.sh $dataSize $symmetric $repeat $task
done

symmetric=1
for task in "${task_list[@]}"; do
    ./scheduling/multitask_mpi_monitor.sh $dataSize $symmetric $repeat $task
done


non-mpi tasks
task_list=("Sort" "Matrix")
dataSize=1048576
symmetric=0
repeat=1

for task in "${task_list[@]}"; do
    ./scheduling/multitask_monitor.sh $dataSize $symmetric $repeat $task
done

symmetric=1
for task in "${task_list[@]}"; do
    ./scheduling/multitask_monitor.sh $dataSize $symmetric $repeat $task
done