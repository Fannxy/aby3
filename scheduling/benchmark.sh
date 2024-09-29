# mpi tasks
dataSize=268435456
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


# # TODO - non-mpi tasks
# task="Sort"
# dataSize=33554432
# symmetric=0
# repeat=1

# ./scheduling/multitask_monitor.sh $dataSize $task $symmetric $repeat
