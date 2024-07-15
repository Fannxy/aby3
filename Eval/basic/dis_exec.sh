args_list=$1
echo ${args_list}

# # synchronize with others
scp ./out/build/linux/frontend/frontend aby31:~/aby3/out/build/linux/frontend/ &
scp ./out/build/linux/frontend/frontend aby32:~/aby3/out/build/linux/frontend/ &
wait;

# ./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list} &
# wait;


./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list} &
ssh aby31 "cd ./aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list}" &
ssh aby32 "cd ./aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list}" &
wait;
