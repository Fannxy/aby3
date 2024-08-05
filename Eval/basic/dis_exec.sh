args_list=$1
echo ${args_list}

USER_FOLDER=/root/aby3

# sync with others.
scp $USER_FOLDER/out/build/linux/frontend/frontend aby31:$USER_FOLDER/out/build/linux/frontend/
scp $USER_FOLDER/out/build/linux/frontend/frontend aby32:$USER_FOLDER/out/build/linux/frontend/
wait;

./out/build/linux/frontend/frontend -role 0 ${args_list} &
ssh aby31 "cd ./aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list}" &
ssh aby32 "cd ./aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list}" &
wait;


# ./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list} &
# ssh aby31 "cd ./aby3/; ./bin/frontend -prog -1 -role 1 ${args_list}" &
# ssh aby32 "cd ./aby3/; ./bin/frontend -prog -1 -role 2 ${args_list}" &
# wait;
