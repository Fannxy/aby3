# compile the main.
cp ./frontend/main.test ./frontend/main.cpp
current_path=$(pwd)
debugFile="${current_path}/debug.txt"
# echo "Current path: ${current_path}"
python build.py --DEBUG_FILE ${debugFile}

# clean debugging files party-*.txt if exist.
for pfile in ./party-*.txt; do
    rm ${pfile};
done

# # synchronize with others
scp ./out/build/linux/frontend/frontend aby31:~/aby3/out/build/linux/frontend/ &
scp ./out/build/linux/frontend/frontend aby32:~/aby3/out/build/linux/frontend/ &
wait;

# run the tests
# current tests: 
# 1) -Bool : boolean share tests; 
# 2) -Arith : arithmetic share tests; 
# 3) -ORAM : ORAM tests; 
# 4) -Init : initialization tests, including the correlated shares; 
# 5) -Shuffle : secure shuffling tests.
# 8) -Comm : test inter-party communication.
# 9) -Sort : test the sort functions.
test_args=" -Bool -Arith -Comm -Sort -ORAM -Init -Shuffle"
./Eval/dis_exec.sh "${test_args}"
wait;

# scp aby31:~/aby3/debug.txt ./debug-p1.txt
# scp aby32:~/aby3/debug.txt ./debug-p2.txt

cat ./debug.txt
rm ./debug.txt
