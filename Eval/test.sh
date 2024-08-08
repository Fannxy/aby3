# prepare the test cpp.
cp ./frontend/main.test ./frontend/main.cpp
python build.py

# clean debugging files party-*.txt if exist.
for pfile in ./party-*.txt; do
    rm ${pfile};
done

# # synchronize with others
# scp ./bin/frontend aby31:~/aby3/bin
# scp ./bin/frontend aby32:~/aby3/bin
# wait;

# run the tests
# current tests: 
# 1) -Bool : boolean share tests; 
# 2) -Arith : arithmetic share tests; 
# 3) -ORAM : ORAM tests; 
# 4) -Init : initialization tests, including the correlated shares; 
# 5) -Shuffle : secure shuffling tests.
# 6) -Graph : basic graph loading tests.
# 7) -GraphQuery : basic graph query tests (block fetching, edge exist & outting edges count.)
# 8) -Comm : test inter-party communication.
# 9) -Sort : test the sort functions.
# test_args=" -ORAM -Graph"
# test_args=" -Sort -Search"
test_args=" -MPI-Search"
# ./Eval/basic/dis_exec.sh "${test_args}"
./Eval/basic/dis_mpi_exec.sh 4 "${test_args}"
wait;

# scp aby31:~/aby3/debug.txt ./debug-p1.txt
# scp aby32:~/aby3/debug.txt ./debug-p2.txt

cat ./debug.txt
rm ./debug.txt
