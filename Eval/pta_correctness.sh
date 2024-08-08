test_args=" -pta_correctness"
cp ./frontend/main.pta ./frontend/main.cpp

rm ./debug.txt

python build.py

# Run the test
./Eval/basic/dis_mpi_exec.sh 4 "${test_args}"

cat ./debug.txt
rm ./debug.txt