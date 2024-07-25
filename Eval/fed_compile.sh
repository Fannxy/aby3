# compile the file.
cp ./frontend/main.fed ./frontend/main.cpp
python build.py

# # synchronize with others
scp ./out/build/linux/frontend/frontend aby31:~/aby3/out/build/linux/frontend/ &
scp ./out/build/linux/frontend/frontend aby32:~/aby3/out/build/linux/frontend/ &
wait;

python ./parameter_aggregation.py