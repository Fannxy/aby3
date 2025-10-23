cp ./frontend/main.test ./frontend/main.cpp
python build.py 

./out/build/linux/frontend/frontend  -role 0 &
./out/build/linux/frontend/frontend  -role 1 &
./out/build/linux/frontend/frontend  -role 2 &

wait;