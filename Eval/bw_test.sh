# 带宽监测校验：仅运行 aby3 三方大量传输（约 100MB），用于验证网络带宽监测
# 用法：在 aby3 根目录执行 ./Eval/bw_test.sh

# 复制 bw.test 到 main.cpp
cp ./frontend/bw.test ./frontend/main.cpp

current_path=$(pwd)
debugFile="${current_path}/debug.txt"
graphFolder="${current_path}/aby3-GORAM/data/"
echo "Current path: ${current_path}"
python build.py --DEBUG_FILE ${debugFile} --GRAPH_FOLDER ${graphFolder}

# 运行三方带宽测试（带监测）
test_args=" -BwTest"
./Eval/dis_exec.sh "${test_args}"
wait

if [ -f ./debug.txt ]; then
  cat ./debug.txt
  rm ./debug.txt
fi
