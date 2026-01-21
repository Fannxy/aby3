#!/bin/bash

# 清理可能存在的旧 frontend 进程，避免端口冲突
pkill -f "frontend.*-role" || true
sleep 1  # 等待进程完全退出

# 测量 build 时间
BUILD_START=$(date +%s%N)
cp ./frontend/main.test ./frontend/main.cpp
python build.py 
BUILD_END=$(date +%s%N)
BUILD_DURATION=$(awk "BEGIN {printf \"%.3f\", ($BUILD_END - $BUILD_START) / 1000000000}")  # 转换为秒，保留3位小数

# 测量运行时间
RUN_START=$(date +%s%N)
./out/build/linux/frontend/frontend  -role 0 &
./out/build/linux/frontend/frontend  -role 1 &
./out/build/linux/frontend/frontend  -role 2 &
wait
RUN_END=$(date +%s%N)
RUN_DURATION=$(awk "BEGIN {printf \"%.3f\", ($RUN_END - $RUN_START) / 1000000000}")  # 转换为秒，保留3位小数

# 清理 frontend 进程
pkill -f "frontend.*-role" || true

# 写入时间结果到文件（如果环境变量设置了）
if [ ! -z "$FEDDB_ABY3_TIME_FILE" ]; then
    cat > "$FEDDB_ABY3_TIME_FILE" << EOF
Aby3 Dispatch Time Measurement
===============================
Build Time: $BUILD_DURATION s
Run Time: $RUN_DURATION s
EOF
fi