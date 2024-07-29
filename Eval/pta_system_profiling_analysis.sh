cp ./frontend/main.pta ./frontend/main.cpp
USER_FOLDER=/root/aby3

times=1

for i in $(seq 1 $times); do
    python ${USER_FOLDER}/PtA_deploy/system_profile.py --system_profile
    mv ${USER_FOLDER}/system_profiling ${USER_FOLDER}/system_profiling_${i}
done

python ${USER_FOLDER}/PtA_deploy/system_analysis.py
mkdir ${USER_FOLDER}/system_profiling/Log
mv ${USER_FOLDER}/system_profiling_* ${USER_FOLDER}/system_profiling/Log/