root_folder=/root/aby3

# rm -rf ${root_folder}/scheduling/Result

NETWORK_NAME="Homo-100M-100ms"
${root_folder}/scheduling/network_setup.sh 100 100 100 25ms 25ms 25ms
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

NETWORK_NAME="Homo-100M-50ms"
${root_folder}/scheduling/network_setup.sh 100 100 100 12.5ms 12.5ms 12.5ms
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

NETWORK_NAME="Homo-100M-10ms"
${root_folder}/scheduling/network_setup.sh 100 100 100 2.5ms 2.5ms 2.5ms
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


NETWORK_NAME="Homo-100M-1ms"
${root_folder}/scheduling/network_setup.sh 100 100 100 0.25ms 0.25ms 0.25ms
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


NETWORK_NAME="Homo-100M-0.1ms"
${root_folder}/scheduling/network_setup.sh 100 100 100 0.1us 0.1us 0.1us
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# mv /ssdshare/fanxy/roundrole1/Result /ssdshare/fanxy/roundrole1/Result_basic

# NETWORK_NAME="Homo-10G"
# ${root_folder}/scheduling/network_setup.sh 10000 10000 10000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


# NETWORK_NAME="Hetero-1G-2G-2G-1"

# ${root_folder}/scheduling/network_setup.sh 1050 2000 2000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# NETWORK_NAME="Hetero-1G-2G-2G-2"

# ${root_folder}/scheduling/network_setup.sh 2000 1050 2000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


# NETWORK_NAME="Hetero-1G-2G-2G-3"

# ${root_folder}/scheduling/network_setup.sh 2000 2000 1050
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


# NETWORK_NAME="Hetero-1G-10G-10G-1"
# ${root_folder}/scheduling/network_setup.sh 1050 10000 10000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# NETWORK_NAME="Hetero-1G-10G-10G-2"
# ${root_folder}/scheduling/network_setup.sh 1050 10000 10000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# # NETWORK_NAME="Hetero-1G-10G-10G-3"
# # ${root_folder}/scheduling/network_setup.sh 10000 10000 1050
# # wait;
# # ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


# NETWORK_NAME="Hetero-1G-4G-5G-1"
# ${root_folder}/scheduling/network_setup.sh 1050 4000 5000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# NETWORK_NAME="Hetero-1G-4G-5G-2"
# ${root_folder}/scheduling/network_setup.sh 1050 5000 4000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# NETWORK_NAME="Hetero-1G-4G-5G-3"
# ${root_folder}/scheduling/network_setup.sh 4000 1050 5000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# NETWORK_NAME="Hetero-1G-4G-5G-4"
# ${root_folder}/scheduling/network_setup.sh 4000 5000 1050
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# NETWORK_NAME="Hetero-1G-4G-5G-5"
# ${root_folder}/scheduling/network_setup.sh 5000 4000 1050
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

# NETWORK_NAME="Hetero-1G-4G-5G-6"
# ${root_folder}/scheduling/network_setup.sh 5000 1050 4000
# wait;
# ${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


# mv ${root_folder}/scheduling/Result ${root_folder}/scheduling/Result_basic



# applications!!!!!

NETWORK_NAME="Hetero-1G-2G-2G-1"

${root_folder}/scheduling/network_setup.sh 1050 2000 2000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-2G-2G-2"

${root_folder}/scheduling/network_setup.sh 2000 1050 2000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-2G-2G-3"

${root_folder}/scheduling/network_setup.sh 2000 2000 1050
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-10G-10G-1"
${root_folder}/scheduling/network_setup.sh 1050 10000 10000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-10G-10G-2"
${root_folder}/scheduling/network_setup.sh 10000 1050 10000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-10G-10G-3"
${root_folder}/scheduling/network_setup.sh 10000 10000 1050
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-4G-5G-1"
${root_folder}/scheduling/network_setup.sh 1050 4000 5000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

NETWORK_NAME="Hetero-1G-4G-5G-2"
${root_folder}/scheduling/network_setup.sh 1050 5000 4000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

NETWORK_NAME="Hetero-1G-4G-5G-3"
${root_folder}/scheduling/network_setup.sh 4000 1050 5000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-4G-5G-4"
${root_folder}/scheduling/network_setup.sh 4000 5000 1050
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

NETWORK_NAME="Hetero-1G-4G-5G-5"
${root_folder}/scheduling/network_setup.sh 5000 4000 1050
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

NETWORK_NAME="Hetero-1G-4G-5G-6"
${root_folder}/scheduling/network_setup.sh 5000 1050 4000
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


# NETWORK_NAME="Homo-100M"
# ${root_folder}/scheduling/network_setup.sh 100 100 100
# wait;
# ${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

# NETWORK_NAME="Homo-1G"
# ${root_folder}/scheduling/network_setup.sh 1050 1050 1050
# wait;
# ${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

# NETWORK_NAME="Homo-10G"
# ${root_folder}/scheduling/network_setup.sh 10000 10000 10000
# wait;
# ${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


# NETWORK_NAME="Hetero-1G-2G-2G"

# ${root_folder}/scheduling/network_setup.sh 1050 2000 2000
# wait;
# ${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


# NETWORK_NAME="Hetero-1G-10G-10G"
# ${root_folder}/scheduling/network_setup.sh 1050 10000 10000
# wait;
# ${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


# NETWORK_NAME="Hetero-1G-4G-5G"
# ${root_folder}/scheduling/network_setup.sh 1050 4000 5000
# wait;
# ${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

mv ${root_folder}/scheduling/Result ${root_folder}/scheduling/Result_application