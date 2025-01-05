root_folder=/root/aby3

NETWORK_NAME="Homo-100M"

${root_folder}/scheduling/network_setup.sh 100 100 100
wait;

${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


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


# NETWORK_NAME="Hetero-1G-11G-10G"

# ${root_folder}/scheduling/network_setup.sh 1050 1000 10000
# wait;

# ${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}