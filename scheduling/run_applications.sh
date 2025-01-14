root_folder=/root/aby3

# rm -rf ${root_folder}/scheduling/Result

NETWORK_NAME="Homo-100M"
${root_folder}/scheduling/network_setup.sh 100 100 100
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

NETWORK_NAME="Homo-1G"
${root_folder}/scheduling/network_setup.sh 1050 1050 1050
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

NETWORK_NAME="Homo-10G"
${root_folder}/scheduling/network_setup.sh 10000 10000 10000
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-2G-2G"

${root_folder}/scheduling/network_setup.sh 1050 2000 2000
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-10G-10G"
${root_folder}/scheduling/network_setup.sh 1050 10000 10000
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}

NETWORK_NAME="Hetero-1G-4G-5G"
${root_folder}/scheduling/network_setup.sh 1050 4000 5000
wait;
${root_folder}/scheduling/benchmark_mininet.sh ${NETWORK_NAME}


mv ${root_folder}/scheduling/Result ${root_folder}/scheduling/Result_basic




NETWORK_NAME="Homo-100M"
${root_folder}/scheduling/network_setup.sh 100 100 100
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

NETWORK_NAME="Homo-1G"
${root_folder}/scheduling/network_setup.sh 1050 1050 1050
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

NETWORK_NAME="Homo-10G"
${root_folder}/scheduling/network_setup.sh 10000 10000 10000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-2G-2G"

${root_folder}/scheduling/network_setup.sh 1050 2000 2000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}


NETWORK_NAME="Hetero-1G-10G-10G"
${root_folder}/scheduling/network_setup.sh 1050 10000 10000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

NETWORK_NAME="Hetero-1G-4G-5G"
${root_folder}/scheduling/network_setup.sh 1050 4000 5000
wait;
${root_folder}/scheduling/benchmark_mininet_applications.sh ${NETWORK_NAME}

mv ${root_folder}/scheduling/Result ${root_folder}/scheduling/Result_application