#!/bin/bash

download_and_extract() {
    local url=$1
    local output_dir=$2
    local file_name=$(basename "$url")
    local output_file="$output_dir/$file_name"

    # create folder
    mkdir -p "$output_dir"

    # download file
    wget -P "$output_dir" "$url"

    # extract file
    if [[ "$output_file" == *.gz ]]; then
        gzip -d "$output_file"
    fi
}

real_world_folder="./aby3-GORAM/data/real_world"
record_folder="./GORAM/results_real_world"
# mkdirs
if [ ! -d "$real_world_folder" ]; then
    mkdir -p "$real_world_folder"
else
    echo "Directory $real_world_folder already exists."
fi

if [ ! -d "$record_folder" ]; then
    mkdir -p "$record_folder"
else
    echo "Directory $record_folder already exists."
fi

echo -e "\e[32mSlashdot\e[0m"

# data prepare
if [ -f "$real_world_folder/soc-Slashdot0902.txt" ]; then
    echo "File $real_world_folder/soc-Slashdot0902.txt already exists."
else
    download_and_extract https://snap.stanford.edu/data/soc-Slashdot0902.txt.gz $real_world_folder/
fi
python ./aby3-GORAM/privGraphQuery/preprocessing.py --target slashdot
python ./GORAM/real_world_graph_analysis.py --target slashdot;

echo -e "\e[32mDBLP\e[0m"

if [ -f "$real_world_folder/com-dblp.all.cmty.txt" ]; then
    echo "File $real_world_folder/com-dblp.all.cmty.txt already exists."
else
    download_and_extract https://snap.stanford.edu/data/bigdata/communities/com-dblp.all.cmty.txt.gz $real_world_folder/
fi
python ./aby3-GORAM/privGraphQuery/preprocessing.py --target dblp
python ./GORAM/real_world_graph_analysis.py --target dblp;

cp -r ./aby3-GORAM/data/real_world/* $record_folder/


# large scale, requiring >3 * 400GB RAM, therefore we do not provide a full version here, the following code is for reference only.
# echo -e "\e[32mTwitter (it will take a long time...)\e[0m"

# if [ -f "$real_world_folder/twitter-2010.txt" ]; then
#     echo "File $real_world_folder/twitter-2010.txt already exists."
# else
#     download_and_extract https://snap.stanford.edu/data/twitter-2010.txt.gz $real_world_folder/
# fi
# python ./aby3-GORAM/privGraphQuery/preprocessing.py --target twitter

# # parallel processing.
# mpi_task=16
# python ./GORAM/mpi_data_organization.py --target twitter --MPI_TASK $mpi_task --origional_folder $real_world_folder/ --mpi_folder $real_world_folder/realworld_mpi/
# python ./GORAM/real_world_graph_analysis.py --target twitter --MPI true --MPI_TASK $mpi_task --data_folder $MAIN_FOLDER/realworld_mpi/
