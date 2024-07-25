import os

if __name__ == "__main__":
    parames_folder = "/root/aby3/split/"
    input_file_prefix = parames_folder + "GCN_params"
    output_file = parames_folder + "parames" + "GCN_params" + ".txt"
    data_providers = 3
    params_len = 162
    
    args = f" -data_providers {data_providers} -input_file_prefix {input_file_prefix} -output_file {output_file} -params_len {params_len}"
    os.system(f"./Eval/dis_exec.sh \"{args}\"; wait;")