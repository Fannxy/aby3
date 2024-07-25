#include <cryptoTools/Common/CLP.h>
#include "aby3-RTR/BuildingBlocks.h"
#include "aby3-GraphQuery/Graph.h"

using namespace oc;
using namespace aby3;

int main(int argc, char** argv)
{
    CLP cmd;
    cmd.parse(argc, argv);

    // std::string input_file;
    int data_providers = 1;
    if(cmd.isSet("data_providers"))
    {
        data_providers = cmd.get<int>("data_providers");
    }
    else{
        THROW_RUNTIME_ERROR("data_providers is not set");
    }
    std::string input_file_prefix;
    if(cmd.isSet("input_file_prefix"))
    {
        input_file_prefix = cmd.get<std::string>("input_file_prefix");
    }
    else{
        THROW_RUNTIME_ERROR("input_file_prefix is not set");
    }
    std::string output_file;
    if(cmd.isSet("output_file"))
    {
        output_file = cmd.get<std::string>("output_file");
    }
    else{
        THROW_RUNTIME_ERROR("output_file is not set");
    }
    std::vector<std::string> input_file_list;
    for(int i = 0; i < data_providers; i++)
    {
        input_file_list.push_back(input_file_prefix + "_" + std::to_string(i+1) + ".txt");
    }
    int param_len = 0;
    if(cmd.isSet("params_len")){
        param_len = cmd.get<int>("params_len");
    }
    else{
        THROW_RUNTIME_ERROR("param_len is not set");
    }

    // get the input and output file.

    // start the aby3 computation engine.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    aby3Info party_info(role, enc, eval, runtime);

    debug_info("before computation!");

    // p0 load the plaintext data.
    i64Matrix plain_params;

    if(role == 0){
        std::vector<std::vector<int64_t>> params(data_providers);
        for(int i=0; i<data_providers; i++)
        {
            std::ifstream input_file(input_file_list[i]);
            if(!input_file.is_open())
            {
                THROW_RUNTIME_ERROR("input file " + input_file_list[i] + " is not open");
            }
            int real_len;
            input_file >> real_len;
            if(real_len != param_len){
                THROW_RUNTIME_ERROR("real length != param length!");
            }

            for(int j=0; j<param_len; j++)
            {
                int64_t param;
                input_file >> param;
                params[i].push_back(param);
            }
            input_file.close();
        }

        int total_size = data_providers * param_len;
        plain_params.resize(total_size, 1);
        for(int i=0; i<data_providers; i++)
        {
            for(int j=0; j<param_len; j++)
            {
                plain_params(i*param_len+j, 0) = params[i][j];
            }
        }
    }

    debug_info("ok loading!");

    // generate the shares.
    si64Matrix enc_params(param_len * data_providers, 1);
    if(role == 0){
        debug_info("plain params size = " + std::to_string(plain_params.size()));
        party_info.enc->localIntMatrix(*(party_info.runtime), plain_params, enc_params).get();
    }
    else{
        party_info.enc->remoteIntMatrix(*(party_info.runtime), enc_params).get();
    }

    debug_info("ok generation!");
    
    // aggregate the parameters.
    si64Matrix agg_params(param_len, 1);
    // init_zeros(agg_params);
    for(int i=0; i<param_len; i++)
    {
        agg_params.mShares[0](i, 0) = 0;
        agg_params.mShares[1](i, 0) = 0;
        for(int j=0; j<data_providers; j++)
        {
            agg_params.mShares[0](i, 0) += enc_params.mShares[0](j*param_len+i, 0);
            agg_params.mShares[1](i, 0) += enc_params.mShares[1](j*param_len+i, 0);
        }
    }

    debug_info("ok aggregation!");

    // all the parties reveal the shares to P0.
    i64Matrix reveal_params(param_len, 1);
    party_info.enc->revealAll(*(party_info.runtime), agg_params, reveal_params).get();

    debug_info("ok reveal!!");

    // P0 saves the result to the output file.
    if(role == 0){
        std::ofstream output_file_fs(output_file);
        if(!output_file_fs.is_open())
        {
            THROW_RUNTIME_ERROR("output file " + output_file + " is not open");
        }
        output_file_fs << param_len << std::endl;
        for(int i=0; i<param_len; i++)
        {
            output_file_fs << reveal_params(i, 0) << std::endl;
        }
        output_file_fs.close();
    }

    return 0;   
}