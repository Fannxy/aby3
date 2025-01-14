#include "GASTest.h"
#include <chrono>
#include <thread>
#include "./Pair_then_Reduce/include/datatype.h"
#include "BuildingBlocks.h"
#include "GASFunction.h"
#include "PTRFunction.h"

using namespace oc;
using namespace aby3;
using namespace std;

int get_value(string key, CLP& cmd, int& val) {
  if (cmd.isSet(key)) {
    auto keys = cmd.getMany<int>(key);
    val = keys[0];
  } else {
    throw std::runtime_error(LOCATION);
  }
  return 0;
}

int get_value(string key, CLP& cmd, string& val) {
  if (cmd.isSet(key)) {
    auto keys = cmd.getMany<string>(key);
    if(keys.size() == 0) cout << key << " value is null" << endl;
    val = keys[0];
  }
  return 0;
}

int get_value(string key, CLP& cmd, size_t& val) {
  if (cmd.isSet(key)) {
    auto keys = cmd.getMany<size_t>(key);
    if(keys.size() == 0) cout << key << "value is null" << endl;
    val = keys[0];
  }
  return 0;
}

int get_mean_and_std(const std::vector<double>& time_list, double& mean, double& std_var){
  mean = 0;
  std_var = 0;
  size_t list_len = time_list.size();

  for(int i=0; i<list_len; i++) mean += (time_list[i]); 
  mean /= list_len;

  for(int i=0; i<list_len; i++) std_var += ((time_list[i] - mean) * (time_list[i] - mean));
  std_var /= (list_len - 1);

  return 0;
}

int dict_to_file(const std::map<string, double>& time_dict, std::string file_path){
  std::ofstream ofs(file_path, std::ios_base::app);
  for(auto it = time_dict.begin(); it != time_dict.end(); ++it){
    ofs << it->first << ": " << it->second << "\n";
  }
  ofs.close();
  return 0;
}


int test_page_rank(CLP& cmd, const std::string& data_folder) {
  // setup the environment
  int role, nodes_num, edges_num, iters;
  size_t vector_size;
  role = -1;
  string log_file, log_folder;

  get_value("role", cmd, role);
  get_value("nodes_num", cmd, nodes_num);
  get_value("edges_num", cmd, edges_num);
  get_value("OPT_BLOCK", cmd, vector_size);
  get_value("iters", cmd, iters);
  get_value("logFile", cmd, log_file);
  get_value("logFolder", cmd, log_folder);

  clock_t start, end;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  log_file = log_folder + log_file + "E=" + to_string(edges_num) + "-V=" + to_string(nodes_num) + "-TASKS=" + to_string(size) + "-OPT_B" + to_string(vector_size) + "-iters=" + to_string(iters);

  start = clock();
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_setup = double((end - start) * 1000) / (CLOCKS_PER_SEC);

  // load the graph data
  start = clock();
  std::vector<sfnode<D8>> cipher_nodes;
  std::vector<sfedge<D8>> cipher_edges;
  load_nodes_data<D8>(role, enc, runtime, data_folder + "/nodes.txt",
                      cipher_nodes, 0, nodes_num);
  load_edges_data<D8>(role, enc, runtime, data_folder + "/edges.txt", cipher_edges, 0, edges_num);

  vector<sf64<D8>> initial_prs;
  for(int i=0; i<nodes_num; i++){
    initial_prs.push_back(cipher_nodes[i].pr);
  }
  end = clock();
  double time_data_loading = double((end - start) * 1000) / (CLOCKS_PER_SEC);

  // init the PtA class.
  auto mpiScatter = new MPIPRScatter<sfedge<D8>, sfnode<D8>, aby3::sf64<D8>, aby3::sf64<D8>, SubPRScatterABY3>(size, vector_size, role, enc, runtime, eval);
  auto mpiGather = new MPIPRScatter<sfnode<D8>, sfedge<D8>, aby3::sf64<D8>, aby3::sf64<D8>, SubPRGatherABY3>(size, vector_size, role, enc, runtime, eval);


  aby3::sf64<D8> dval;
  dval.mShare.mData[0] = 0, dval.mShare.mData[1] = 0;
  mpiScatter->set_default_value(dval); mpiGather->set_default_value(dval);

  mpiScatter->circuit_construct({(size_t) edges_num}, {(size_t) nodes_num});
  mpiGather->circuit_construct({(size_t) nodes_num}, {(size_t) edges_num});

  vector<double> gather_time_list(iters, 0.0);
  vector<double> scatter_time_list(iters, 0.0);
  vector<double> middle_edges_sharing(iters, 0.0);
  vector<double> middle_nodes_sharing(iters - 1, 0.0);

  start = clock();
  int iter = 0;
  while(iter < iters){
    // scatter stage.
    clock_t mid_start, mid_end;
    mid_start = clock();
    vector<sf64<D8>> edges_data; init_zeros(role, enc, runtime, edges_data, edges_num);
    mpiScatter->circuit_evaluate(cipher_edges.data(), cipher_nodes.data(), nullptr, edges_data.data());
    mid_end = clock();
    double time_scatter = double((mid_end - mid_start) * 1000) / (CLOCKS_PER_SEC);
    
    mid_start = clock();
    // data sharing and synchronize.
    sf64<D8>* edges_data_ptr = edges_data.data();
    mpiScatter->data_sharing<sf64<D8>>(edges_data_ptr, edges_num, 0);
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i=0; i<edges_num; i++) cipher_edges[i].data = edges_data[i];
    mid_end = clock();
    double time_edge_sharing = double((mid_end - mid_start) * 1000) / (CLOCKS_PER_SEC);

    // gather stage
    mid_start = clock();
    vector<sf64<D8>> nodes_data; init_zeros(role, enc, runtime, nodes_data, nodes_num);
    mpiGather->circuit_evaluate(cipher_nodes.data(), cipher_edges.data(), nullptr, nodes_data.data());
    mid_end = clock();
    double time_gather = double((mid_end - mid_start) * 1000) / (CLOCKS_PER_SEC);

    mid_start = clock();
    if(iter < (iters - 1)){
      // data sharing and synchonize.
      sf64<D8>* nodess_data_ptr = nodes_data.data();
      mpiScatter->data_sharing<sf64<D8>>(nodess_data_ptr, nodes_num, 0);
      MPI_Barrier(MPI_COMM_WORLD);
      for(int i=0; i<nodes_num; i++) cipher_nodes[i].pr = nodes_data[i];
    }
    mid_end = clock();
    double time_node_sharing = double((mid_end - mid_start) * 1000) / (CLOCKS_PER_SEC);

    gather_time_list[iter] = time_gather;
    scatter_time_list[iter] = time_scatter;
    middle_edges_sharing[iter] = time_edge_sharing;
    if(iter < (iters - 1)) middle_nodes_sharing[iter] = time_node_sharing;

    iter += 1;
  }
  end = clock();
  double time_eval = double((end - start) * 1000) / (CLOCKS_PER_SEC);

  map<string, double> time_dict;
  time_dict["time_setup"] = time_setup;
  time_dict["time_data_prepare"] = time_data_loading;
  time_dict["time_eval"] = time_eval;
  double mean_gather, mean_scatter, mean_edge_sharing, mean_node_sharing;
  double std_gather, std_scatter, std_edge_sharing, std_node_sharing;
  get_mean_and_std(gather_time_list, mean_gather, std_gather);
  get_mean_and_std(scatter_time_list, mean_scatter, std_scatter);
  get_mean_and_std(middle_edges_sharing, mean_edge_sharing, std_gather);
  get_mean_and_std(middle_nodes_sharing, mean_node_sharing, std_node_sharing);
  time_dict["time_mean_iter"] = mean_gather + mean_scatter + mean_edge_sharing + mean_node_sharing;
  time_dict["time_mean_gather"] = mean_gather;
  time_dict["time_mean_scatter"] = mean_scatter;
  time_dict["time_mean_edge_sharing"] = mean_edge_sharing;
  time_dict["time_mean_node_sharing"] = mean_node_sharing;

  if(rank == 0){
    dict_to_file(time_dict, log_file);
  }
  return 0;
}