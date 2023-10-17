#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh.h"
//#include <chrono>

int main() {
  std::string coord_path = "coord_csv/";
  std::ofstream myfile;
  std::vector<std::vector<std::string>> pdb_chain_data;
  read_csv(pdb_chain_data, "pdb_12_with_hinges.csv");
  std::vector<std::pair<std::string, std::string>> pdb_pair_vec;
  for (const auto &pdb_chain : pdb_chain_data) {
    pdb_pair_vec.push_back(std::make_pair(pdb_chain[0], pdb_chain[1]));
  }
  myfile.open("rmsdhk_monge_check_result.csv");
  myfile << "p_pdb_id"
         << ","
         << "q_pdb_id"
         << ","
         << "monge_rate_1"
         << ","
         << "monge_rate_2"
         << ","
         << "monotonicity_rate_1"
         << ","
         << "monotonicity_rate_2" << std::endl;
  for (int i = 0; i < (int)pdb_pair_vec.size(); i++) {
    std::string p_pdb_id, q_pdb_id;
    std::string p_pdb_chain_id = pdb_pair_vec[i].first;
    std::string q_pdb_chain_id = pdb_pair_vec[i].second;
    std::transform(p_pdb_chain_id.begin(), p_pdb_chain_id.begin() + 4,
                   std::back_inserter(p_pdb_id), ::tolower);
    std::string p_chain_id = p_pdb_chain_id.substr(5, 1);
    std::transform(q_pdb_chain_id.begin(), q_pdb_chain_id.begin() + 4,
                   std::back_inserter(q_pdb_id), ::tolower);
    std::string q_chain_id = q_pdb_chain_id.substr(5, 1);
    Eigen::MatrixXd p, q;
    myfile << p_pdb_id << "," << q_pdb_id << ",";
    p = openMatrixData(coord_path + "coord_" + p_pdb_id + "_" + p_chain_id +
                       "_" + q_pdb_id + "_" + q_chain_id + ".csv");
    q = openMatrixData(coord_path + "coord_" + q_pdb_id + "_" + q_chain_id +
                       "_" + p_pdb_id + "_" + p_chain_id + ".csv");
    int total_residue_length = p.cols();
    std::cout << total_residue_length << std::endl;
    std::cout << p_pdb_id << ", " << q_pdb_id << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<std::vector<double>, std::vector<double>>
        monge_monotonicity_rate =
            CalcRMSDhKMongeCheck(p, q, total_residue_length, 2);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_time_ms = end - start;
    std::cout << exec_time_ms.count() << " ms" << std::endl;
    myfile << "," << monge_monotonicity_rate.first[0];
    myfile << "," << monge_monotonicity_rate.first[1];
    myfile << "," << monge_monotonicity_rate.second[0];
    myfile << "," << monge_monotonicity_rate.second[1] << std::endl;
  }
  myfile.close();
}
