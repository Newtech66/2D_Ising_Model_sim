#include "reader.h"
#include <iostream>
#include <string>
#include <map>

void analyze(const char* filename){
    std::ifstream fin(filename,std::ios_base::binary);
    if(!fin.is_open()){
        std::cout << "Failed to open " << std::string(filename) << '\n';
        return;
    }
    std::cout << "Checking data for " << std::string(filename) << "..." << std::endl;
    Params p = read_params(fin);
    std::cout << "N, t_eq, t_a: " << p.N << ' ' << p.t_eq << ' ' << p.t_a << std::endl;
    std::cout << "n_s, T_min, T_max, n_T: " << p.n_s << ' ' << p.T_min << ' ' << p.T_max << ' ' << p.n_T << std::endl;
    std::map<int,int> freq;
    for(int s_i = 0;s_i < p.n_s * p.n_T;++s_i){
        auto [T_i,lattice] = read_next(p.N, fin);
        freq[T_i]++;
    }
    fin.close();
    int sum = 0;
    for(auto [T_i,cnt]:freq){
        float cur_T = p.T_min + T_i * (p.T_max - p.T_min) / (p.n_T - 1);
        std::cout << "Number of snapshots for temperature " << cur_T << ": " << cnt << '\n';
        sum += cnt;
    }
    std::cout << "Total number of snapshots: " << sum << std::endl;
}

int main(int argc,char* argv[]){
    if(argc == 1){
        std::cout << "No arguments were passed. Please pass one or more filenames to analyze.";
        return 0;
    }else if(argc > 2){
        for(int argi = 1;argi < argc;++argi){
            analyze(argv[argi]);
        }
    }
    return 0;
}