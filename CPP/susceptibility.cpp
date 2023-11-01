#include "reader.h"
#include <iostream>
#include <string>
#include <map>
#include <cmath>

void analyze(const char* filename){
    std::ifstream fin(filename,std::ios_base::binary);
    if(!fin.is_open()){
        std::cout << "Failed to open " << std::string(filename) << '\n';
        return;
    }
    std::cout << "Calculating susceptibility for " << std::string(filename) << "..." << std::endl;
    Params p = read_params(fin);
    std::map<int,double> tot, sq_tot;
    for(int s_i = 0;s_i < p.n_s * p.n_T;++s_i){
        auto [T_i,lattice] = read_next(p.N, fin);
        double cur_T = p.T_min + T_i * (p.T_max - p.T_min) / (p.n_T - 1);
        double sum = 0.0;
        for(int i = 0;i < p.N;++i){
            for(int j = 0;j < p.N;++j){
                sum += lattice[i][j];
            }
        }
        tot[T_i] += sum;
        sq_tot[T_i] += sum * sum;
    }
    fin.close();
    std::cout << "Calculated susceptibility data for " << std::string(filename) << ": " << std::endl;
    for(auto& [T_i,S]:tot)  S /= p.n_s;
    for(auto& [T_i,S]:sq_tot)  S /= p.n_s;
    for(auto [T_i,S1]:sq_tot){
        double cur_T = p.T_min + T_i * (p.T_max - p.T_min) / (p.n_T - 1);
        double S2 = tot[T_i];
        double X = (S1 - S2 * S2) / (2 * cur_T); //We defined k = 2 initially
        std::cout << cur_T << ' ' << X / (p.N * p.N) << '\n';
    }
}

int main(int argc,char* argv[]){
    if(argc == 1){
        std::cout << "No arguments were passed. Please pass a filename to analyse.";
        return 0;
    }else if(argc > 2){
        for(int argi = 1;argi < argc;++argi){
            analyze(argv[argi]);
        }
    }

    return 0;
}