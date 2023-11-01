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
    std::cout << "Calculating energy for " << std::string(filename) << "..." << std::endl;
    Params p = read_params(fin);
    std::map<int,double> H;
    for(int s_i = 0;s_i < p.n_s * p.n_T;++s_i){
        auto [T_i,lattice] = read_next(p.N, fin);
        double cur_T = p.T_min + T_i * (p.T_max - p.T_min) / (p.n_T - 1);
        double sum = 0.0;
        for(int i = 0;i < p.N;++i){
            for(int j = 0;j < p.N;++j){
                //Only take left and up to avoid double counting
                sum += lattice[i][j] * (lattice[(i-1+p.N)%p.N][j] + lattice[i][(j-1+p.N)%p.N]);
            }
        }
        H[T_i] += (-1)*sum / (p.N * p.N);   //We defined J = 1 initially
    }
    fin.close();
    std::cout << "Calculated energy data for " << std::string(filename) << ": " << std::endl;
    for(auto& [T_i,energy]:H)  energy /= p.n_s;
    for(auto [T_i,energy]:H){
        double cur_T = p.T_min + T_i * (p.T_max - p.T_min) / (p.n_T - 1);
        std::cout << cur_T << ' ' << energy << '\n';
    }
}

int main(int argc,char* argv[]){
    if(argc == 1){
        std::cout << "No arguments were passed. Please pass a filename to analyse.";
        return 0;
    }else{
        for(int argi = 1;argi < argc;++argi){
            analyze(argv[argi]);
        }
    }

    return 0;
}