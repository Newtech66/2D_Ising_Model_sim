#include "reader.h"
#include <iostream>
#include <string>
#include <map>
#include <cmath>

int main(int argc,char* argv[]){
    if(argc == 1){
        std::cout << "No arguments were passed. Please pass a filename to analyse.";
        return 0;
    }else if(argc > 2){
        std::cout << "More than one argument was passed. Only the first one will be kept.\n";
    }
    std::ifstream fin(argv[1],std::ios_base::binary);
    if(!fin.is_open()){
        std::cout << "Failed to open " << std::string(argv[1]) << '\n';
        return 0;
    }
    Params p = read_params(fin);
    std::map<float,float> magnetisation;
    for(int s_i = 0;s_i < p.n_s * p.n_T;++s_i){
        auto [T_i,lattice] = read_next(p.N, fin);
        float cur_T = p.T_min + T_i * (p.T_max - p.T_min) / (p.n_T - 1);
        float sum = 0.0;
        for(int i = 0;i < p.N;++i){
            for(int j = 0;j < p.N;++j){
                sum += lattice[i][j];
            }
        }
        sum /= (p.N * p.N);
        magnetisation[cur_T] += std::abs(sum);
    }
    fin.close();
    for(auto [T,m]:magnetisation){
        std::cout << T << ' ' << m / p.n_s << '\n';
    }
    return 0;
}