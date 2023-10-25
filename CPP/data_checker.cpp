#include "reader.h"
#include <iostream>
#include <string>
#include <map>

int main(int argc,char* argv[]){
    if(argc == 1){
        std::cout << "No arguments were passed. Please pass a filename to analyse.";
        return 0;
    }else if(argc > 2){
        std::cout << "More than one arguments were passed. Only first one will be kept.\n";
    }
    std::ifstream fin(argv[1],std::ios_base::binary);
    if(!fin.is_open()){
        std::cout << "Failed to open " << std::string(argv[1]) << '\n';
        return 0;
    }
    Params p = read_params(fin);
    std::cout << p.N << ' ' << p.t_eq << ' ' << p.t_a << ' ';
    std::cout << p.n_s << ' ' << p.T_min << ' ' << p.T_max << ' ';
    std::cout << p.n_T << '\n';
    std::map<int,int> freq;
    for(int s_i = 0;s_i < p.n_s * p.n_T;++s_i){
        auto [T_i,lattice] = read_next(p.N, fin);
        freq[T_i]++;
    }
    fin.close();
    int sum = 0;
    for(auto [T_i,cnt]:freq){
        std::cout << T_i << ": " << cnt << '\n';
        sum += cnt;
    }
    std::cout << sum;
    return 0;
}