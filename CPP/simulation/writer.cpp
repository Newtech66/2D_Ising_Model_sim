#include "writer.h"
#include <iostream>

void write_params(const Params& sim_params, std::ofstream& fout){
    if(!fout.is_open()){
        std::cout << "File is not open!";
        exit(1);
    }
    //Write list of parameters
    fout.write(reinterpret_cast<const char*>(&sim_params.N),sizeof(int));
    fout.write(reinterpret_cast<const char*>(&sim_params.t_eq),sizeof(int));
    fout.write(reinterpret_cast<const char*>(&sim_params.t_a),sizeof(int));
    fout.write(reinterpret_cast<const char*>(&sim_params.n_s),sizeof(int));
    fout.write(reinterpret_cast<const char*>(&sim_params.T_min),sizeof(float));
    fout.write(reinterpret_cast<const char*>(&sim_params.T_max),sizeof(float));
    fout.write(reinterpret_cast<const char*>(&sim_params.n_T),sizeof(int));
}

void write_next(const int T_i, const int N, const Grid& lattice, std::ofstream& fout){
    if(!fout.is_open()){
        std::cout << "File is not open!";
        exit(1);
    }
    LatticeType flattened[N * N];
    for(int i = 0;i < N;++i){
        for(int j = 0;j < N;++j){
            flattened[N * i + j] = lattice[i][j];
        }
    }
    //Calculate quotient and remainder of N * N by 8
    int qNN8 = N * N / 8;
    int rNN8 = N * N % 8;
    //Write temperature index
    fout.write(reinterpret_cast<const char*>(&T_i),sizeof(int));
    //Write lattice by mapping bools to blocks of 8 (1 byte = 8 bits)
    for(int i = 0;i < N * N ;i += 8){
        unsigned char sto = 0;
        for(int k = 0;k < 8;++k){
            sto |= ((flattened[i + k] + 1) >> 1) << k;
        }
        fout.write(reinterpret_cast<const char*>(&sto),sizeof(sto));
    }
    //Possible remainder block
    if(rNN8 > 0){
        unsigned char sto = 0;
        for(int k = 0;k < rNN8;++k){
            sto |= ((flattened[8 * qNN8 + k] + 1) >> 1) << k;
        }
        fout.write(reinterpret_cast<const char*>(&sto),sizeof(sto));
    }
}