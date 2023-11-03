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
    //Calculate quotient and remainder of N by 8
    int qN8 = N / 8;
    int rN8 = N % 8;
    //Write temperature index
    fout.write(reinterpret_cast<const char*>(&T_i),sizeof(int));
    //Write lattice by mapping bools to blocks of 8 (1 byte = 8 bits)
    for(int i = 0;i < N;++i){
        for(int j = 0;j < qN8;++j){
            unsigned char sto = 0;
            for(int k = 0;k < 8;++k){
                sto |= ((lattice[i][8 * j + k] + 1) >> 1) << k;
            }
            fout.write(reinterpret_cast<const char*>(&sto),sizeof(sto));
        }
        //Possible remainder block
        if(rN8 > 0){
            unsigned char sto = 0;
            for(int k = 0;k < rN8;++k){
                sto |= ((lattice[i][8 * qN8 + k] + 1) >> 1) << k;
            }
            fout.write(reinterpret_cast<const char*>(&sto),sizeof(sto));
        }
    }
}