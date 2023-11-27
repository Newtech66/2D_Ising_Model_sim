#include "reader.h"
#include <iostream>

Params read_params(std::ifstream& fin){
    /*
    *Read parameters from input file.
    *Input: ifstream to opened input file.
    */
    if(!fin.is_open()){
        std::cout << "File is not open!";
        exit(1);
    }
    Params p;
    fin.read(reinterpret_cast<char*>(&p.N),sizeof(p.N));
    fin.read(reinterpret_cast<char*>(&p.t_eq),sizeof(p.t_eq));
    fin.read(reinterpret_cast<char*>(&p.t_a),sizeof(p.t_a));
    fin.read(reinterpret_cast<char*>(&p.n_s),sizeof(p.n_s));
    fin.read(reinterpret_cast<char*>(&p.T_min),sizeof(p.T_min));
    fin.read(reinterpret_cast<char*>(&p.T_max),sizeof(p.T_max));
    fin.read(reinterpret_cast<char*>(&p.n_T),sizeof(p.n_T));
    return p;
}

std::pair<int,Grid> read_next(const int N, std::ifstream& fin){
    /*
    *Read (i, lattice) from input file.
    *Input: Lattice size N, ifstream to opened input file.
    */
    if(!fin.is_open()){
        std::cout << "File is not open!";
        exit(1);
    }
    LatticeType flattened[N * N]{};
    Grid lattice(N,Row(N));
    //read temp value index i
    int T_i;
    fin.read(reinterpret_cast<char*>(&T_i),sizeof(int));
    //Calculate quotient and remainder of N * N by 8
    int qNN8 = N * N / 8;
    int rNN8 = N * N % 8;
    //read snapshot
    for(int i = 0;i < N * N ;i += 8){
        unsigned char sto;
        fin.read(reinterpret_cast<char*>(&sto),sizeof(unsigned char));
        for(int k = 0;k < 8;++k){
            flattened[i + k] = 2 * ((sto & (1 << k)) >> k) - 1;
        }
    }
    //Possible remainder block
    if(rNN8 > 0){
        unsigned char sto;
        fin.read(reinterpret_cast<char*>(&sto),sizeof(unsigned char));
        for(int k = 0;k < rNN8;++k){
            flattened[8 * qNN8 + k] = 2 * ((sto & (1 << k)) >> k) - 1;
        }
    }
    //test if successfully read
    if(fin.eof()){
        std::cout << "No more snapshots to read!";
        exit(1);
    }
    //un-flatten array
    for(int i = 0;i < N;++i){
        for(int j = 0;j < N;++j){
            lattice[i][j] = flattened[N * i + j];
        }
    }
    return make_pair(T_i,lattice);
}