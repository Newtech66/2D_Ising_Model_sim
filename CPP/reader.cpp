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
    std::pair<int,Grid> snapshot;
    //read temp value index i
    int T_i;
    fin.read(reinterpret_cast<char*>(&T_i),sizeof(T_i));
    //read snapshot
    int qN8 = N / 8;
    int rN8 = N % 8;
    Grid lattice(N,Row(N));
    for(int i = 0;i < N;++i){
        for(int j = 0;j < qN8;++j){
            unsigned char sto;
            fin.read(reinterpret_cast<char*>(&sto),sizeof(sto));
            for(int k = 0;k < 8;++k){
                lattice[i][8 * j + k] = 2 * ((sto & (1 << k)) >> k) - 1;
            }
        }
        if(rN8 > 0){
            unsigned char sto;
            fin.read(reinterpret_cast<char*>(&sto),sizeof(sto));
            for(int k = 0;k < rN8;++k){
                lattice[i][8 * qN8 + k] = 2 * ((sto & (1 << k)) >> k) - 1;
            }
        }
    }
    //test if successfully read
    if(fin.eof()){
        std::cout << "No more snapshots to read!";
        exit(1);
    }
    return make_pair(T_i,lattice);
}