#ifndef COMMONDEFS_H
#define COMMONDEFS_H

#include <vector>

struct Params{
    int N;  //lattice size
    int t_eq;   //equilibriation timesteps
    int t_a;    //autocorrelation timesteps
    int n_s;    //number of snapshots for fixed temp
    float T_min;    //min temp
    float T_max;    //max temp
    int n_T;    //number of temp values taken (including endpoints)
};

//For convenience
using LatticeType = int;
using Row = std::vector<LatticeType>;
using Grid = std::vector<Row>;
using GridList = std::vector<Grid>;

#endif /*COMMONDEFS_H*/