#include <fstream>
#include <vector>
#include <utility>

struct Params{
    int N;  //lattice size
    int t_eq;   //equilibriation timesteps
    int t_a;    //autocorrelation timesteps
    int n_s;    //number of snapshots for fixed temp
    float T_min;    //min temp
    float T_max;    //max temp
    int n_T;    //number of temp values taken (including endpoints)
};

Params read_params(std::ifstream& fin);
std::pair<int,std::vector<std::vector<int>>> read_next(const int N, std::ifstream& fin);