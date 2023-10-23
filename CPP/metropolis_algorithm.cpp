#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <omp.h>

//For convenience
using LatticeType = int;
using Grid = std::vector<std::vector<LatticeType>>;
using GridList = std::vector<Grid>;

//Mersenne Twister engine seeded using system time
std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
//Uniform real distribution
std::uniform_real_distribution<float> urd(1);

void mc_timestep(LatticeType lattice[][40],float vexp[][9],int T_i,int N){
    for(int i = 0;i < N;++i){
        for(int j = 0;j < N;++j){
            LatticeType u = lattice[i][j];
            LatticeType s = lattice[(i-1+N)%N][j] + lattice[(i+1)%N][j] + lattice[i][(i+1)%N] + lattice[i][(j-1+N)%N];
            float accept = urd(rng) - vexp[T_i][u * s + 4];
            lattice[i][j] = std::copysign(u,accept);
        }
    }
}

int main(){
    //This is a C++ program that implements a naive implementation of the Metropolis algorithm.
    //I'm making this to test how much faster C++ is than Python
    //Conclusion: Much faster. I've also used OpenMP to parallelize this.
    //Compiled with flags -fopenmp and -O2.

    //Input data (part 1)
    int N, t_equilibrium, snapshot_interval, snapshot_count;
    std::cout<<"Input size of lattice (N): ";
    std::cin>>N;
    std::cout<<"Input equilibriation time: ";
    std::cin>>t_equilibrium;
    std::cout<<"Input autocorrelation time: ";
    std::cin>>snapshot_interval;
    std::cout<<"Input number of snapshots to take: ";
    std::cin>>snapshot_count;

    //Input data (part 2)
    float T_min, T_max;
    int T_count;
    std::cout<<"Input minimum value of temperature: ";
    std::cin>>T_min;
    std::cout<<"Input maximum value of temperature: ";
    std::cin>>T_max;
    std::cout<<"Input number of values to take (including endpoints): ";
    std::cin>>T_count;

    int pwidth = std::to_string(T_count).length();
    int progress = 0;
    std::cout << "\nProgress: " << std::setw(pwidth) << progress << '/' << T_count;
    double start = omp_get_wtime();

    //Create lattice
    LatticeType lattices[T_count][40][40];
    for(int k = 0;k < T_count;++k){
        for(int i = 0;i < N;++i){
            for(int j = 0;j < N;++j){
                lattices[k][i][j] = 1;
            }
        }
    }

    //Precalculate exponentials...because exp is VERY expensive!
    //This achieved a speedup from ~20s to ~5s (!!) for N = 16
    //and t_eq = 10^5 and T = [0.1,4] in steps of 0.1
    float vexp[T_count][9];
    for(int T_i = 0;T_i < T_count;++T_i){
        for(int j = -4;j <= 4;++j){
            float cur_T = T_min + T_i * (T_max - T_min) / (T_count - 1);
            vexp[T_i][j + 4] = std::exp(j / cur_T);
        }
    }

    //For each value of temperature...
    #pragma omp parallel for
    for(int T_i = 0;T_i < T_count;++T_i){
        //Calculate current temperature
        float cur_T = T_min + T_i * (T_max - T_min) / (T_count - 1);
        
        //Allow lattice to equilibriate
        for(int t_i = 0;t_i < t_equilibrium;++t_i){
            mc_timestep(lattices[T_i], vexp, T_i, N);
        }

        //Start taking snapshots
        //GridList snaps(snapshot_count);
        for(int snap = 0;snap < snapshot_count;snap++){
            for(int t_i = 0;t_i < snapshot_interval;++t_i){
                mc_timestep(lattices[T_i], vexp, T_i, N);
            }
            //snaps[snap] = lattice;
        }

        #pragma omp critical
        {
            ++progress;
            std::cout << "\rProgress: " << std::setw(pwidth) << progress << '/' << T_count;
        }
    }

    double end = omp_get_wtime();
    std::cout << "\nTime taken = " << end - start <<" s";

    return 0;
}