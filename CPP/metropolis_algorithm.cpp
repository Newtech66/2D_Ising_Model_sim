#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <omp.h>

//Mersenne Twister engine seeded using system time
std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
//Uniform real distribution
std::uniform_real_distribution<float> urd(0.0,1.0);

void mc_step(std::vector<std::vector<int>>& lattice, float T,int N, int i, int j){
    //Compute delE
    int delE = lattice[i][j] * (lattice[(i-1+N)%N][j] + lattice[(i+1)%N][j] + lattice[i][(j-1+N)%N] + lattice[i][(j+1)%N]);

    //Check if can accept
    if(urd(rng) < std::exp(-delE/T)){
        //If accepted then flip
        lattice[i][j] = -lattice[i][j];
    }
}

void mc_timestep(std::vector<std::vector<int>>& lattice, float T,int N){
    //Compute a step for each cell
    for(int i = 0;i < N;i++){
        for(int j = 0;j < N;j++){
            mc_step(lattice, T, N, i, j);
        }
    }
}

int main(){
    //This is a C++ program that implements a naive implementation of the Metropolis algorithm.
    //I'm making this to test how much faster C++ is.

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

    // auto start = std::chrono::steady_clock::now();

    //For each value of temperature...
    #pragma omp parallel for
    for(int T_i = 0;T_i < T_count;T_i++){
        //Calculate current temperature
        float cur_T = T_min + T_i * (T_max - T_min) / (T_count - 1);
        
        //Create lattice
        std::vector<std::vector<int>> lattice(N,std::vector<int>(N,1));

        //Allow lattice to equilibriate
        for(int t_i = 0;t_i < t_equilibrium;t_i++){
            mc_timestep(lattice, cur_T, N);
        }

        //Start taking snapshots
        for(int snap = 0;snap < snapshot_count;snap++){
            for(int t_i = 0;t_i < snapshot_interval;t_i++){
                mc_timestep(lattice, cur_T, N);
            }
        }
    }

    // auto end = std::chrono::steady_clock::now();
    // auto diff = end - start;

    // std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms";

    return 0;
}