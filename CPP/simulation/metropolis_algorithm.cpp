#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include "common_defs.h"
#include "writer.h"

//Uniform real distribution
std::uniform_real_distribution<float> urd(0.0,1.0);

//Monte - Carlo timestep
void mc_timestep(Grid& lattice,float vexp[][9],int T_i,int N, std::mt19937& rng){
    //Execute one step for cell
    for(int i = 0;i < N;++i){
        for(int j = 0;j < N;++j){
            LatticeType u = -lattice[i][j];
            LatticeType s = lattice[(i-1+N)%N][j] + lattice[(i+1)%N][j] + lattice[i][(j+1)%N] + lattice[i][(j-1+N)%N];
            float accept = urd(rng) - vexp[T_i][u * s + 4];
            lattice[i][j] *= 1 - 2 * std::signbit(accept);
        }
    }
}

int main(){
    //This is a C++ program that implements a naive implementation of the Metropolis algorithm.
    //I'm making this to test how much faster C++ is than Python.
    //Conclusion: Much faster. I've also used OpenMP to parallelize this.

    Params sim_params;
    //Input data (part 1)
    std::cout << "Input size of lattice (N): ";
    std::cin >> sim_params.N;
    std::cout << "Input equilibriation time: ";
    std::cin >> sim_params.t_eq;
    std::cout << "Input autocorrelation time: ";
    std::cin >> sim_params.t_a;
    std::cout << "Input number of snapshots to take: ";
    std::cin >> sim_params.n_s;

    //Input data (part 2)
    std::cout << "Input minimum value of temperature: ";
    std::cin >> sim_params.T_min;
    std::cout << "Input maximum value of temperature: ";
    std::cin >> sim_params.T_max;
    std::cout << "Input number of values to take (including endpoints): ";
    std::cin >> sim_params.n_T;

    //Input data (part 3)
    char create_file;
    std::cout << "Create data file? (y/n): ";
    std::cin >> create_file;

    std::ofstream fout;
    if(create_file == 'y'){
        //Create output file
        std::string filename_noext;
        std::cout << "Name of file to create (without .dat extension): ";
        std::cin >> filename_noext;
        std::string filename = filename_noext + ".dat";
        fout.open(filename, std::ios_base::binary);
        write_params(sim_params, fout);
    }

    //Set up progress printing
    int pwidth = std::to_string(sim_params.n_T).length();
    int progress = 0;
    std::cout << "\nProgress: " << std::setw(pwidth) << progress << '/' << sim_params.n_T;

    //Create lattices
    GridList lattices(sim_params.n_T,Grid(sim_params.N,Row(sim_params.N,1)));

    //Precalculate exponentials...because exp is VERY expensive!
    //This achieved a speedup from ~20s to ~5s (!!) for N = 16
    //and t_eq = 10^5 and T = [0.1,4] in steps of 0.1
    float vexp[sim_params.n_T][9];
    for(int T_i = 0;T_i < sim_params.n_T;++T_i){
        float cur_T = sim_params.T_min + T_i * (sim_params.T_max - sim_params.T_min) / (sim_params.n_T - 1);
        
        for(int j = -4;j < 0;++j){
            vexp[T_i][j + 4] = std::exp(j / cur_T);
        }
        for(int j = 0;j <= 4;++j){
            vexp[T_i][j + 4] = 1.0;
        }
    }

    //Start timing
    double start = omp_get_wtime();

    //For each value of temperature...
    #pragma omp parallel for
    for(int T_i = 0;T_i < sim_params.n_T;++T_i){
        //Mersenne Twister engine seeded using system time
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        
        //Allow lattice to equilibriate
        for(int t_i = 0;t_i < sim_params.t_eq;++t_i){
            mc_timestep(lattices[T_i], vexp, T_i, sim_params.N, rng);
        }

        //Start taking snapshots
        for(int snap = 0;snap < sim_params.n_s;snap++){
            for(int t_i = 0;t_i < sim_params.t_a;++t_i){
                mc_timestep(lattices[T_i], vexp, T_i, sim_params.N, rng);
            }
            if(create_file == 'y'){
                //Write snapshot to file
                #pragma omp critical
                write_next(T_i, sim_params.N, lattices[T_i], fout);
            }
        }

        //Update progress
        #pragma omp critical
        {
            ++progress;
            std::cout << "\rProgress: " << std::setw(pwidth) << progress << '/' << sim_params.n_T;
        }
    }

    if(create_file == 'y'){
        //Close output file
        fout.close();
    }

    //End timing and print time taken
    double end = omp_get_wtime();
    const std::chrono::duration<double> time_taken{end - start};
    const auto hrs = std::chrono::duration_cast<std::chrono::hours>(time_taken);
    const auto mins = std::chrono::duration_cast<std::chrono::minutes>(time_taken - hrs);
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(time_taken - hrs - mins);
    const auto millisecs = std::chrono::duration_cast<std::chrono::milliseconds>(time_taken - hrs - mins - secs);
    std::cout << "\nTime taken = " << hrs.count() << "h " << mins.count() << "m " << secs.count() << "s " << millisecs.count() << "ms";

    return 0;
}