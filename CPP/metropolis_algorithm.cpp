#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <omp.h>

//For convenience
using LatticeType = int;
using Grid = std::vector<std::vector<LatticeType>>;
using GridList = std::vector<Grid>;

//Mersenne Twister engine seeded using system time
std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
//Uniform real distribution
std::exponential_distribution<float> urd(1);

void mc_timestep_slow(Grid& lattice,float T,int N){
    for(int i = 0;i < N;++i){
        for(int j = 0;j < N;++j){
            int delE = lattice[i][j] * (lattice[(i-1+N)%N][j] + lattice[(i+1)%N][j] + lattice[i][(j-1+1)%N] + lattice[i][(j+1)%N]);
            lattice[i][j] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
        }
    }
}

void mc_timestep(Grid& lattice,float T,int N){
    //Compute a step for each cell
    //I am doing this in several parts to get rid of modulo operations
    //for imposing periodic boundary conditions. Modulo operations are
    //slow so I expect some improvement...

    //interior
    for(int i = 1;i < N - 1;++i){
        for(int j = 1;j < N - 1;++j){
            int delE = lattice[i][j] * (lattice[i-1][j] + lattice[i+1][j] + lattice[i][j-1] + lattice[i][j+1]);
            lattice[i][j] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
        }
    }
    
    //top edge
    for(int i = 1;i < N - 1;++i){
        int delE = lattice[0][i] * (lattice[N-1][i] + lattice[1][i] + lattice[0][i-1] + lattice[0][i+1]);
        lattice[0][i] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
    }

    //bottom edge
    for(int i = 1;i < N - 1;++i){
        int delE = lattice[N-1][i] * (lattice[N-2][i] + lattice[0][i] + lattice[N-1][i-1] + lattice[N-1][i+1]);
        lattice[N-1][i] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
    }

    //left edge
    for(int i = 1;i < N - 1;++i){
        int delE = lattice[i][0] * (lattice[i][N-1] + lattice[i][1] + lattice[i-1][0] + lattice[i+1][0]);
        lattice[i][0] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
    }

    //right edge
    for(int i = 1;i < N - 1;++i){
        int delE = lattice[i][N-1] * (lattice[i][N-2] + lattice[i][0] + lattice[i-1][N-1] + lattice[i+1][N-1]);
        lattice[i][N-1] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
    }

    //top left corner
    {
        int delE = lattice[0][0] * (lattice[N-1][0] + lattice[1][0] + lattice[0][N-1] + lattice[0][1]);
        lattice[0][0] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
    }

    //top right corner
    {
        int delE = lattice[0][N-1] * (lattice[N-1][N-1] + lattice[1][N-1] + lattice[0][N-2] + lattice[0][0]);
        lattice[0][N-1] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
    }

    //bottom left corner
    {
        int delE = lattice[N-1][0] * (lattice[N-2][0] + lattice[0][0] + lattice[N-1][N-1] + lattice[N-1][1]);
        lattice[N-1][0] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
    }

    //bottom right corner
    {
        int delE = lattice[N-1][N-1] * (lattice[N-2][N-1] + lattice[0][N-1] + lattice[N-1][N-2] + lattice[N-1][0]);
        lattice[N-1][N-1] *= 1 - 2 * std::signbit(T * urd(rng) - delE);
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

    double start = omp_get_wtime();

    //For each value of temperature...
    #pragma omp parallel for
    for(int T_i = 0;T_i < T_count;++T_i){
        //Calculate current temperature
        float cur_T = T_min + T_i * (T_max - T_min) / (T_count - 1);
        
        //Create lattice
        Grid lattice(N,std::vector<LatticeType>(N,1));

        //Allow lattice to equilibriate
        for(int t_i = 0;t_i < t_equilibrium;++t_i){
            mc_timestep_slow(lattice, cur_T, N);
        }

        //Start taking snapshots
        //GridList snaps(snapshot_count);
        for(int snap = 0;snap < snapshot_count;snap++){
            for(int t_i = 0;t_i < snapshot_interval;++t_i){
                mc_timestep_slow(lattice, cur_T, N);
            }
            //snaps[snap] = lattice;
        }
    }

    double end = omp_get_wtime();
    std::cout << "Time taken = " << end - start <<" s";

    return 0;
}