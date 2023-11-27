#include "reader.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cmath>
#include <chrono>
#include <omp.h>

void analyze(const char* filename){
    std::ifstream fin(filename,std::ios_base::binary);
    if(!fin.is_open()){
        std::cout << "Failed to open " << std::string(filename) << '\n';
        return;
    }
    int seek_T_i;
    std::cout << "Enter T_i: ";
    std::cin >> seek_T_i;
    std::cout << "Calculating histogram data for " << std::string(filename) << "..." << std::endl;
    Params p = read_params(fin);
    std::vector<double> mag,energy;
    //Start timing
    double start = omp_get_wtime();
    #pragma omp parallel for
    for(int s_i = 0;s_i < p.n_s * p.n_T;++s_i){
        int T_i;
        Grid lattice;
        #pragma omp critical
        {
            auto [T,lat] = read_next(p.N, fin);
            T_i=T;
            lattice=lat;
        }
        if(T_i!=seek_T_i)  continue;
        //Magnetisation
        double sum = 0.0;
        for(int i = 0;i < p.N;++i){
            for(int j = 0;j < p.N;++j){
                sum += lattice[i][j];
            }
        }
        sum /= (p.N * p.N);
        mag.push_back(std::abs(sum));
        //Energy
        double sum2 = 0.0;
        for(int i = 0;i < p.N;++i){
            for(int j = 0;j < p.N;++j){
                //Only take left and up to avoid double counting
                sum2 += lattice[i][j] * (lattice[(i-1+p.N)%p.N][j] + lattice[i][(j-1+p.N)%p.N]);
            }
        }
        energy.push_back((-1)*sum2);   //We defined J = 1 initially
    }
    fin.close();
    //End timing and print time taken
    double end = omp_get_wtime();
    const std::chrono::duration<double> time_taken{end - start};
    const auto hrs = std::chrono::duration_cast<std::chrono::hours>(time_taken);
    const auto mins = std::chrono::duration_cast<std::chrono::minutes>(time_taken - hrs);
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(time_taken - hrs - mins);
    const auto millisecs = std::chrono::duration_cast<std::chrono::milliseconds>(time_taken - hrs - mins - secs);
    std::cout << "\nTime taken = " << hrs.count() << "h " << mins.count() << "m " << secs.count() << "s " << millisecs.count() << "ms";
    std::cout << std::endl;
    std::string filename_noext = std::string(filename);
    filename_noext.erase(filename_noext.size()-4);
    std::string outfile_name = "hist_" + filename_noext + ".csv";
    std::ofstream fout(outfile_name);
    std::cout << "Calculated histogram data for " << std::string(filename);
    std::cout << "\nWriting data to " << outfile_name << "...";
    fout << "T,Magnetization,Total Energy\n";
    for(int i = 0;i < p.n_s;i++){
        float cur_T = p.T_min + seek_T_i * (p.T_max - p.T_min) / (p.n_T - 1);
        fout << cur_T << ',' << mag[i] << ',' << energy[i] << '\n';
    }
    std::cout << "\nFinished";
    fout.close();
}

int main(int argc,char* argv[]){
    if(argc == 1){
        std::cout << "No arguments were passed. Please pass a filename to analyse.";
        return 0;
    }else{
        for(int argi = 1;argi < argc;++argi){
            analyze(argv[argi]);
        }
    }

    return 0;
}