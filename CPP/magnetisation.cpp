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
    std::cout << "Calculating magnetization for " << std::string(filename) << "..." << std::endl;
    Params p = read_params(fin);
    std::map<int,double> magnetisation;
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
        double sum = 0.0;
        for(int i = 0;i < p.N;++i){
            for(int j = 0;j < p.N;++j){
                sum += lattice[i][j];
            }
        }
        sum /= (p.N * p.N);
        magnetisation[T_i] += std::abs(sum);
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
    std::string outfile_name = "magnetisation_" + filename_noext + ".csv";
    std::ofstream fout(outfile_name);
    std::cout << "Calculated magnetization data for " << std::string(filename);
    std::cout << "\nWriting data to " << outfile_name << "...";
    fout << "T,Magnetization\n";
    for(auto [T_i,m]:magnetisation){
        float cur_T = p.T_min + T_i * (p.T_max - p.T_min) / (p.n_T - 1);
        fout << cur_T << ',' << m / p.n_s << '\n';
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