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
    int r_min,r_max,seek_T_i;
    std::cout << "Enter minimum radius: ";
    std::cin >> r_min;
    std::cout << "Enter maximum radius: ";
    std::cin >> r_max;
    std::cout << "Enter T_i: ";
    std::cin >> seek_T_i;
    std::cout << "Calculating two point function for " << std::string(filename) << "..." << std::endl;
    Params p = read_params(fin);
    //Set up progress printing
    int pwidth = std::to_string(p.n_T).length();
    int progress = 0;
    std::cout << "\nProgress: " << std::setw(pwidth) << progress << '/' << p.n_s;
    //stores s[i]*s[i+j] in s12[j][i]
    std::vector<double> s12(r_max-r_min+1,0);
    //stores s[i] in s11[i]
    double s1=0;
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
        if(T_i!=seek_T_i)  continue;   //only take T=1.2
        for(int j=0;j<=r_max-r_min;j++){
            s12[j] += lattice[0][0]*lattice[0][j%p.N];
        }
        s1+=lattice[0][0];
        //Update progress
        #pragma omp critical
        {
            ++progress;
            std::cout << "\rProgress: " << std::setw(pwidth) << progress << '/' << p.n_s;
        }
    }
    //take average
    for(int j=0;j<r_max-r_min+1;j++){
        s12[j] /= p.n_s;
    }
    s1 /= p.n_s;
    //compute G[j]
    std::vector<double> G(r_max-r_min+1,0);
    for(int j=0;j<r_max-r_min+1;j++){
        G[j]+=s12[j]-s1*s1;
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
    std::string outfile_name = "two_point_func_" + std::to_string(p.N) + ".csv";
    std::ofstream fout(outfile_name);
    std::cout << "\nCalculated two point function data for " << std::string(filename);
    std::cout << "\nWriting data to " << outfile_name << "...";
    fout << "r,G(r)\n";
    for(int j=0;j<r_max-r_min+1;j++){
        fout << j+r_min << ',' << G[j] << '\n';
    }
    std::cout << "\nFinished";
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