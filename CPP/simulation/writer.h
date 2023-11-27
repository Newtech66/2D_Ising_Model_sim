#include <fstream>
#include <vector>
#include <utility>
#include "common_defs.h"

void write_params(const Params& sim_params, std::ofstream& fout);
void write_next(const int T_i, const int N, const Grid& lattice, std::ofstream& fout);