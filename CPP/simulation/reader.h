#include <fstream>
#include <vector>
#include <utility>
#include "common_defs.h"

Params read_params(std::ifstream& fin);
std::pair<int,Grid> read_next(const int N, std::ifstream& fin);