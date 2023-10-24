Compiled using GCC 11.2 with flags ```-fopenmp``` and ```-O2```. To run the simulation, simply run ```metropolis_algorithm.exe```.

~WARNING: Data file generated may be very large.~ I packed bools into chars to get around this. See below.

Data file format (```.dat```):

The data is written in binary format.
First, the simulation parameters are written:
- Lattice size $N$ (int)
- Equilibriation timesteps $t_{eq}$ (int)
- Autocorrelation timesteps $t_a$ (int)
- Number of snapshots $n_s$ (int)
- Min temperature $T_{min}$ (float)
- Max temperature $T_{max}$ (float)
- Number of temperature values to take (including endpoints) $n_T$ (int)

The temperature values are taken at a constant spacing determined by the last 3 parameters:

$$T_i = T_{min} + \frac{i}{n_T-1} \cdot \left(T_{max} - T_{min}\right)$$ where $0 \leq i \leq n_T-1$.

Next, a series of snapshots follow. The total number of snapshots in the file can be calculated as $n_s \cdot n_T$. Note that snapshots are not written into the file in order of temperature. So, a snapshot starts with the value of $i$ in $T_i$ (int). Then, the lattice is written in row-major order as a sequence of (char)s. The bits of the lattice are written in blocks of 8 bits = 1 bytes, that is, a (char). See the source code for details.

Analysis:

C++ does not have a plotting library, so what one could do is write the data to a portable format like CSV and then use Python to do the analysis. This is what I plan to do next.
