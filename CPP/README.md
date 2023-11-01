Compiled using GCC 11.2 with flags ```-fopenmp``` and ```-O3```. To run the simulation, simply run ```metropolis_algorithm.exe```.

WARNING: Data file generated may be very large.

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

I have written a data reader interface so that one doesn't have to re-implement the data decoding process. There are two functions ```read_start``` and ```read_next```. First you open the file you want to read and pass it to ```read_start``` which returns a ```Params``` struct containing  the simulation parameters. Next, you can loop $n_s \cdot n_T$ times and each time you call ```read_next``` which returns a pair of $i$ and the corresponding decoded lattice. Now you can process the data as you like.

A simple data checker program is provided that takes a list of data file names as command-line arguments and prints a summary of the contents of each of the data files.

I have written programs for finding:

- Spontaneous magnetisation per unit spin
- Susceptibility per unit spin
- Energy per unit spin
- Specific heat per unit spin

To use the programs, just pass the names of the data files to be analysed as command-line arguments. C++ does not have a plotting library, so what one could do is write the data to a portable format like CSV and then use Python to do the analysis. For now, the programs just output the values to the console.
