# Kernel Density Plot Tool

This is a Python implementation of the Kernel Plot Tool that generates SNP density and closest neighbor plots from aligned FASTA files.

## Installation

### Setting Up Conda Environment

1. Ensure [Conda](https://docs.conda.io/en/latest/miniconda.html) is installed on your system.

2. Create a new environment:

```bash
conda create -n kernel_density_plots kernel_density_plots
```

3. Activate the environment:

```bash
conda activate kernel_density_plots
```

## Usage

The tool provides a command-line interface:

```bash
kernel_plot.py [options] -f input.fasta
```

## Test

```bash
cd ${HOME}
mkdir kernel_test
cd kernel_test
```

Copy test file from conda install
```bash
cp -v $CONDA_PREFIX/share/kernel_density_plots/test/*fasta .
```

Run test file
```bash
kernel_plot.py --lineage="Lineage L2" --bin-size=20 -f La2_test.fasta
```
## Options

See `kernel_plot.py -h` for full list of options

## Output Files

The tool generates the following output files in the specified output directory:

1. `*_distances.tab` - Tab-delimited distance matrix
2. `*_no_root_YYYY-MM-DD.tab` - Distance matrix with root sequence removed (if present)
3. `*_lowertriangle.txt` - Lower triangle values from the distance matrix
4. `*_density_plot.pdf` - SNP density plot
5. `*_closest_neighbor.pdf` - Closest neighbor histogram
6. `*_combined_figure.pdf` - Combined figure with both plots
7. `*_density_bins.csv` - CSV file with density histogram bin counts
8. `*_neighbor_bins.csv` - CSV file with neighbor histogram bin counts
9. `*_density_stats.txt` - Text file with density statistics
10. `*_neighbor_stats.txt` - Text file with neighbor statistics
