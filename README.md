# Kernel Density Plot Tool

This is a Python implementation of the Kernel Plot Tool that generates SNP density and closest neighbor plots from FASTA files.

## Installation

### Setting Up Conda Environment

1. Make sure you have [Conda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.

2. Create a new conda environment using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
```

3. Activate the environment:

```bash
conda activate kernel_density_plots
```

## Usage

The tool provides a command-line interface:

```bash
python kernel_plot.py [options] input.fasta
```

### Options:

```
--lineage=NAME       Specify lineage (optional)
--outputdir=DIR      Output directory (default: creates folder next to input file)
--density-xlim=N     Density plot X-axis limit (default: 1400)
--neighbor-xlim=N    Closest neighbor X-axis limit (default: 600)
--bin-size=N         Histogram bin size (default: 25)
--width=N            Plot width in inches (default: 7)
--height=N           Plot height in inches (default: 5)
--title-size=N       Title text size (default: 18)
--axis-title-size=N  Axis title text size (default: 14)
--axis-text-size=N   Axis text size (default: 12)
--annotation-size=N  Annotation text size (default: 6)
--no-annotations     Hide text annotations on plots
--theme=NAME         Plot theme (light, minimal, classic, gray, dark, bw) (default: light)
```

### Example:

```bash
python kernel_plot.py --lineage="Lineage L2" --bin-size=20 La2_test.fasta
```

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