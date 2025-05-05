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

# Interpreting the Kernel Density Plot

## Basic Interpretation

 Kernel Density Plot shows the distribution of SNP distances between isolates in a dataset. The x-axis represents the number of SNP differences, while the y-axis shows the density or frequency of those differences.

Each peak represents a cluster of isolates with similar genetic distances. Multiple distinct peaks suggest population structure within the dataset.  Well-separated peaks indicate genetically distinct subpopulations or clades within the species. Broad peaks suggest higher genetic diversity within a cluster, while narrow peaks indicate genetic homogeneity. Deep valleys between peaks suggest clear genetic boundaries between subpopulations.

## Practical Applications

When comparing a specific dataset:

- **Closely related isolates** will form peaks at lower SNP distances (left side of plot)
- **Distantly related isolates** will appear in peaks at higher SNP distances (right side)
- **Outbreak scenarios** often show tight clusters (narrow peaks) at very low SNP distances
- **Endemic circulation** typically shows broader peaks across a range of distances

### Interpretation Pitfalls

- **Sample bias**: Ensure your sample collection is representative
- **Sequencing quality**: Low-quality sequences can artificially inflate SNP counts
- **Recombination**: Areas of high recombination can skew SNP distributions
- **Reference genome choice**: Using a distant reference genome can shift the entire distribution rightward

Only looking at parsimony SNP analysis  (such as with [vSNP](https://github.com/USDA-VS/vSNP3)) focuses on the minimum number of evolutionary changes needed to explain the observed variation, which may underestimate true evolutionary distances in some cases.

# Interpreting the Closest Neighbor SNP Distance Histogram

Histogram reveals several important aspects about your dataset population:

## Basic Interpretation

- **Dominant Clonal Group**: A large peak near zero indicates a predominant cluster of closely related isolates, suggesting recent transmission or a clonal expansion event.

- **Population Structure**: A distinct separation between the main cluster and scattered distant isolates suggests this dataset has a structured population with multiple deep clades.

- **Outliers**: Isolates with very large SNP distances (300-450 SNPs) likely represent genetically distinct strains or subpopulations that have evolved separately.

- **Transmission Patterns**: A tight clustering near zero could indicate recent host-to-host transmission or point source outbreaks.

## Practical Applications

- **Relationship Analysis**: Isolates within the main cluster might represent closely related variants with potential epidemiological connections.

- **Genetic Diversity Assessment**: A wide range of distances demonstrates the extent of genetic variation within the analyzed population.

- **Reference Threshold Establishment**: The median value can serve as a reference point for establishing thresholds that define genetic relatedness in monitoring and surveillance efforts.

- **Evolutionary Group Delineation**: The distinct clusters observed can aid in defining evolutionary groups or lineages for taxonomic and phylogenetic purposes.

Unlike kernel density plots which show smoothed distributions, this histogram clearly quantifies the frequency at each genetic distance interval, making it particularly useful for identifying distinct genetic clusters and inferring evolutionary relationships.
