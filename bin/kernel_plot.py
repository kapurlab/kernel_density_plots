#!/usr/bin/env python3

# Kernel Plot Tool
# A command-line tool for SNP density and closest neighbor analysis
# Accepts a FASTA file as input and generates kernel density and closest neighbor plots

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from Bio import SeqIO
from scipy import signal
from scipy import stats
import seaborn as sns
from pathlib import Path
import textwrap
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Kernel Plot Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
        Example:
          python kernel_plot.py --lineage=Lineage_4 --bin-size=20 input.fasta
        ''')
    )
    parser.add_argument('input_fasta', type=str, help='Input FASTA file')
    parser.add_argument('--lineage', type=str, default='Unknown',
                        choices=['Lineage_1', 'Lineage_2', 'Lineage_3', 'Lineage_4', 'Bovis', 'Caprae', 'Orygis', 'Unknown'],
                        help='Specify lineage (optional)')
    parser.add_argument('--outputdir', type=str, help='Output directory (default: creates folder next to input file)')
    parser.add_argument('--density-xlim', type=int, default=1400, help='Density plot X-axis limit (default: 1400)')
    parser.add_argument('--neighbor-xlim', type=int, default=600, help='Closest neighbor X-axis limit (default: 600)')
    parser.add_argument('--bin-size', type=int, default=25, help='Histogram bin size (default: 25)')
    parser.add_argument('--width', type=float, default=7, help='Plot width in inches (default: 7)')
    parser.add_argument('--height', type=float, default=5, help='Plot height in inches (default: 5)')
    parser.add_argument('--title-size', type=int, default=18, help='Title text size (default: 18)')
    parser.add_argument('--axis-title-size', type=int, default=14, help='Axis title text size (default: 14)')
    parser.add_argument('--axis-text-size', type=int, default=12, help='Axis text size (default: 12)')
    parser.add_argument('--annotation-size', type=int, default=6, help='Annotation text size (default: 6)')
    parser.add_argument('--no-annotations', action='store_true', help='Hide text annotations on plots')
    parser.add_argument('--theme', type=str, default='light',
                      choices=['light', 'minimal', 'classic', 'gray', 'dark', 'bw'],
                      help='Plot theme (default: light)')
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)
        
    return parser.parse_args()

def detect_lineage(filename):
    """Detect lineage from filename"""
    filename = filename.lower()
    
    lineage_patterns = {
        "Lineage_1": ["lineage_1", "lineage1", "l1", "indo-oceanic", "eai"],
        "Lineage_2": ["lineage_2", "lineage2", "l2", "beijing", "east-asian"],
        "Lineage_3": ["lineage_3", "lineage3", "l3", "delhi", "central-asian"],
        "Lineage_4": ["lineage_4", "lineage4", "l4", "euro-american"],
        "Bovis": ["bovis", "bov"],
        "Caprae": ["caprae", "cap"],
        "Orygis": ["orygis", "ory"]
    }
    
    for lineage_name, patterns in lineage_patterns.items():
        for pattern in patterns:
            if pattern in filename:
                return lineage_name
    
    # Default to unknown
    return "Unknown"

def calculate_snp_distances(fasta_file):
    """Calculate SNP distances from FASTA file"""
    print("Reading FASTA file...")
    
    # Read sequences from FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    seq_count = len(sequences)
    
    print(f"Processing {seq_count} sequences...")
    
    # Initialize the distance matrix
    dist_matrix = np.zeros((seq_count, seq_count), dtype=int)
    seq_names = [seq.id for seq in sequences]
    
    # Calculate pairwise distances
    for i in range(seq_count-1):
        for j in range(i+1, seq_count):
            # Get the two sequences
            seq1 = str(sequences[i].seq).lower()
            seq2 = str(sequences[j].seq).lower()
            
            # Check if sequences have the same length
            if len(seq1) != len(seq2):
                print(f"Warning: Sequences {i} and {j} have different lengths. Using minimum length.")
                min_length = min(len(seq1), len(seq2))
                seq1 = seq1[:min_length]
                seq2 = seq2[:min_length]
            
            # Count positions where sequences differ, excluding positions with gaps or N
            diff_count = 0
            valid_positions = 0
            
            for pos in range(len(seq1)):
                base1 = seq1[pos]
                base2 = seq2[pos]
                
                # Skip if either position has a gap ('-') or ambiguous base (N)
                if base1 == "-" or base2 == "-" or base1 == "n" or base2 == "n":
                    continue
                
                valid_positions += 1
                if base1 != base2:
                    diff_count += 1
            
            # Store distance in matrix (symmetric)
            dist_matrix[i, j] = diff_count
            dist_matrix[j, i] = diff_count
    
    # Create DataFrame with sequence names
    dist_df = pd.DataFrame(dist_matrix, index=seq_names, columns=seq_names)
    
    print("SNP distance calculation complete.")
    return dist_df

def calculate_bin_counts(data, breaks, include_stats=True):
    """Calculate and format histogram bin counts"""
    # Create the histogram data
    hist, bin_edges = np.histogram(data, bins=breaks)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Create a DataFrame with bin information
    bin_df = pd.DataFrame({
        'bin_start': bin_edges[:-1],
        'bin_end': bin_edges[1:],
        'bin_mid': bin_mids,
        'count': hist,
        'percentage': np.round(hist / np.sum(hist) * 100, 1)
    })
    
    # Format the output
    output_lines = []
    
    # Add statistics if requested
    if include_stats:
        output_lines.extend([
            f"Min: {min(data)}",
            f"Max: {max(data)}",
            f"Median: {np.median(data)}",
            f"Mean: {round(np.mean(data), 1)}",
            f"Standard Deviation: {round(np.std(data), 1)}",
            f"Total count: {len(data)}",
            ""
        ])
    
    # Add bin header
    output_lines.extend([
        "Bin Counts:",
        f"{'Bin Range':<20}{'Count':<10}{'Percentage':<12}",
        "-" * 42
    ])
    
    # Add bin data
    for _, row in bin_df.iterrows():
        bin_range = f"[{row['bin_start']}, {row['bin_end']})"
        output_lines.append(
            f"{bin_range:<20}{row['count']:<10}{row['percentage']}%"
        )
    
    # Return the formatted output and the bin data frame
    return {
        'text': '\n'.join(output_lines),
        'data': bin_df
    }

def process_fasta(input_fasta, lineage, output_folder=None):
    """Process a FASTA file and calculate SNP distances"""
    # Expand the path to handle tilde (~) in file paths
    input_fasta = os.path.expanduser(input_fasta)
    
    # Create output file names based on input
    input_file_name = os.path.splitext(os.path.basename(input_fasta))[0]
    output_file = f"{input_file_name}_distances.tab"
    
    # Create output folder in same directory as input file if not provided
    if output_folder is None:
        input_dir = os.path.dirname(input_fasta)
        output_folder = os.path.join(input_dir, f"KDP_{input_file_name}")
    else:
        # Also expand output folder path if provided
        output_folder = os.path.expanduser(output_folder)
    
    os.makedirs(output_folder, exist_ok=True)
    
    # Define the output file path
    full_output_file = os.path.join(output_folder, output_file)
    
    # Calculate SNP distances
    print(f"Calculating SNP distances...")
    print(f"Input file: {input_fasta}")
    
    # Call our SNP distance calculator
    dist_matrix = calculate_snp_distances(input_fasta)
    
    # Save the distance matrix as a tab-delimited file
    print(f"Saving distance matrix to {full_output_file}")
    dist_matrix.to_csv(full_output_file, sep='\t')
    
    # Process the output file
    print("Processing distance matrix...")
    distances = dist_matrix.copy()
    
    # Remove root sequence if present
    if "root" in distances.index:
        distances_noroot = distances.drop("root", axis=0).drop("root", axis=1)
    else:
        distances_noroot = distances
    
    # Save the modified table
    from datetime import date
    output_noroot_file = os.path.join(output_folder, f"{input_file_name}_no_root_{date.today().strftime('%Y-%m-%d')}.tab")
    distances_noroot.to_csv(output_noroot_file, sep='\t')
    
    # Copy/move files to output folder
    from shutil import copy
    copy(input_fasta, os.path.join(output_folder, os.path.basename(input_fasta)))
    
    # Extract lower triangle values
    data = distances_noroot.values
    lower_triangle_values = data[np.tril_indices(len(data), k=-1)]
    
    # Save lower triangle values
    output_filename = os.path.join(output_folder, f"{input_file_name}_lowertriangle.txt")
    pd.DataFrame(lower_triangle_values).to_csv(output_filename, sep='\t', index=True, header=False)
    
    print(f"Lower triangle values saved to {output_filename}")
    
    return {
        'lower_triangle_file': output_filename,
        'output_folder': output_folder,
        'input_file_name': input_file_name,
        'data': distances_noroot,
        'lower_triangle_values': lower_triangle_values,
        'lineage': lineage
    }

def create_density_plot(result, xlim_max=1400, show_annotations=True,
                       title_size=20, axis_title_size=16, axis_text_size=12,
                       annotation_size=6, plot_theme='light', output_path=None):
    """Create the kernel density plot with styling based on lineage"""
    # Read the lower triangle values
    my_data_num = result['lower_triangle_values']
    
    # Set xlims
    xlims = (0, xlim_max)
    
    # Set colors based on lineage
    lineage_colors = {
        "Lineage_1": ["#F0E1B3", "#E7BD42", "black"],  # light pastel yellow, pastel yellow, black
        "Lineage_2": ["#C1F0C1", "#6BCE7D", "black"],  # light mint green, mint green, black
        "Lineage_3": ["#B3D9FF", "#4996D5", "black"],  # light sky blue, sky blue, black
        "Lineage_4": ["#FEBAB3", "#E3687C", "black"],  # light coral pink, medium coral, black
        "Caprae": ["#ADD8E6", "#1E90FF", "black"],     # light blue, dodger blue, black
        "Orygis": ["#90EE90", "#32CD32", "black"],     # light green, lime green, black
        "Bovis": ["#E6E6FA", "#9370DB", "red"]         # lavender, medium purple, red
    }
    
    # Get colors for current lineage, default to first one if not found
    if result['lineage'] not in lineage_colors:
        colors = list(lineage_colors.values())[0]
    else:
        colors = lineage_colors[result['lineage']]
    
    polygon_fill_color = colors[0]
    iqr_highlight_color = colors[1]
    median_line_color = colors[2]
    
    # For dark theme, ensure text is visible
    text_color = "white" if plot_theme == "dark" else "black"
    
    # Apply theme
    set_plot_theme(plot_theme)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Calculate KDE
    kde = stats.gaussian_kde(my_data_num)
    x = np.linspace(0, xlim_max, 1000)
    y = kde(x)
    
    # Calculate IQR
    q1, q3 = np.percentile(my_data_num, [25, 75])
    iqr_x = x[(x >= q1) & (x <= q3)]
    iqr_y = kde(iqr_x)
    
    # Find peaks in the density distribution
    peaks, _ = signal.find_peaks(y, height=max(y) * 0.1, distance=5)
    peak_values = x[peaks]
    
    # Plot the density
    ax.fill_between(x, y, alpha=0.7, color=polygon_fill_color, edgecolor="#2c2c2c", linewidth=0.7)
    
    # Add IQR highlight
    ax.fill_between(iqr_x, iqr_y, alpha=0.3, color=iqr_highlight_color)
    
    # Add median line
    median_val = np.median(my_data_num)
    ax.axvline(median_val, color="red", linestyle="solid", linewidth=1.5)
    
    # Add peaks if they exist and annotations are enabled
    if len(peak_values) > 0 and show_annotations:
        for peak in peak_values:
            ax.axvline(peak, color='red', linestyle='dashed')
            ax.text(peak, max(y) * 0.9, f"{int(round(peak, 0))}", 
                    va='center', ha='center', color='red', fontsize=annotation_size*2)
    
    # Add text annotation with statistics if enabled
    if show_annotations:
        stats_text = f"Min: {min(my_data_num)}\nMax: {max(my_data_num)}\nMedian: {round(median_val, 1)}"
        ax.text(xlim_max * 0.75, max(y) * 0.8, stats_text,
                ha='left', va='top', fontsize=annotation_size*2, color=text_color, style='italic')
    
    # Set limits and labels
    ax.set_xlim(xlims)
    ax.set_xlabel("SNP distance", fontsize=axis_title_size, fontweight='bold')
    ax.set_ylabel("Density", fontsize=axis_title_size, fontweight='bold')
    
    # Set title
    number_of_sequences = len(result['data'])
    ax.set_title(f"{result['lineage']} SNP density plot (n = {number_of_sequences})", 
                fontsize=title_size, fontweight='bold')
    
    # Style the plot
    apply_panel_styling(ax, axis_text_size)
    
    # Save the plot if output path is provided
    if output_path:
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
    
    return fig, ax

def create_closest_neighbor_plot(result, bin_size=25, xlim_max=600, show_annotations=True,
                                title_size=20, axis_title_size=16, axis_text_size=12,
                                annotation_size=6, plot_theme='light', output_path=None):
    """Create a closest neighbor histogram"""
    # Get the data matrix and convert to float to handle NaN values
    data = result['data'].values.astype(float)
    
    # Replace diagonal elements with NaN
    np.fill_diagonal(data, np.nan)
    
    # Find the lowest pairwise distance for each row
    lowest_distances = np.nanmin(data, axis=1)
    
    # Define breaks
    breaks = np.arange(0, xlim_max + bin_size, bin_size)
    
    # Calculate median
    median_val = np.median(lowest_distances)
    
    # Set colors based on lineage (same as density plot)
    lineage_colors = {
        "Lineage_1": "#F0E1B3",  # light pastel yellow
        "Lineage_2": "#C1F0C1",  # light mint green
        "Lineage_3": "#B3D9FF",  # light sky blue
        "Lineage_4": "#FEBAB3",  # light coral pink
        "Caprae": "#ADD8E6",     # light blue
        "Orygis": "#90EE90",     # light green
        "Bovis": "#E6E6FA"       # lavender
    }
    
    # Get color for current lineage, default to first one if not found
    if result['lineage'] not in lineage_colors:
        fill_color = list(lineage_colors.values())[0]
    else:
        fill_color = lineage_colors[result['lineage']]
    
    # For dark theme, ensure text is visible
    text_color = "white" if plot_theme == "dark" else "black"
    
    # Apply theme
    set_plot_theme(plot_theme)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Create histogram
    n, bins, patches = ax.hist(lowest_distances, bins=breaks, color=fill_color, 
                              edgecolor='black', alpha=0.7)
    
    # Add median line
    ax.axvline(median_val, color="red", linestyle="solid", linewidth=1.5)
    
    # Add text annotation with median value if enabled
    if show_annotations:
        ax.text(xlim_max * 0.75, max(n), f"Median = {round(median_val, 1)}",
                ha='left', va='top', fontsize=annotation_size*2, color=text_color, style='italic')
    
    # Set limits and labels
    ax.set_xlim(0, xlim_max)
    ax.set_xlabel("SNP distance to closest neighbor", fontsize=axis_title_size, fontweight='bold')
    ax.set_ylabel("Frequency", fontsize=axis_title_size, fontweight='bold')
    
    # Set title
    ax.set_title(f"{result['lineage']} Closest Neighbor (n = {len(result['data'])})", 
                fontsize=title_size, fontweight='bold')
    
    # Style the plot
    apply_panel_styling(ax, axis_text_size)
    
    # Save the plot if output path is provided
    if output_path:
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
    
    return fig, ax

def create_combined_figure(result, density_xlim, neighbor_xlim, bin_size, show_annotations=True,
                         title_size=20, axis_title_size=16, axis_text_size=12,
                         annotation_size=6, plot_theme='light', plot_width=7, plot_height=5):
    """Create a combined figure with density plot and closest neighbor plot side by side"""
    # Apply theme
    set_plot_theme(plot_theme)
    
    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(plot_width * 2, plot_height))
    
    # Determine text color based on theme
    text_color = "white" if plot_theme == "dark" else "black"
    
    # DENSITY PLOT in left subplot
    # Get data
    my_data_num = result['lower_triangle_values']
    
    # Calculate KDE
    kde = stats.gaussian_kde(my_data_num)
    x = np.linspace(0, density_xlim, 1000)
    y = kde(x)
    
    # Calculate IQR
    q1, q3 = np.percentile(my_data_num, [25, 75])
    iqr_x = x[(x >= q1) & (x <= q3)]
    iqr_y = kde(iqr_x)
    
    # Get colors for lineage
    lineage_colors = {
        "Lineage_1": ["#F0E1B3", "#E7BD42", "black"],
        "Lineage_2": ["#C1F0C1", "#6BCE7D", "black"],
        "Lineage_3": ["#B3D9FF", "#4996D5", "black"],
        "Lineage_4": ["#FEBAB3", "#E3687C", "black"],
        "Caprae": ["#ADD8E6", "#1E90FF", "black"],
        "Orygis": ["#90EE90", "#32CD32", "black"],
        "Bovis": ["#E6E6FA", "#9370DB", "red"]
    }
    
    if result['lineage'] not in lineage_colors:
        colors = list(lineage_colors.values())[0]
    else:
        colors = lineage_colors[result['lineage']]
        
    # Plot density
    ax1.fill_between(x, y, alpha=0.7, color=colors[0], edgecolor="#2c2c2c", linewidth=0.7)
    ax1.fill_between(iqr_x, iqr_y, alpha=0.3, color=colors[1])
    
    # Add median line
    median_val = np.median(my_data_num)
    ax1.axvline(median_val, color="red", linestyle="solid", linewidth=1.5)
    
    # Find peaks in the density distribution
    peaks, _ = signal.find_peaks(y, height=max(y) * 0.1, distance=5)
    peak_values = x[peaks]
    
    # Add peaks if they exist and annotations are enabled
    if len(peak_values) > 0 and show_annotations:
        for peak in peak_values:
            ax1.axvline(peak, color='red', linestyle='dashed')
            ax1.text(peak, max(y) * 0.9, f"{int(round(peak, 0))}", 
                    va='center', ha='center', color='red', fontsize=annotation_size*2)
    
    # Add annotation
    if show_annotations:
        stats_text = f"Min: {min(my_data_num)}\nMax: {max(my_data_num)}\nMedian: {round(median_val, 1)}"
        ax1.text(density_xlim * 0.75, max(y) * 0.8, stats_text,
                ha='left', va='top', fontsize=annotation_size*2, color=text_color, style='italic')
    
    # Set limits and labels
    ax1.set_xlim(0, density_xlim)
    ax1.set_xlabel("SNP distance", fontsize=axis_title_size, fontweight='bold')
    ax1.set_ylabel("Density", fontsize=axis_title_size, fontweight='bold')
    ax1.set_title(f"{result['lineage']} SNP density plot (n = {len(result['data'])})", 
                fontsize=title_size*0.8, fontweight='bold')
    
    # Style the plot
    apply_panel_styling(ax1, axis_text_size)
    
    # CLOSEST NEIGHBOR PLOT in right subplot
    # Get data for closest neighbor plot
    data_float = result['data'].values.astype(float)
    np.fill_diagonal(data_float, np.nan)
    lowest_distances = np.nanmin(data_float, axis=1)
    
    # Set colors for histogram
    if result['lineage'] not in lineage_colors:
        fill_color = list(lineage_colors.values())[0][0]
    else:
        fill_color = lineage_colors[result['lineage']][0]
    
    # Create breaks for histogram
    breaks = np.arange(0, neighbor_xlim + bin_size, bin_size)
    
    # Plot histogram
    n, bins, patches = ax2.hist(lowest_distances, bins=breaks, color=fill_color, 
                              edgecolor='black', alpha=0.7)
    
    # Add median line
    neighbor_median = np.median(lowest_distances)
    ax2.axvline(neighbor_median, color="red", linestyle="solid", linewidth=1.5)
    
    # Add annotation
    if show_annotations:
        ax2.text(neighbor_xlim * 0.75, max(n), f"Median = {round(neighbor_median, 1)}",
                ha='left', va='top', fontsize=annotation_size*2, color=text_color, style='italic')
    
    # Set limits and labels
    ax2.set_xlim(0, neighbor_xlim)
    ax2.set_xlabel("SNP distance to closest neighbor", fontsize=axis_title_size, fontweight='bold')
    ax2.set_ylabel("Frequency", fontsize=axis_title_size, fontweight='bold')
    ax2.set_title(f"{result['lineage']} Closest Neighbor (n = {len(result['data'])})", 
                 fontsize=title_size*0.8, fontweight='bold')
    
    # Style the plot
    apply_panel_styling(ax2, axis_text_size)
    
    # Add main title
    plt.suptitle(f"{result['lineage']} - SNP Analysis", 
                 fontsize=title_size, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    return fig

def set_plot_theme(theme_name):
    """Set the plot theme based on user input"""
    if theme_name == 'light':
        plt.style.use('seaborn-v0_8-whitegrid')
    elif theme_name == 'minimal':
        plt.style.use('seaborn-v0_8-white')
    elif theme_name == 'classic':
        plt.style.use('classic')
    elif theme_name == 'gray':
        plt.style.use('seaborn-v0_8-darkgrid')
    elif theme_name == 'dark':
        # Start with default style
        plt.style.use('default')
        
        # Create a custom dark theme that resembles the R plot
        plt.rcParams['figure.facecolor'] = '#4d4d4d'  # Medium-dark gray background
        plt.rcParams['axes.facecolor'] = '#4d4d4d'    # Same gray for plot area
        plt.rcParams['axes.edgecolor'] = 'white'      # White edges
        plt.rcParams['axes.labelcolor'] = 'white'     # White labels
        plt.rcParams['axes.grid'] = True              # Show grid
        plt.rcParams['grid.color'] = 'white'          # White grid
        plt.rcParams['grid.alpha'] = 0.3              # Subtle grid
        plt.rcParams['grid.linestyle'] = '-'          # Solid grid lines
        plt.rcParams['xtick.color'] = 'white'         # White ticks
        plt.rcParams['ytick.color'] = 'white'         # White ticks
        plt.rcParams['text.color'] = 'white'          # White text
        plt.rcParams['savefig.facecolor'] = '#4d4d4d' # Same gray for saved files
    elif theme_name == 'bw':
        plt.style.use('bmh')
    else:
        plt.style.use('seaborn-v0_8-whitegrid')  # Default fallback

def apply_panel_styling(ax, axis_text_size):
    """Apply consistent panel styling to axis"""
    ax.tick_params(labelsize=axis_text_size)
    
    # Format the grid and border
    if plt.rcParams['axes.facecolor'] == '#4d4d4d':  # If using our custom dark theme
        ax.grid(color='white', linestyle='-', linewidth=0.5, alpha=0.3)
        border_color = 'white'
    else:
        ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
        # Determine if we're using a dark background
        is_dark_theme = plt.rcParams['axes.facecolor'] == 'black'
        # Set the border color appropriately
        border_color = 'white' if is_dark_theme else 'black'
    
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(border_color)
        spine.set_linewidth(1)

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Auto-detect lineage from filename if not specified
    if args.lineage == "Unknown":
        detected_lineage = detect_lineage(os.path.basename(args.input_fasta))
        if detected_lineage != "Unknown":
            args.lineage = detected_lineage
            print(f"Auto-detected lineage: {detected_lineage}")
    
    # Display selected parameters
    print("\nKernel Plot Tool - Command Line Version")
    print("=============================================")
    print(f"Input file: {args.input_fasta}")
    print(f"Lineage: {args.lineage}")
    print(f"Output directory: {args.outputdir if args.outputdir else 'Auto (next to input file)'}")
    print(f"Density plot X-limit: {args.density_xlim}")
    print(f"Neighbor plot X-limit: {args.neighbor_xlim}")
    print(f"Bin size: {args.bin_size}")
    print(f"Plot dimensions: {args.width}x{args.height} inches")
    print(f"Show annotations: {not args.no_annotations}")
    print(f"Theme: {args.theme}")
    print("=============================================")
    
    try:
        print("Starting SNP distance calculation...")
        
        # Process FASTA file
        result = process_fasta(args.input_fasta, args.lineage, args.outputdir)
        
        print("Creating density plot...")
        # Create density plot
        density_fig, density_ax = create_density_plot(
            result,
            xlim_max=args.density_xlim,
            show_annotations=not args.no_annotations,
            title_size=args.title_size,
            axis_title_size=args.axis_title_size,
            axis_text_size=args.axis_text_size,
            annotation_size=args.annotation_size,
            plot_theme=args.theme
        )
        
        print("Creating closest neighbor plot...")
        # Create neighbor plot
        neighbor_fig, neighbor_ax = create_closest_neighbor_plot(
            result,
            xlim_max=args.neighbor_xlim,
            bin_size=args.bin_size,
            show_annotations=not args.no_annotations,
            title_size=args.title_size,
            axis_title_size=args.axis_title_size,
            axis_text_size=args.axis_text_size,
            annotation_size=args.annotation_size,
            plot_theme=args.theme
        )
        
        print("Creating combined figure...")
        # Create combined figure using our new function
        combined_fig = create_combined_figure(
            result,
            density_xlim=args.density_xlim,
            neighbor_xlim=args.neighbor_xlim,
            bin_size=args.bin_size,
            show_annotations=not args.no_annotations,
            title_size=args.title_size,
            axis_title_size=args.axis_title_size,
            axis_text_size=args.axis_text_size,
            annotation_size=args.annotation_size,
            plot_theme=args.theme,
            plot_width=args.width,
            plot_height=args.height
        )
        
        # Calculate bin counts and save to CSV
        print("Calculating bin statistics...")
        my_data_num = result['lower_triangle_values']
        density_breaks = np.arange(0, args.density_xlim + args.bin_size, args.bin_size)
        density_bin_info = calculate_bin_counts(my_data_num, density_breaks, include_stats=True)
        
        # Convert to float type to handle NaN values
        data = result['data'].values.astype(float)
        np.fill_diagonal(data, np.nan)
        lowest_distances = np.nanmin(data, axis=1)
        neighbor_breaks = np.arange(0, args.neighbor_xlim + args.bin_size, args.bin_size)
        neighbor_bin_info = calculate_bin_counts(lowest_distances, neighbor_breaks, include_stats=True)
        
        # Save bin counts as CSV
        density_bins_file = os.path.join(result['output_folder'], f"{result['input_file_name']}_density_bins.csv")
        density_bin_info['data'].to_csv(density_bins_file, index=False)
        
        neighbor_bins_file = os.path.join(result['output_folder'], f"{result['input_file_name']}_neighbor_bins.csv")
        neighbor_bin_info['data'].to_csv(neighbor_bins_file, index=False)
        
        # Save statistics as text files
        density_stats_file = os.path.join(result['output_folder'], f"{result['input_file_name']}_density_stats.txt")
        with open(density_stats_file, 'w') as f:
            f.write(density_bin_info['text'])
        
        neighbor_stats_file = os.path.join(result['output_folder'], f"{result['input_file_name']}_neighbor_stats.txt")
        with open(neighbor_stats_file, 'w') as f:
            f.write(neighbor_bin_info['text'])
        
        # Print density statistics
        print("\nSNP Distance Statistics:")
        print(density_bin_info['text'])
        print("\n")
        
        # Print neighbor statistics
        print("Closest Neighbor Statistics:")
        print(neighbor_bin_info['text'])
        print("\n")
        
        # Save the plots
        print("Saving plots...")
        
        # Save individual plots
        pdf_file = os.path.join(result['output_folder'], f"{result['input_file_name']}_density_plot.pdf")
        density_fig.savefig(pdf_file, dpi=300, bbox_inches='tight')
        
        hist_file = os.path.join(result['output_folder'], f"{result['input_file_name']}_closest_neighbor.pdf")
        neighbor_fig.savefig(hist_file, dpi=300, bbox_inches='tight')
        
        # Save combined figure
        combined_file = os.path.join(result['output_folder'], f"{result['input_file_name']}_combined_figure.pdf")
        combined_fig.savefig(combined_file, dpi=300, bbox_inches='tight')
        
        print(f"All plots saved to {result['output_folder']}")
        print(f"Combined figure saved as: {combined_file}")
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Exit with success
    sys.exit(0)

if __name__ == "__main__":
    main()