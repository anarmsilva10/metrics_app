import pandas as pd
import argparse
import os
import json
from backend import depth, samtools_stats, plot
from log import shared_log

# Metrics Initialization
def initialize_metrics(thresholds):
    """
    Initialize a dictionary of default coverage metrics.

    Args:
        thresholds (list[int]): Coverage thresholds to include in metrics (e.g., [1, 10, 20, 30]).

    Returns:
        dict: Dictionary with initialized coverage metrics set to zero.
    """

    metrics = {
        'Size Coding': 0,
        'Size Covered': 0,
        'Breadth of Coverage %': 0,
        'Average Read Depth': 0,
        'Min Read Depth': 0,
        'Max Read Depth': 0,
        'Median Read Depth': 0
    }

    # Initialize coverage thresholds metrics
    for t in thresholds:
        metrics[f'Depth of Coverage {t}x'] = 0

    metrics[f'Depth of Coverage (> {thresholds[-1]}x)'] = 0
    metrics[f'Depth of Coverage (0-{thresholds[0]}x)'] = 0

    for i in range(len(thresholds) - 1):
        lower = thresholds[i]
        upper = thresholds[i + 1]

        if upper:
            metrics[f'Depth of Coverage ({lower + 1}-{upper}x)'] = 0

    return metrics

# Coverage Metrics Calculation
def calculate_coverage_metrics(depth_series, size_coding, thresholds):
    """
    Calculate coverage statistics from a depth series.

    Args:
        depth_series (pd.Series): Depth values across coding regions.
        size_coding (int): Total size of coding regions in bases.
        thresholds (list[int]): Coverage thresholds for evaluation.

    Returns:
        pd.DataFrame: DataFrame containing computed coverage metrics, including depth distribution and summary statistics.
    """

    metrics_df = pd.DataFrame([initialize_metrics(thresholds)])

    total_positions = len(depth_series)
    
    if total_positions == 0 or size_coding == 0:
        return metrics_df

    # Metric Calculation and Storage in a Dataframe
    metrics_df['Size Coding'] = size_coding
    # Calculate covered size
    metrics_df['Size Covered'] = depth_series.count()
    # Average read depth
    metrics_df['Average Read Depth'] = depth_series.mean().round(2)
    # Min, max and median read depth
    metrics_df['Min Read Depth'] = depth_series.min()
    metrics_df['Max Read Depth'] = depth_series.max()
    metrics_df['Median Read Depth'] = depth_series.median() # Was add because of new metrics implementation
    # Breadth of coverage (%)
    metrics_df['Breadth of Coverage %'] = (metrics_df['Size Covered'] / metrics_df['Size Coding']) * 100

    # Calculate depth of coverage for each threshold
    for t in thresholds:
        metrics_df[f'Depth of Coverage {t}x'] = (depth_series >= t).sum() / total_positions * 100

    for i in range(len(thresholds)):
        lower = thresholds[i]
        upper = thresholds[i + 1] if i + 1 < len(thresholds) else None

        if upper:
            metrics_df[f'Depth of Coverage ({lower + 1}-{upper}x)'] = ((depth_series > lower) & (depth_series <= upper)).sum()
        else:
            metrics_df[f'Depth of Coverage (> {lower}x)'] = (depth_series > lower).sum()

    if metrics_df.empty:
        shared_log.logger.error('Error calculating metrics.')
    
    return metrics_df

def calculate_gene_metrics(bed_df, depth_df, thresholds):
    """
    Compute coverage metrics per gene.

    Args:
        bed_df (pd.DataFrame): BED file data containing genomic regions, with columns ['GENE', 'START', 'END', 'SIZE'].
        depth_df (pd.DataFrame): Depth data with columns['POSITION', 'DEPTH'].
        thresholds (list[int]): Coverage thresholds for evaluation.

    Returns:
        pd.DataFrame: DataFrame with metrics for each gene.
    """

    gene_data = []
    genes = bed_df['GENE'].unique()

    for gene in genes:
        gene_bed_df = bed_df[bed_df['GENE'] == gene]
        size_coding_gene = gene_bed_df['SIZE'].sum()

        gene_depths = pd.concat([
            depth_df[(depth_df['POSITION'] >= row['START']) & (depth_df['POSITION'] <= row['END'])]['DEPTH']
            for _, row in gene_bed_df.iterrows()
        ], ignore_index=True)

        gene_metrics = calculate_coverage_metrics(gene_depths, size_coding_gene, thresholds)
        gene_metrics.insert(0, 'Gene', gene)
        gene_data.append(gene_metrics)

    gene_metrics_df = pd.concat(gene_data, ignore_index=True) if gene_data else pd.DataFrame()

    return gene_metrics_df

def calculate_exon_metrics(bed_df, depth_df, thresholds):
    """
    Compute coverage metrics per exon.

    Args:
        bed_df (pd.DataFrame): BED file data containing genomic regions, with columns ['EXON', 'START', 'END', 'SIZE'].
        depth_df (pd.DataFrame): Depth data with columns ['POSITION', 'DEPTH'].
        thresholds (list[int]): Coverage thresholds for evaluation.

    Returns:
        pd.DataFrame: DataFrame with metrics for each exon.
    """

    exon_data = []

    for exon_id in bed_df['EXON'].unique():
        exon_bed_df = bed_df[bed_df['EXON'] == exon_id]
        size_coding_exon = exon_bed_df['SIZE'].sum()

        exon_depths = pd.concat([
            depth_df[(depth_df['POSITION'] >= row['START']) & (depth_df['POSITION'] <= row['END'])]['DEPTH']
            for _, row in exon_bed_df.iterrows()
        ], ignore_index=True)

        exon_metrics = calculate_coverage_metrics(exon_depths, size_coding_exon, thresholds)
        exon_metrics.insert(0, 'Exon', exon_id)
        exon_data.append(exon_metrics)

    exon_metrics_df = pd.concat(exon_data, ignore_index=True) if exon_data else pd.DataFrame()

    return exon_metrics_df

def calculate_metrics(cram_path, bed_df, gene_selection=None, exon_selection=None, method='run_pysam_coverage', thresholds = [1, 10, 15, 20, 30, 50, 100, 500], stats=False, ref_cache=None, method_stats='run_samtools_stats'):
    """
    Compute coverage and statistical (optional) metrics using a CRAM/BAM file and a BED file.

    Args:
        cram_path (str): Path to the CRAM/BAM input file.
        bed_df (pd.DataFrame): BED dataframe containing genomic regions.
        gene_selection (list[str], optional): List of selected genes to analyze.
        exon_selection (list[int], optional): List of selected exons to analyze.
        method (str, optional): Method for extracting depth ('run_samtools_depth', 'run_pysam_depth', 'run_pysam_coverage'). Default is 'run_pysam_coverage'.
        thresholds (list[int], optional): Coverage thresholds. Default is [1, 10, 15, 20, 30, 50, 100, 500].
        stats (bool, optional): Whether to include additional samtools stats metrics. Default is False.
        ref_cache (str, optional): Path to htslib ref cache.
        method_stats (str, optional): Method for extracting stats metrics ('run_samtools_stats', 'run_pysam_stats'). Default is 'run_samtools_stats'.

    Returns:
        pd.DataFrame | None: DataFrame containing merged coverage metrics and optional statistics, or None if no depth data found.
    """

    if ref_cache:
        os.environ['REF_CACHE'] = os.path.join(ref_cache, '%2s/%2s/%s')

    shared_log.logger.info(f"Using BAM/CRAM file: {cram_path} to calculate metrics.")

    if stats:
        shared_log.logger.info(f'Using stats extraction method: {method_stats}')
        metrics_df =  samtools_stats.calculate_stats(cram_path, method_stats=method_stats, ref_cache=ref_cache)
    else:
        metrics_df = pd.DataFrame()

    bed_df['SIZE'] = bed_df['END'] - bed_df['START']
    bed_df = depth.handling_files(cram_path, bed_df)

    shared_log.logger.info(f'Using depth extraction method: {method}')
    
    call =  getattr(depth, method)
    depth_df = call(cram_path, bed_df)

    if depth_df.empty:
        shared_log.logger.warning('No depth data found.')
        return None
    
    if exon_selection:
        coverage_df = calculate_exon_metrics(bed_df, depth_df, thresholds)
    elif gene_selection:
        coverage_df = calculate_gene_metrics(bed_df, depth_df, thresholds)
    else:
        coverage_df = calculate_coverage_metrics(depth_df['DEPTH'], bed_df['SIZE'].sum(), thresholds)

    all_metrics_df = pd.concat([coverage_df.reset_index(drop=True), metrics_df.reset_index(drop=True)], axis=1) 

    return all_metrics_df

# Parsing arguments
def main():
    parser = argparse.ArgumentParser(description="Calculate depth of coverage from CRAM/BAM and BED files.")

    # Mandatory arguments
    parser.add_argument('cram_path', type=str, help="Path to the input CRAM/BAM file.")
    parser.add_argument('bed_path', type=str, help="Path to the input BED file.")
    parser.add_argument('output', type=str, help='Path to save output file.')

    # Optional arguments
    parser.add_argument('--json', help='Output in JSON format.')
    parser.add_argument('--gene_selection', type=str, nargs='+', help="List of genes to include (optional).")
    parser.add_argument('--exon_selection', type=int, nargs='+', help="List of exons to include (optional).")
    parser.add_argument('--log_name', type=str, help="Base name of thr log file (default: log to stdout).")
    parser.add_argument('--output_plot', type=str, help="Path to save coverage plot.")
    parser.add_argument('--ref_cache', type=str, help="Path to htslib ref cache.")
    parser.add_argument('--stats', action='store_true', help="Add statistics metrics.")

    # Method argument - It isn't mandatory, it has a default method
    parser.add_argument('--method', choices=['run_samtools_depth', 'run_pysam_depth', 'run_pysam_coverage'], default='run_pysam_coverage', help="Choose the method to extract depth (default: run_pysam_coverage).")
    parser.add_argument('--method_stats', choices=['run_samtools_stats', 'run_pysam_stats'], default='run_samtools_stats', help="Choose the method to extract stats metrics (default: run_samtools_stats).")

    # Parse arguments
    args = parser.parse_args()

    if args.log_name:
        shared_log.logger = shared_log.Logging(basename=args.log_name, foldername=__path__).get_logger()
    
    shared_log.logger.info(f"Using BED file: {args.bed_path} and BAM/CRAM file: {args.cram_path} for coverage calculation metrics.")
    
    bed_df = depth.filter_bed(args.bed_path, args.gene_selection, args.exon_selection)

    # Calculate metrics
    metrics_df = calculate_metrics(args.cram_path, bed_df, args.gene_selection, args.exon_selection, method=args.method, stats=args.stats, ref_cache=args.ref_cache, method_stats=args.method_stats)
    
    if metrics_df is not None:

        if args.json is not None:
            cram_name = os.path.splitext(os.path.basename(args.cram_path))[0]
            bed_name = os.path.splitext(os.path.basename(args.bed_path))[0]
            sample = f"{cram_name}_{bed_name}"

            with open(args.output, 'w') as file:
                json.dump({sample: metrics_df.iloc[0].to_dict()}, file, indent=4)
        else:
            metrics_df.to_excel(args.output, index=False)
        shared_log.logger.info(f"Metrics saved successfully to {args.output}")
    else:
        shared_log.logger.error("Error: Metrics calculation failed.")

    if args.output_plot:
        fig = plot.coverage_plot(args.cram_path, bed_df, args.method)
        fig.write_image(args.output_plot)
        shared_log.logger.info(f"Plot saved in {args.output_plot}.")

if __name__ == '__main__':
    main()