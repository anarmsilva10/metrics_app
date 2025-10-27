import pysam
import pysam.samtools
import argparse
import pandas as pd
import numpy as np
import os
import json
from log import shared_log

def initialize_metrics():
    """
    Initialize a dictionary with default quality metrics.

    Returns:
        dict: Dictionary containing default metrics for sequencing statistics set to zero.
    """

    metrics = {
        'Total Reads': 0,
        'Total Bases': 0, # the giving metric doesn't match any sample of the validation dataset
        'Average Base Quality': 0,
        'Read Length (Max)': 0,
        'Reads Mapped (%)': 0,
        'Mapping Quality 0 (%)': 0,
        'Error rate': 0,
        'Insert Size (Avg)': 0, # the giving metric doesn't match any sample of the validation dataset
        'Insert Size (Std)': 0,  # the giving metric doesn't match any sample of the validation dataset
        'Reads paired (%)': 0,
        'Reads properly paired (%)': 0,
        'Duplicate Reads (%)': 0,
        'Proportion of duplicated reads': 0,
    }

    return metrics

def run_samtools_stats(cram_path, ref_cache=None, insert_size_filter=None, most_inserts=None):
    """
    Extract quality metrics using pysam.samtools.stats.

    Args:
        cram_path: path to the input CRAM/BAM file.
        ref_cache: path to htslib ref cache (optional).
        insert_size_filter: tuple (min, max) to filter insert sizes (optional)
        most_inserts: float to limit the number of insert sizes (optional)

    Returns:
        pd.DataFrame: DataFrame containing sequencing statistics metrics.
    """
    metrics = initialize_metrics()
    # Fields that are in the pysam.samtools.stats that are needed to calculate/retrieve metrics
    sn_fields = {
    'sequences': 'Total Reads',
    'bases mapped': 'Total Bases',
    'average quality': 'Average Base Quality',
    'error rate': 'Error rate',
    'percentage of properly paired reads (%)': 'Reads properly paired (%)'
    }
    raw_values = {}
    max_read_length = 0
    insert_sizes = []

    if ref_cache:
        os.environ['REF_CACHE'] = os.path.join(ref_cache, '%2s/%2s/%s')

    output = pysam.samtools.stats(cram_path)

    if output:
        for line in output.splitlines():
            if line.startswith('SN'):
                columns = line.strip().split('\t')
                if len(columns) >= 3:
                    raw_data = columns[1].strip().rstrip(':')
                    value = columns[2].strip()

                    # Ensure that the output is storage as numeric instead as text
                    try:
                        value = float(value) if '.' in value else int(value)
                    except ValueError:
                        continue

                    # Store for raw fields
                    raw_values[raw_data] = value

                    # Store metrics that can be retrieve directly
                    if raw_data in sn_fields and sn_fields[raw_data]:
                        metrics[sn_fields[raw_data]] = value

            elif line.startswith('RL'):
                columns = line.strip().split('\t')
                if len(columns) >= 2:
                    try:
                        read_length = int(columns[1])
                        if read_length > max_read_length:
                            max_read_length = read_length
                    except ValueError:
                        continue

            elif line.startswith('IS'):
                columns = line.strip().split('\t')
                if len(columns) >= 3:
                    try:
                        insert_size = int(columns[1])
                        count = int(columns[2])
                        if insert_size_filter:
                            min_sz, max_sz = insert_size_filter
                            if min_sz <= insert_size <= max_sz:
                                insert_sizes.extend([insert_size] * count)
                        else:
                            insert_sizes.extend([insert_size] * count)
                    except ValueError:
                        continue

        # Apply insert size percentile filtering if requested
        if insert_sizes and most_inserts is not None:
            insert_sizes = np.array(insert_sizes)
            cutoff = np.quantile(insert_sizes, most_inserts)
            insert_sizes = insert_sizes[insert_sizes <= cutoff]
            insert_sizes = insert_sizes.tolist()

        metrics['Read Length (Max)'] = max_read_length
        metrics['Error rate'] = round(metrics['Error rate'], 6)
        if insert_sizes:
            metrics['Insert Size (Avg)'] = round(np.mean(insert_sizes), 2)
            metrics['Insert Size (Std)'] = round(np.std(insert_sizes), 2)

    # Derived metrics
    try:
        if 'reads mapped' in raw_values and 'raw total sequences' in raw_values:
            metrics['Reads Mapped (%)'] = round(
                raw_values['reads mapped'] / raw_values['raw total sequences'] * 100, 2
            )
        if 'reads MQ0' in raw_values and 'raw total sequences' in raw_values:
            metrics['Mapping Quality 0 (%)'] = round(
                raw_values['reads MQ0'] / raw_values['raw total sequences'] * 100, 4
            )
        if 'reads duplicated' in raw_values and 'raw total sequences' in raw_values:
            metrics['Duplicate Reads (%)'] = round(
                raw_values['reads duplicated'] / raw_values['raw total sequences'] * 100, 2
            )
            metrics['Proportion of duplicated reads'] = round(
                raw_values['reads duplicated'] / raw_values['raw total sequences'], 2
            )
        if 'reads paired' in raw_values and 'raw total sequences' in raw_values:
            metrics['Reads paired (%)'] = round(
                raw_values['reads paired'] / raw_values['raw total sequences'] * 100, 2
            )

    except ZeroDivisionError:
        pass

    return pd.DataFrame([metrics])


def run_pysam_stats(cram_path, ref_cache=None, min_insert_size=1, max_insert_size=8000):
    """
    [DEBUG/TESTING ONLY] Extract quality metrics using pure Python iteration over reads.
    This function relies on pysam.AlignmentFile. Therefore it is much slower.
    Intended only for debugging and validation, not for production use.

    Args:
        cram_path: path to the input CRAM/BAM file.
        ref_cache: path to htslib ref cache (optional).
        min_insert_size_: float for minimum insert sizes (optional). Default 1.
        max_insert_size: float for maximum insert sizes (optional). Default 8000.

    Returns:
        pd.DataFrame: DataFrame containing sequencing statistics metrics.
    """

    metrics = initialize_metrics()

    # Counters
    total_reads = 0
    mapped_reads = 0
    mq0_reads = 0
    mapq_considered = 0
    duplicate_reads = 0
    paired_reads = 0
    proper_paired_reads = 0
    total_bases = 0
    mismatches = 0
    mapped_bases = 0
    max_read_length = 0

    # Base quality accumulators
    baseq_sum = 0
    baseq_count = 0

    # Insert size accumulators (Welford’s algorithm)
    insert_count = 0
    insert_mean = 0.0
    insert_M2 = 0.0

    with pysam.AlignmentFile(cram_path, "rc" if cram_path.endswith(".cram") else "rb", reference_filename=ref_cache) as aln:
        for read in aln.fetch(until_eof=True):
            # Skip secondary/supplementary like samtools stats
            if read.is_secondary or read.is_supplementary:
                continue

            total_reads += 1
            read_length = read.query_length or 0
            total_bases += read_length
            max_read_length = max(max_read_length, read_length)

            if not read.is_unmapped:
                mapped_reads += 1

                # Count mapped bases (excluding softclips/hardclips)
                mapped_bases += read.query_alignment_length or 0

                # Count mismatches using NM tag if present
                try:
                    mismatches += read.get_tag("NM")
                except KeyError:
                    pass

            # MAPQ0 calculation: only for mapped, not duplicate, not QC-fail
            if (not read.is_unmapped) and (not read.is_duplicate) and (not read.is_qcfail):
                mapq_considered += 1
                if read.mapping_quality == 0:
                    mq0_reads += 1

            if read.is_duplicate:
                duplicate_reads += 1

            if read.is_paired:
                paired_reads += 1
                if read.is_proper_pair:
                    proper_paired_reads += 1
                # Only count template lengths within allowed insert_size
                tlen = abs(read.template_length)
                if min_insert_size <= tlen < max_insert_size:
                    insert_count += 1
                    delta = tlen - insert_mean
                    insert_mean += delta / insert_count
                    insert_M2 += delta * (tlen - insert_mean)

            if read.query_qualities is not None:
                baseq_sum += sum(read.query_qualities)
                baseq_count += len(read.query_qualities)

    # Fill in metrics
    if total_reads > 0:
        metrics['Total Reads'] = total_reads
        metrics['Total Bases'] = total_bases
        metrics['Reads Mapped (%)'] = round(mapped_reads / total_reads * 100, 2)
        metrics['Duplicate Reads (%)'] = round(duplicate_reads / total_reads * 100, 2)
        metrics['Proportion of duplicated reads'] = round(duplicate_reads / total_reads, 2)
        metrics['Reads paired (%)'] = round(paired_reads / total_reads * 100, 2)
        metrics['Reads properly paired (%)'] = round(proper_paired_reads / total_reads * 100, 2)
        metrics['Read Length (Max)'] = max_read_length

    if mapq_considered > 0:
        metrics['Mapping Quality 0 (%)'] = round(mq0_reads / mapq_considered * 100, 4)

    if baseq_count > 0:
        metrics['Average Base Quality'] = round(baseq_sum / baseq_count, 2)

    if insert_count > 1:
        metrics['Insert Size (Avg)'] = round(insert_mean, 2)
        metrics['Insert Size (Std)'] = round((insert_M2 / (insert_count - 1))**0.5, 2)

    if mapped_bases > 0:
        metrics['Error rate'] = round(mismatches / mapped_bases, 6)

    return pd.DataFrame([metrics])


def calculate_stats(cram_path, ref_cache=None, method_stats="run_samtools_stats"):
    """
    Calculate statistical metrics.

    Args:
        cram_path: path to the input CRAM/BAM file.
        ref_cache: path to htslib ref cache (optional).
        method_stats: method to extract statistical metrics. Default is run_samtools_stats.

    Returns:
        pd.DataFrame: DataFrame containing sequencing statistics metrics.    
    """
    
    shared_log.logger.info(f"Using BAM/CRAM file: {cram_path} with method {method_stats}.")

    if method_stats == "run_pysam_stats":
        return run_pysam_stats(cram_path, ref_cache, min_insert_size=1, max_insert_size=8000)
    else:
        return run_samtools_stats(cram_path, ref_cache)


def main():

    parser = argparse.ArgumentParser(description="Calculate quality metrics from a CRAM/BAM file.")

    # Mandatory arguments
    parser.add_argument('cram_path', type=str, help="Path to the input CRAM/BAM file.")
    parser.add_argument('output', type=str, help='Path to save output file (Excel or JSON).')

    # Optional arguments
    parser.add_argument('--json', action='store_true', help='Output in JSON format (default is Excel).')
    parser.add_argument('--ref_cache', type=str, help="Path to htslib ref cache.")
    parser.add_argument('--method_stats', choices=['run_samtools_stats', 'run_pysam_stats'], default='run_samtools_stats', help="Choose the method to extract stats metrics (default: run_samtools_stats).")

    args = parser.parse_args()

    # Calculate metrics
    metrics_df = calculate_stats(args.cram_path, args.ref_cache, args.method_stats)
    if metrics_df is not None:
        if args.json:
            cram_name = os.path.splitext(os.path.basename(args.cram_path))[0]
            with open(args.output, 'w') as file:
                json.dump({cram_name: metrics_df.iloc[0].to_dict()}, file, indent=4)
        else:
            metrics_df.to_excel(args.output, index=False)
        shared_log.logger.info(f"Metrics saved successfully to {args.output}")
    else:
        shared_log.logger.error("Error: Metrics calculation failed.")

if __name__ == '__main__':
    main()