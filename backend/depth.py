import argparse
import os
import subprocess
import pandas as pd
import pysam
import io
import tempfile
import re
from log import shared_log

def validate_paths(cram_path, bed_path):
    """
    Validate that the provided CRAM/BAM and BED file paths exist and have correct extensions.

    Args:
        cram_path (str): Path to CRAM or BAM file.
        bed_path (str): Path to BED file.

    Returns:
        bool: True if paths and extensions are valid, False otherwise.
    """

    # Ensure valide file path
    if not os.path.isfile(cram_path) or not os.path.isfile(bed_path):
        shared_log.logger.error(f'Invalid file path provided for CRAM/BAM ({cram_path}) or BED ({bed_path}) file.')
        return False

    # Ensure valid file extensions
    if not (cram_path.endswith('.cram') or cram_path.endswith('.bam')):
        shared_log.logger.error(f'Unsupported file type for {cram_path}. Must be .cram or .bam file.')
        return False

    if not bed_path.endswith('.bed'):
        shared_log.logger.error(f'Invalid BED file: {bed_path}. Must be a .bed file.')
        return False

    return True

def validate_samtools():
    """
    Validate that Samtools is installed and accessible in the system PATH.

    Raises:
        FileNotFoundError: If Samtools is not found.
        subprocess.CalledProcessError: If Samtools fails to run.
    """

    try:
        subprocess.run(["samtools", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    except FileNotFoundError:
        shared_log.logger.critical("Error: 'samtools' is not installed or not found in PATH.")
        raise
    except subprocess.CalledProcessError:
        shared_log.logger.critical("Error: Unable to validate 'samtools' installation.")
        raise

def handling_files(cram_path, bed_df):
    """
    Adjust chromosome naming in BED file to match CRAM/BAM header.

    Args:
        cram_path (str): Path to CRAM or BAM file.
        bed_df (pd.DataFrame): DataFrame containing BED file data with 'CHROM' column.

    Returns:
        pd.DataFrame | None: Adjusted BED DataFrame or None if empty.
    """

    # Validate file extension
    if cram_path.endswith('.cram'):
        mode = 'rc'
    elif cram_path.endswith('.bam'):
        mode = 'rb'
    else:
        shared_log.logger.error('Unnsuported file type. Must be .cram or .bam.')

    with pysam.AlignmentFile(cram_path, mode) as cram_file:
        header = cram_file.header.to_dict()

        # Check if chr are in CRAM/BAM header
        match = any(re.search(r'^chr[0-9XY]+', sq['SN']) for sq in header.get('SQ', []))

        if match:
            # Add chr to bed_df if missing
            bed_df['CHROM'] = bed_df['CHROM'].astype(str).str.replace(r'^(?!chr)(\d+|Y|X)', r'chr\1', regex=True)
        else:
            # Remove chr from bed_df if present
            bed_df['CHROM'] = bed_df['CHROM'].astype(str).str.replace(r'^chr', '', regex=True)

    return bed_df

def filter_bed(bed_path, gene_selection=None, exon_selection=None):
    """
    Filter BED file based on selected genes or exons.

    Args:
        bed_path (str): Path to BED file.
        gene_selection (list[str] | str, optional): Genes to include.
        exon_selection (list[int] | int, optional): Exons to include.

    Returns:
        pd.DataFrame | None: Filtered BED DataFrame, or None if no matching regions.
    """

    try:

        if gene_selection or exon_selection:
            columns = ['CHROM', 'START', 'END', 'GENE', 'EXON', 'SIZE']
        else:
            columns = ['CHROM', 'START', 'END']
        bed_df = pd.read_csv(bed_path, sep='\t', header=None, usecols=list(range(len(columns))))
        bed_df.columns = columns

        if gene_selection:
            shared_log.logger.info(f'Filtering genes: {gene_selection}')

            if isinstance(gene_selection, str):
                gene_selection = [gene_selection]
            bed_df = bed_df[bed_df['GENE'].isin(gene_selection)]
        
        if exon_selection:
            shared_log.logger.info(f'Filtering exons: {exon_selection}')

            if isinstance(exon_selection, str):
                exon_selection = [exon_selection]
            bed_df = bed_df[bed_df['EXON'].isin(exon_selection)]

        if bed_df.empty:
            shared_log.logger.warning('No matching regions found for the provided gene or exon selection.')
            return None
        else:
            shared_log.logger.info(f'Bed successfully filtered. Total rows: {len(bed_df)}.')

        return bed_df
    except Exception as e:
        shared_log.logger.error(f'Error filtering BED file: {e}')
        raise

def run_samtools_depth(cram_path, bed_df):
    """
    Runs samtools depth on a given CRAM/BAM file using a filtered BED file.
    Returns samtools output as a dataframe.

    Args:
        cram_path (str): Path to CRAM/BAM file.
        bed_df (pd.DataFrame): Filtered BED DataFrame.

    Returns:
        pd.DataFrame: DataFrame with columns ['CHROM', 'POSITION', 'DEPTH'].
    """

    validate_samtools()

    try:
        samtools_command = ['samtools', 'depth', '-aa', '-b', '-', cram_path]
        samtools_process = subprocess.Popen(samtools_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Pass filtered BED content as input
        filtered_bed = bed_df.to_csv(sep='\t', index=False, header=False)
        samtools_output, samtools_error = samtools_process.communicate(input=filtered_bed)

        if samtools_error:
            shared_log.logger.error(f'Samtools error: {samtools_error.strip()}.')
        
        if samtools_output:
            # Convert the output string to a DataFrame directly
            depth_df = pd.read_csv(io.StringIO(samtools_output), sep='\t', header=None, names=['CHROM', 'POSITION', 'DEPTH'])
            return depth_df
        else:
            shared_log.logger.warning('No output from samtools.')
            return pd.DataFrame()

    except Exception as e:
        shared_log.logger.error(f'Error running samtools: {e}')
        raise

def run_pysam_depth(cram_path, bed_df):
    """
    Runs pysam.depth on a given CRAM/BAM file using a filtered BED file.
    Returns depth output as a dataframe.

    Args:
        cram_path (str): Path to CRAM/BAM file.
        bed_df (pd.DataFrame): Filtered BED DataFrame.

    Returns:
        pd.DataFrame: DataFrame with columns ['CHROM', 'POSITION', 'DEPTH'].
    """
   
    try:
     
        with tempfile.NamedTemporaryFile('w', suffix='.bed') as bed_file:
            bed_df.to_csv(bed_file.name, sep='\t', index=False, header=False)

            output = pysam.depth('-a', '-b', bed_file.name, cram_path)

        if output:
            depth_df = pd.read_csv(io.StringIO(output), sep='\t', header=None, names=['CHROM', 'POSITION', 'DEPTH'])
            return depth_df
        else:
            shared_log.logger.warning('No depth data returned from pysam.')
            return pd.DataFrame()

    except Exception as e:
        shared_log.logger.error(f'Error extracting depth with pysam: {e}')
        raise

def run_pysam_coverage(cram_path, bed_df):
    """
    Run pysam coverage on the regions defined in the BED file and return the depth output.

    Args:
        cram_path (str): Path to CRAM/BAM file.
        bed_df (pd.DataFrame): Filtered BED DataFrame.

    Returns:
        pd.DataFrame: DataFrame with columns ['CHROM', 'POSITION', 'DEPTH'].
    """

    if cram_path.endswith('.cram'):
        index_path = cram_path + '.crai'
        mode = 'rc'
    elif cram_path.endswith('.bam'):
        index_path = cram_path + '.bai'
        mode = 'rb'
    else:
        shared_log.logger.error('Unsupported file type. File must be .cram or .bam')
        
    # Check if index file exists, if not, generate it
    if not os.path.exists(index_path):
        shared_log.logger.warning(f'Index file {index_path} not found. Generating index...')
        pysam.index(cram_path)
        shared_log.logger.info(f'Index file generated: {index_path}')
    else:
        shared_log.logger.info(f'Index file already exist: {index_path}')
        
    with pysam.AlignmentFile(cram_path, mode) as cram_file:
        output = []
            
        for _, row in bed_df.iterrows():
            chrom, start, end = str(row['CHROM']), row['START'], row['END']
                
            coverage = cram_file.count_coverage(chrom, start, end, quality_threshold = 0)
            total_depth = [sum(base_counts) for base_counts in zip(*coverage)]

            for pos, depth in enumerate(total_depth, start=start):
                output.append((chrom, pos+1, depth))

    if output:
        depth_df = pd.DataFrame(output, columns=['CHROM', 'POSITION', 'DEPTH'])
        return depth_df
    else:
        shared_log.logger.warning('No depth data returned from pysam.')
        return pd.DataFrame()

def get_depth_data(cram_path, bed_df, output_path = None, method='run_pysam_coverage', ref_cache=None):
    """
    Calculate depth of coverage using either samtools or pysam.
    Returns a DataFrame with depth output.

    Args:
        cram_path (str): Path to CRAM/BAM file.
        bed_df (pd.DataFrame): Filtered BED DataFrame.
        output_path (str, optional): Path to save depth data as TSV.
        method (str, optional): Method to use ('run_samtools_depth', 'run_pysam_depth', 'run_pysam_coverage'). Default is 'run_pysam_coverage'.
        ref_cache (str, optional): Path to reference cache for htslib.

    Returns:
        pd.DataFrame | None: Depth DataFrame, or None if no data found.
    """
    bed_df = handling_files(cram_path, bed_df)

    if ref_cache:
        os.environ['REF_CACHE'] = os.path.join(ref_cache, '%2s/%2s/%s')

    shared_log.logger.info(f'Using method: {method}.')

    call = globals().get(method)
    depth_df = call(cram_path, bed_df)
    
    if depth_df is None or depth_df.empty:
        shared_log.logger.warning('No depth data found.')
        return None
    
    if output_path:
        depth_df.to_csv(output_path, sep='\t', index=False, header=False)
        shared_log.logger.info(f'Depth data saved to {output_path}')

    shared_log.logger.info(f'Depth data successfully generated. Total rows: {len(depth_df)}')

    return depth_df

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate depth of coverage from CRAM/BAM and BED files.")

    # Mandatory arguments
    parser.add_argument('cram_path', type=str, help="Path to the input CRAM/BAM file.")
    parser.add_argument('bed_path', type=str, help="Path to the input BED file.")
    parser.add_argument('output_path', type=str, help="Path to the output file.")
    
    # Optional arguments
    parser.add_argument("--log_file", type=str, help="Path to a log file (default: log to stdout).")
    parser.add_argument('--gene_selection', type=str, nargs='+', help="List of genes to include (optional).")
    parser.add_argument('--exon_selection', type=int, nargs='+', help="List of exons to include (optional).")

    # Method argument - It isn't mandatory, it has a default method
    parser.add_argument('--method', choices=['run_samtools_depth', 'run_pysam_depth','run_pysam_coverage'], default='run_pysam_coverage', help="Choose the method to extract depth (default: run_pysam_coverage).")
    parser.add_argument('--ref_cache', type=str, help="Path to htslib ref cache")

    # Parse arguments
    args = parser.parse_args()

    if args.log_file:
        shared_log.logger = shared_log.Logging(filename=args.log_file, foldername=__path__).get_logger()

    if not validate_paths(args.cram_path, args.bed_path):
        raise
    
    shared_log.logger.info(f"Using BED file: {args.bed_path} and BAM/CRAM file: {args.cram_path}.")

    bed_df = filter_bed(args.bed_path, args.gene_selection, args.exon_selection)    

    get_depth_data(args.cram_path, bed_df, output_path=args.output_path, method=args.method, ref_cache=args.ref_cache)

if __name__ == "__main__":
    main()