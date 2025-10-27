from pathlib import Path
import pandas as pd
from backend.helper import extract_chr_value, convert_numeric_to_chr, load_json, exon_number_json, utr_info, padding_adding
from log import shared_log

data_dir = Path('data')
    
# Extract all the genes stored in the JSON file
def get_genes(genome):
    """
    Extract all gene symbols stored in a genome JSON file.

    Args:
        genome (str): Genome identifier (e.g., "hg38").

    Returns:
        list[str]: List of gene symbols present in the genome file.
    """

    shared_log.logger.info("Getting Gene Data.")

    data = load_json(genome)

    # gene symbol are the keys from each json
    return list(data.keys())

# For the gene selected, it builds the list of exons presented in the JSON file
def get_exons(genome, gene):
    """
    Retrieve exon information for a given gene.

    Args:
        genome (str): Genome identifier (e.g., "hg38").
        gene (str): Gene symbol to extract exons for.

    Returns:
        list[dict]: List of exons for the gene, where each exon dict contains:
            - gene (str): Gene symbol
            - chrom (str): Chromosome
            - start (int): Start position coordinate
            - end (int): End position coordinate
            - number (str): Exon number (e.g., "Exon 1")
            - size (int): Exon length in base pairs
    """

    shared_log.logger.info("Getting Exon Data.")

    data = load_json(genome)
    gene_data = data.get(gene)

    if not gene_data:
        return []

    chrom = gene_data['chromosome']

    exons = []

    for key, transcripts_list in gene_data.items():
        if key.startswith('NM_'):
            for transcript in transcripts_list:
                exons.extend(transcript['exons'])

    return [
        {
            'gene': gene,
            'chrom': chrom,
            'start': exon['start'],
            'end': exon['end'],
            'number': f"Exon {exon['exon_number']}",
            'size': exon['end'] - exon['start'] + 1
        }
        for exon in exons
    ]

# Creates a DataFrame to accomodate BED filtering
def create_bed(genome, gene, select_exons=None, trim_utr=False, padding=0):
    """
    Create a BED format DataFrame for selected exons of a given gene.

    Args:
        genome (str): Genome identifier (e.g., "hg38").
        gene (str): Gene symbol to generate BED regions for.
        select_exons (list[str], optional): Specific exon numbers to include (e.g., ["Exon 1", "Exon 2"]). If None or matches all,
        includes all exons. Default is None.
        trim_utr (bool, optional): Whether to trim 5' and 3' UTRs to coding regions. Default is False.
        padding (int, optional): Number of bases to add upstream and downstream. Default is 0.

    Returns:
        pd.DataFrame: BED-like DataFrame with columns:
            - CHROM (str): Chromosome
            - START (int): Start coordinate
            - END (int): End coordinate
            - GENE (str): Gene symbol
            - EXON (str): Exon number
            - SIZE (int): Exon size
    """

    shared_log.logger.info(f'Gene Selected: {gene}')

    shared_log.logger.info(f'Using triming: {trim_utr}; and Padding: {padding}')
    
    exons = get_exons(genome, gene)

    extracted = [exon['number'] for exon in exons]

    if set(select_exons) == set(extracted) or not select_exons:
        selected = [exon for exon in exons]
        utr_dict = utr_info(genome, gene)
        processed_exons = exon_number_json(selected, trim_utr, padding, utr_dict)
    else:
        selected = [exon for exon in exons if exon['number'] in select_exons]
        processed_exons = padding_adding(selected, padding)

    bed_lines = [
        {
            'CHROM': exon['chrom'],
            'START': exon['start']-1,
            'END': exon['end'],
            'GENE': gene,
            'EXON': exon['number'],
            'SIZE': exon['size']
        }
        for exon in processed_exons
    ]

    shared_log.logger.info('Creating BED DataFrame.')

    bed_df = pd.DataFrame(bed_lines)

    bed_df['chrom'] = bed_df['CHROM'].apply(extract_chr_value).astype(int)
    bed_df.sort_values(by=['chrom', 'START'], ascending=[True, True], inplace=True)

    bed_df['CHROM'] = bed_df['chrom'].apply(convert_numeric_to_chr).astype(str)

    bed_df.drop(columns=['chrom'], inplace=True)
    print(bed_df)

    return bed_df