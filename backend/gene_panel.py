import json
import pandas as pd
from pathlib import Path
from backend.helper import extract_chr_value, convert_numeric_to_chr, load_json, exon_number_json, utr_info, padding_adding
from log import shared_log

data_dir = Path('data/gene_panels')

# Load Gene Panel by user selection
def load_gene_panel(name):
    """
    Load a gene panel JSON file by its name.

    Args:
        name (str): Name of the gene panel (without `.json` extension).

    Returns:
        dict: Dictionary with panel information including genes and metadata.
    """

    filename = data_dir / f'{name}.json'

    with open(filename) as file:
        return json.load(file)

# Extract all the gene panels available
def get_gene_panel():

    panels = []

    for file in data_dir.glob("*.json"):
        if not file.stem.startswith("hg_"):
            panels.append(file.stem)

    return panels

# Extract all genes in gene panel
def get_panel_genes(name):
    """
    Retrieve all available gene panel names.

    Returns:
        list[str]: List of gene panel names.
    """

    shared_log.logger.info("Getting Gene Data For Each Gene In The Gene Panel.")

    data = load_gene_panel(name)

    return data.get('genes')

# Extract all the exons info for each gene in the gene panel
def get_panel_exons(genome, select_genes):
    """
    Extract the list of genes from a specific gene panel.

    Args:
        name (str): Name of the gene panel.

    Returns:
        list[str] | None: List of gene names in the panel, or None if not found.
    """

    shared_log.logger.info("Getting Exon Data For Each Gene In The Gene Panel.")

    data = load_json(genome)

    all_exons = []

    for gene in select_genes:
        gene_data = data.get(gene)
        if not gene_data:
            continue

        chrom = gene_data['chromosome']

        exons = []

        for key, value in gene_data.items():
            if key.startswith('NM_'):
                for transcript in value:
                    exons.extend(transcript['exons'])

        for exon in exons:
            all_exons.append({
                'chrom': chrom,
                'start': exon['start'],
                'end': exon['end'],
                'gene': gene,
                'number': f'Exon {exon['exon_number']}',
                'size': exon['end'] - exon['start']
            })

    return all_exons

# Creates a DataFrame to accomodate BED filtering
def create_panel_bed(genome, select_genes=None, select_exons=None, panel_name=None, trim_utr=False, padding=0):
    """
    Extract exon information for selected genes from a genome reference.

    Args:
        genome (str): Identifier for genome reference JSON.
        select_genes (list[str]): List of gene names to extract exons from.

    Returns:
        list[dict]: List of exon dictionaries, each containing:
            - chrom (str): Chromosome identifier.
            - start (int): Exon start coordinate.
            - end (int): Exon end coordinate.
            - gene (str): Gene name.
            - number (str): Exon number (e.g., 'Exon 1').
            - size (int): Exon length in base pairs.
    """

    if select_exons is None:
        select_exons = []

    if not panel_name:
        shared_log.logger.error("Required panel name.")

    shared_log.logger.info(f'Using Gene Panel: {panel_name}')

    if not select_genes:
        select_genes = get_panel_genes(panel_name)
        
    shared_log.logger.info(f'Selected genes: {select_genes}')

    shared_log.logger.info(f'Using triming: {trim_utr}; and Padding: {padding}')

    all_exons = get_panel_exons(genome, select_genes)

    extracted = [exon['number'] for exon in all_exons]

    if set(select_exons) == set(extracted) or not select_exons:
        select_exons = all_exons
        utr_dict = utr_info(genome, select_genes)
        processed_exons = exon_number_json(select_exons, trim_utr, padding, utr_dict)        
    else:
        select_exons = [exon for exon in all_exons if exon['number'] in select_exons]
        processed_exons = padding_adding(select_exons, padding)        
    
    bed_lines = [
        {
            'CHROM': exon['chrom'],
            'START': exon['start']-1,
            'END': exon['end'],
            'GENE': exon['gene'],
            'EXON': exon['number'],
            'SIZE': exon['size']
        }
        for exon in processed_exons
    ]
    
    shared_log.logger.info('Creating BED dataframe.')

    bed_df = pd.DataFrame(bed_lines)

    bed_df['chrom'] = bed_df['CHROM'].apply(extract_chr_value).astype(int)
    bed_df.sort_values(by=['chrom', 'START'], ascending=[True, True], inplace=True)

    bed_df['CHROM'] = bed_df['chrom'].apply(convert_numeric_to_chr).astype(str)

    bed_df.drop(columns=['chrom'], inplace=True)
    print(bed_df)

    return bed_df