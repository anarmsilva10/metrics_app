import re
import shutil
import json
from pathlib import Path
from datetime import datetime
from backend.retrieve_data import build_genome
from log import shared_log

data_dir = Path('data')

# Update JSON files
def update_file():
    """
    Update local genome JSON files if outdated.

    Checks for JSON files matching the pattern 'hg*.json' inside the `data` directory.
    If none exist or they are older than 7 days, removes outdated files and regenerates new genome files using `build_genome()`.

    Returns:
        None
    """

    files = list(data_dir.glob("hg*.json"))

    if not files:
        shared_log.logger.warning("No matching file found. Generating new file.")
        build_genome()
        return
    
    today_date = datetime.today().date()
    outdated_files = []
    needs_update = False
    
    for file in files:
        fname = file.stem

        try:
            name, date = fname.rsplit('_', 1)
            file_date = datetime.strptime(date, "%Y-%m-%d").date()
        except ValueError:
            shared_log.logger.warning(f'Invalid filename format: {fname}. Rebulding file.')
            outdated_files.append(file)
            needs_update = True
            continue
        
        days_old = (today_date - file_date).days
        if days_old >= 7:
            shared_log.logger.info(f'{fname} is {days_old} days old. Updating file.')
            outdated_files.append(file)
            needs_update = True

    if needs_update:
        shared_log.logger.info('Updating genome files.')

        for file in outdated_files:
            if file.exists():
                shared_log.logger.info(f"Deleting old file:{file.name}.")
                file.unlink()

        build_genome()
    else:
        shared_log.logger.info(f"All files are already updated.")

# Load JSON file choosed by the user
def load_json(genome):
    """
    Load genome data from a JSON file.

    Args:
        genome (str): Genome identifier (e.g., 'hg38').

    Returns:
        dict: Parsed JSON data containing genome annotation.
    """

    files = list(data_dir.glob(f"{genome}_*.json"))

    filename = files[0]

    with open(filename) as file:
        return json.load(file)

# Transforming chr X and Y into numeric to allow sorting
def extract_chr_value(chr_str):
    """
    Convert chromosome string into a numeric value for sorting.

    Args:
        chr_str (str): Chromosome string (e.g., 'chr1', 'X', 'Y').

    Returns:
        int: Numeric chromosome value.
            - from 1 to 22: Autosomes
            - 23: X
            - 24: Y
            - 25: Any non-standard chromosome
    """

    chr_str = str(chr_str).strip()

    match = re.match(r'^(chr)?(\d+|Y|X)', chr_str, re.IGNORECASE)

    if match:
        value = match.group(2).upper()

        if value.isdigit():
            return value
        elif value == 'X':
            return 23
        elif value == 'Y':
            return 24
        
    return 25

# Restablish chr X and Y after sorting
def convert_numeric_to_chr(value):
    """
    Convert numeric chromosome value back to string representation.

    Args:
        value (int | str): Numeric chromosome value or string.

    Returns:
        str: Chromosome string (from '1' to '22', 'X', 'Y').
    """

    if value == 23:
        return 'X'
    elif value == 24:
        return 'Y'
    elif isinstance(value, int):
        return f'{value}'
    
    return value

# Deleting temporary folders
def clean_temp_dirs():
    """
    Delete all temporary directories inside the data directory.
    Removes directories matching the pattern 'tmp*'.

    Raises:
        Exception: If directory removal fails.
    """

    for temp_folder in data_dir.glob('tmp*'):
        if temp_folder.is_dir():
            try:
                shutil.rmtree(temp_folder)
            except Exception:
                raise

# Getting UTR information
def utr_info(genome, selected_genes):
    """
    Extract UTR (untranslated region) information for selected genes.

    Args:
        genome (str): Genome identifier.
        selected_genes (str | list[str]): Gene name or list of genes.

    Returns:
        dict: Dictionary mapping each gene to UTR coordinates with keys:
            - '5_start': 5' UTR start
            - '5_end': 5' UTR end
            - '3_start': 3' UTR start
            - '3_end': 3' UTR end
    """

    data = load_json(genome)
    utr_dict = {}

    # Needed to garantee that when just one gene is selected that comes as a list
    if isinstance(selected_genes, str):
        selected_genes = [selected_genes]
    
    for gene in selected_genes:
        gene_data = data.get(gene)

        if not gene_data:
            shared_log.logger.warning(f'No data found for gene: {gene}')
            continue

        utr_dict[gene] = {
            "5_start": gene_data.get("5' UTR start"),
            "5_end": gene_data.get("5' UTR end"),
            "3_start": gene_data.get("3' UTR start"),
            "3_end": gene_data.get("3' UTR end")
        }

    return utr_dict

# Adding Padding
def padding_adding(exons, padding):
    """
    Apply padding to exon regions.

    Args:
        exons (list[dict]): List of exon dictionaries with keys ['chrom', 'start', 'end', 'number', 'gene'].
        padding (int): Number of bases to add upstream and downstream.

    Returns:
        list[dict]: List of updated exon dictionaries with padding applied, including recalculated size.
    """

    processed_exons = []

    for exon in exons:
        exon_start = exon['start']
        exon_end = exon['end']

        # Apply padding
        exon_start = max(1, exon_start - padding)
        exon_end = exon_end + padding

        processed_exons.append({
            "chrom": exon['chrom'],
            "start": exon_start,
            "end": exon_end,
            "number": exon['number'],
            "size": exon_end - exon_start + 1,
            "gene": exon['gene']
        })

    return processed_exons

# Trimming and Padding 
def exon_number_json(exons, trim_utr, padding, utr_dict):
    """
    Numbering of exons in order of their appearance in the genome.
    Optionally trims 5' and 3' UTRs to include only coding regions,
    and applies padding to include splice site regions.

    Args:
        exons (list[dict]): List of exon dictionaries with keys ['chrom', 'start', 'end', 'gene'].
        trim_utr (bool): Whether to trim 5' and 3' UTRs.
        padding (int): Number of bases to add upstream and downstream.
        utr_dict (dict): Dictionary with UTR coordinates for genes.

    Returns:
        list[dict]: List of processed exon dictionaries with updated coordinates, padding applied, UTR trimming applied 
        if requested, and sequential exon numbering per gene.
    """

    exon_counter = {}
    processed_exons = []

    for exon in exons:
        exon_start = exon['start']
        exon_end = exon['end']
        gene = exon['gene']

        if trim_utr:
            utr = utr_dict.get(gene, {})
            utr_5_start = utr.get("5_start")
            utr_5_end = utr.get("5_end")
            utr_3_start = utr.get("3_start")
            utr_3_end = utr.get("3_end")

            # Trim 5' UTR
            if utr_5_start is not None and utr_5_end is not None:
                exon_start = max(exon_start, utr_5_end + 1)

            # Trim 3' UTR
            if utr_3_start is not None and utr_3_end is not None:
                exon_end = min(exon_end, utr_3_start - 1)

            if exon_start > exon_end:
                continue

        # Apply padding
        exon_start = max(1, exon_start - padding)
        exon_end = exon_end + padding

        exon_counter[gene] = exon_counter.get(gene, 0) + 1
        exon_number = exon_counter[gene]

        processed_exons.append({
            "chrom": exon['chrom'],
            "start": exon_start,
            "end": exon_end,
            "number": f"Exon {exon_number}",
            "size": exon_end - exon_start + 1,
            "gene": exon['gene']
        })

    return processed_exons