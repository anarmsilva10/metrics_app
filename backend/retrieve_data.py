import requests
import pandas as pd
from io import StringIO
import re
import json
import os
import argparse
from datetime import datetime
from log import shared_log

# BioMart service URL
BIOMART_URL_HG38 = "https://www.ensembl.org/biomart/martservice"

BIOMART_URL_HG19 = "https://grch37.ensembl.org/biomart/martservice"

#BioMart Query in xml
BIOMART_SELECT_REFSEQ = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
           
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
        <Filter name = "mane_select" excluded = "0"/>
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "hgnc_symbol" />
        <Attribute name = "hgnc_id" />
        <Attribute name = "transcript_mane_select" />
        <Attribute name = "transcript_mane_plus_clinical" />
    </Dataset>
</Query>"""

BIOMART_PC_REFSEQ = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
           
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
        <Filter name = "mane_plus_clinical" excluded = "0"/>
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "hgnc_symbol" />
        <Attribute name = "hgnc_id" />
        <Attribute name = "transcript_mane_select" />
        <Attribute name = "transcript_mane_plus_clinical" />
    </Dataset>
</Query>"""

BIOMART_HG19_REFSEQ = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
           
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "hgnc_symbol" />
        <Attribute name = "hgnc_id" />
        <Attribute name = "refseq_mrna" />
    </Dataset>
</Query>"""

BIOMART_MANE_SELECT = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "mane_select" excluded = "0"/>
        <Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "chromosome_name" />
        <Attribute name = "exon_chrom_start" />
        <Attribute name = "exon_chrom_end" />
        <Attribute name = "5_utr_start" />
        <Attribute name = "5_utr_end" />
        <Attribute name = "3_utr_start" />
        <Attribute name = "3_utr_end" />
        <Attribute name = "strand" />
    </Dataset>
</Query>"""

BIOMART_PLUS = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
        <Filter name = "mane_plus_clinical" excluded = "0" />
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "chromosome_name" />
        <Attribute name = "exon_chrom_start" />
        <Attribute name = "exon_chrom_end" />
        <Attribute name = "5_utr_start" />
        <Attribute name = "5_utr_end" />
        <Attribute name = "3_utr_start" />
        <Attribute name = "3_utr_end" />
        <Attribute name = "strand" />
    </Dataset>
</Query>"""

BIOMART_HG19 = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "chromosome_name" value = "CHROM_RANGE"/>
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "chromosome_name" />
        <Attribute name = "exon_chrom_start" />
        <Attribute name = "exon_chrom_end" />
        <Attribute name = "5_utr_start" />
        <Attribute name = "5_utr_end" />
        <Attribute name = "3_utr_start" />
        <Attribute name = "3_utr_end" />
        <Attribute name = "strand" />
        <Attribute name = "external_gene_name" />
    </Dataset>
</Query>"""

# Data Retrieving

# Removing spaces from biomart query
def compact_biomart_query(xml_query: str) -> str:
    """
    Removes extra whitespace and tabs from a BioMart XML query to make it more compact.

    Args:
        xml_query (str): The raw XML BioMart query.

    Returns:
        str: The compacted XML query string.
    """
    return re.sub(r">\s+<", "><", xml_query.strip().replace("\t", ""))

# Retrieving RefSeq ID for Mane Select data
def refseq_select():
    """
    Retrieves RefSeq IDs for MANE Select transcripts from BioMart using xml.
   
    Args:
        None
   
    Returns:
        pd.DataFrame: DataFrame containing transcript IDs, gene symbols, HGNC IDs, and RefSeq MANE Select matches.
    """

    shared_log.logger.info("Retrieving RefSeq ID - Mane Select data.")

    query = compact_biomart_query(BIOMART_SELECT_REFSEQ)

    response = requests.get(BIOMART_URL_HG38, params={"query":query})

    if response.status_code != 200:
        shared_log.logger.error("Error retrieving RefSeq ID for Mane Select data.")

    # Convert TSV text response to DataFrame
    refseq_select_df = pd.read_csv(StringIO(response.text), sep="\t")

    refseq_select_df.rename(columns={'RefSeq match transcript (MANE Select)': 'mane_select'}, inplace=True)
    refseq_select_df.rename(columns={'RefSeq match transcript (MANE Plus Clinical)': 'plus_clinical'}, inplace=True)
   
    return refseq_select_df

# Retrieve MANE Select data
def mane_select():
    """
    Retrieves MANE Select exon data from BioMart using xml.

    Args:
        None

    Returns:
        pd.DataFrame: DataFrame containing chromosome, exon start, end, UTRs, strand, and transcript information for MANE Select.
    """

    shared_log.logger.info("Retrieving MANE Select Data from Biomart.")

    query = compact_biomart_query(BIOMART_MANE_SELECT)

    response = requests.get(BIOMART_URL_HG38, params={"query": query})

    if response.status_code != 200:
        shared_log.logger.error("Error retrieving MANE Select data!")
        return None
   
    # Convert TSV text response to DataFrame
    mane_select_df = pd.read_csv(StringIO(response.text), sep="\t", low_memory=False)

    return mane_select_df

# Retrieving RefSeq ID for Mane Plus Clinical data
def refseq_pc():
    """
    Retrieves RefSeq IDs for MANE Plus Clinical transcripts from BioMart using xml.

    Args:
        None

    Returns:
        pd.DataFrame: DataFrame containing transcript IDs, gene symbols, HGNC IDs, and RefSeq MANE Plus Clinical matches.
    """

    shared_log.logger.info("Retrieving RefSeq ID - Mane Plus Clinical data.")

    query = compact_biomart_query(BIOMART_PC_REFSEQ)

    response = requests.get(BIOMART_URL_HG38, params={"query":query})

    if response.status_code != 200:
        shared_log.logger.error("Error retrieving RefSeq ID for Mane Plus Clinical data.")

    # Convert TSV text response to DataFrame
    refseq_pc_df = pd.read_csv(StringIO(response.text), sep="\t")

    refseq_pc_df.rename(columns={'RefSeq match transcript (MANE Select)': 'mane_select'}, inplace=True)
    refseq_pc_df.rename(columns={'RefSeq match transcript (MANE Plus Clinical)': 'plus_clinical'}, inplace=True)

    return refseq_pc_df

# Retrieve MANE Plus Clinical data
def mane_plus_clinical():
    """
    Retrieves MANE Plus Clinical exon data from BioMart using xml.

    Args:
        None

    Returns:
        pd.DataFrame: DataFrame containing chromosome, exon start, end, UTRs, strand, and transcript information from MANE Plus Clinical.
    """

    shared_log.logger.info("Retrieving MANE Plus Clinical Data from Biomart.")

    query = compact_biomart_query(BIOMART_PLUS)

    response = requests.get(BIOMART_URL_HG38, params={"query": query})

    if response.status_code != 200:
        shared_log.logger.error("Error retrieving MANE Plus Clinical data!")

    # Convert TSV text response to DataFrame
    mane_pc_df = pd.read_csv(StringIO(response.text), sep="\t")

    return mane_pc_df

# Retrieving RefSeq ID for hg19 data
def hg19_refseq():
    """
    Retrieves RefSeq IDs for hg19 transcripts from BioMart using xml.

    Args:
        None

    Returns:
        pd.DataFrame: DataFrame with transcript IDs, gene symbols, HGNC IDs, and RefSeq IDs (renamed from “RefSeq mRNA ID”).
    """

    shared_log.logger.info("Retrieving RefSeq ID - hg19 data.")

    query = compact_biomart_query(BIOMART_HG19_REFSEQ)

    response = requests.get(BIOMART_URL_HG19, params={"query":query})

    if response.status_code != 200:
        shared_log.logger.error("Error retrieving RefSeq ID for hg19 data.")

    # Convert TSV text response to DataFrame
    refseq_hg19_df = pd.read_csv(StringIO(response.text), sep="\t")

    refseq_hg19_df.rename(columns={'RefSeq mRNA ID': 'RefSeq ID'}, inplace=True)

    return refseq_hg19_df

# Retrieve hg19 data
def hg19():
    """
    Retrieving hg19 data from Biomart using xml in chunks to avoid timeouts or truncation.

    Args:
        None
   
    Returns:
        pd.DataFrame: Concatenated DataFrame with all chromosomes.
        None: If all requests fail.
    """

    chromosomes = list(map(str, range(1,23))) + ['X', 'Y']
    chunk_size = 5

    dfs = []
    for i in range(0, len(chromosomes), chunk_size):
        chrom_chunk = chromosomes[i:i + chunk_size]
        chrom_str = ",".join(chrom_chunk)

        # Change CHROM_RANGE to chromosome range
        query = BIOMART_HG19.replace("CHROM_RANGE", chrom_str)
        query = compact_biomart_query(query)
       
        response = requests.get(BIOMART_URL_HG19, params={"query": query})

        if response.status_code != 200:
            shared_log.logger.error("Error retrieving hg19 data!")
            return None
       
        df_chunk = pd.read_csv(StringIO(response.text), sep="\t", low_memory=False)
        dfs.append(df_chunk)

    if dfs:
        shared_log.logger.info("Concatenating hg19 dataframes")
        return pd.concat(dfs, ignore_index=True)
    else:
        shared_log.logger.error("No data retrieved from BioMart.")
        return pd.DataFrame()

# Data Processing

# Sort DataFrame by chromossome and start position
def sort(df):
    """
    Sorts a BioMart DataFrame by Transcript ID and exon start position.

    Args:
        df (pd.DataFrame): DataFrame containing exon and transcript data.

    Returns:
        pd.DataFrame: Sorted DataFrame.
        None: If input DataFrame is empty or None.
        """

    if df is None or df.empty:
        shared_log.logger.warning("Warning: Dataframe is empty!")
        return None

    # Sorted by Transcript ID because one gene can have more than one transcript, this way it garantees that every exon of each
    # transcript it sorted together; and then inside the trancript, its sorted by exon start
    sorted_df = df.sort_values(
            by=["Transcript stable ID", "Exon region start (bp)"],
            ascending=[True, True]
        ).reset_index(drop=True)

    return sorted_df

# Assing exon number
def exon_number(df, trim_utr, padding):
    """
    Numbering of exons in order of their appearance in the genome.
    Optionally trims 5' and 3' UTRs to include only coding regions,
    and applies padding to include splice site regions.
    Note: there is no need to look at the strand because the positions are swaped from biomart.

    Args:
        df (pd.DataFrame): BioMart-style DataFrame with exon/transcript data.
        trim_utr (bool): If True, trims UTRs to retain only coding regions.
        padding (int): Number of base pairs to extend exon boundaries (default = 10).

    Returns:
        pd.DataFrame: DataFrame with 'Exon number' and adjusted coordinates.
    """

    current_transcript = None
    exon_counter = 0
    processed_rows = []

    for idx, row in df.iterrows():
        transcript = row['Transcript stable ID']
        exon_start = row['Exon region start (bp)']
        exon_end = row['Exon region end (bp)']

        if trim_utr:
            # Trim 5' UTR
            if not pd.isna(row["5' UTR start"]) and not pd.isna(row["5' UTR end"]):
                exon_start = max(exon_start, row["5' UTR end"] + 1)

            # Trim 3' UTR
            if not pd.isna(row["3' UTR start"]) and not pd.isna(row["3' UTR end"]):
                exon_end = min(exon_end, row["3' UTR start"] - 1)

            # Skip if trimming removes the whole exon
            if exon_start > exon_end:
                continue

        # Apply padding
        exon_start = max(1, exon_start - padding)  # avoid positions < 1
        exon_end = exon_end + padding

        # Update row
        row['Exon region start (bp)'] = exon_start
        row['Exon region end (bp)'] = exon_end

        # Number the exon
        if transcript != current_transcript:
            current_transcript = transcript
            exon_counter = 1
        else:
            exon_counter += 1

        row['Exon number'] = exon_counter
        processed_rows.append(row)

    return pd.DataFrame(processed_rows)

def choose_mane(row):
    """
    Selects the MANE Plus Clinical transcript ID if available, otherwise falls back to MANE Select.

    Args:
        row (pd.Series): Row from a DataFrame containing plus_clinical and mane_select columns.

    Returns:
        str | None: Preferred transcript ID.
    """

    if pd.notna(row['plus_clinical']):
        return row['plus_clinical']
    else:
        return row['mane_select']

# Merge RefSeq data from Mane Select and Mane Plus Clinical
def merge_refseq(refseq_select_df, refseq_pc_df):
    """
    Merges RefSeq data from MANE Select and MANE Plus Clinical sources.

    Args:
        refseq_select_df (pd.DataFrame): DataFrame from MANE Select.
        refseq_pc_df (pd.DataFrame): DataFrame from MANE Plus Clinical.

    Returns:
        pd.DataFrame: Combined RefSeq dataset with unique transcript IDs and cleaned RefSeq IDs.
    """

    refseq_df = pd.merge(refseq_select_df[['HGNC symbol', 'Transcript stable ID', 'mane_select']],
                         refseq_pc_df[['HGNC symbol', 'Transcript stable ID', 'plus_clinical']],
                         on=['HGNC symbol', 'Transcript stable ID'], how='outer')

    refseq_df['RefSeq ID'] = refseq_df.apply(lambda row: choose_mane(row), axis=1)
    refseq_df['RefSeq ID'] = [t.rsplit('.', 1)[0] for t in refseq_df['RefSeq ID']]

    refseq_df = refseq_df.groupby('Transcript stable ID', as_index=False).first()

    refseq_df = refseq_df.drop(columns={'mane_select', 'plus_clinical'})

    return refseq_df

# Merge MANE data by chromossome, start and end position
def merge_mane(mane_select_df, mane_pc_df, trim_utr, padding):
    """
    Merges MANE Select and MANE Plus Clinical exon datasets by chromosome, start, and end position.

    Args:
        mane_select_df (pd.DataFrame): DataFrame from MANE Select.
        mane_pc_df (pd.DataFrame): DataFrame from MANE Plus Clinical.
        trim_utr (bool): If True, trims UTRs.
        padding (int): Number of base pairs to extend exon boundaries.

    Returns:
        pd.DataFrame: Combined MANE dataset with sorted exons and exon numbers.
        None: If input datasets are missing.
    """
   
    if mane_select_df is None or mane_pc_df is None:
        shared_log.logger.warning("Missing data for merge.")
        return None

    # Define keys for identifying unique exons
    merge_keys = ['Chromosome/scaffold name', 'Exon region start (bp)', 'Exon region end (bp)']

    # Which MANE Plus Clinical rows that are NOT in MANE Select
    pc_df = mane_pc_df.merge(
        mane_select_df[merge_keys],
        on=merge_keys,
        how='left',
        indicator=True
    ).query('_merge == "left_only"').drop(columns=['_merge'])

    # Combine MANE Select with exclusive MANE Plus
    mane_df = pd.concat([mane_select_df, pc_df], ignore_index=True)
    mane_df = mane_df.drop_duplicates(subset=['Transcript stable ID', 'Chromosome/scaffold name',
                                              'Exon region start (bp)', 'Exon region end (bp)'])

    shared_log.logger.info("Sorting positions and chromossomes in MANE dataframe.")
    mane_df = sort(mane_df)

    shared_log.logger.info("Numbering exons in MANE dataframe")
    mane_df = exon_number(mane_df, trim_utr, padding)

    return mane_df

def hg38(refseq_select_df, refseq_pc_df, mane_select_df, mane_pc_df, trim_utr, padding):
    """
    Builds a complete hg38 dataset by merging RefSeq and MANE exon data.

    Args:
        refseq_select_df (pd.DataFrame): RefSeq from MANE Select.
        refseq_pc_df (pd.DataFrame): RefSeq from MANE Plus Clinical.
        mane_select_df (pd.DataFrame): MANE Select exon data.
        mane_pc_df (pd.DataFrame): MANE Plus Clinical exon data.
        trim_utr (bool): If True, trims UTRs.
        padding (int): Exon boundary padding.

    Returns:
        pd.DataFrame: Merged hg38 dataset with exon data.
    """

    # Merge RefSeq data
    refseq_df = merge_refseq(refseq_select_df, refseq_pc_df)

    # Merge MANE exon data
    mane_df = merge_mane(mane_select_df, mane_pc_df, trim_utr, padding)    

    # Final DataFrame with all info about hg38 data from MANE
    hg38_df = pd.merge(mane_df, refseq_df, on='Transcript stable ID', how='inner')

    return hg38_df

# HG19 filtering by MANE data
def hg19_filter(refseq_hg19_df, hg19_df, hg38_df, trim_utr, padding):
    """
    Filters hg19 dataset to include only transcripts that match hg38 (MANE) RefSeq IDs.

    Args:
        refseq_hg19_df (pd.DataFrame): hg19 RefSeq IDs.
        hg19_df (pd.DataFrame): hg19 exon dataset.
        hg38_df (pd.DataFrame): hg38 MANE dataset.
        trim_utr (bool): If True, trims UTRs.
        padding (int): Exon boundary padding.

    Returns:
        pd.DataFrame: Filtered hg19 dataset with sorted and numbered exons.
        None: If input datasets are missing.
        """

    if hg19_df is None or hg38_df is None:
        shared_log.logger.warning("Missing data for filtering.")
        return None

    # Merging RefSeq and exon data from hg19 data
    hg19_merged_df = pd.merge(refseq_hg19_df, hg19_df, on='Transcript stable ID', how='inner')

    # Filter hg19 to keep only rows with transcript IDs present in hg38 (MANE data)
    hg19_filtered_df = hg19_merged_df[hg19_merged_df['RefSeq ID'].isin(hg38_df['RefSeq ID'])].copy()

    shared_log.logger.info("Sorting positions and chromossomes in hg19 dataframe.")
    hg19_filtered_df = sort(hg19_filtered_df)

    shared_log.logger.info("Numbering exons in hg19 dataframe.")
    hg19_filtered_df = exon_number(hg19_filtered_df, trim_utr, padding)

    return hg19_filtered_df

# Data Transformation

# Convert DataFrame to JSON
def df_to_json(df):
    """
    Converts a transcript DataFrame into a structured JSON-like dictionary.

    Args:
        df (pd.DataFrame): DataFrame with exon and transcript data.

    Returns:
        dict: JSON-style dictionary grouped by gene name and RefSeq ID.
    """

    result = {}

    for gene_name, group in df.groupby('HGNC symbol'):
        group = group.sort_values('Exon region start (bp)')

        chromosome = str(group['Chromosome/scaffold name'].iloc[0])
        strand = int(group['Strand'].iloc[0])

        # Dealing with NA
        def_convert = lambda x: int(x) if pd.notna(x) else None

        utr5_start = def_convert(group["5' UTR start"].iloc[0])
        utr5_end   = def_convert(group["5' UTR end"].iloc[0])
        utr3_start = def_convert(group["3' UTR start"].iloc[0])
        utr3_end   = def_convert(group["3' UTR end"].iloc[0])

        refseq_id = group['RefSeq ID'].iloc[0]
        group = group[group['RefSeq ID'] == refseq_id]

        exons = []

        for _, row in group.iterrows():
            exons.append({
                    "exon_number": int(row['Exon number']),
                    "start": int(row['Exon region start (bp)']),
                    "end": int(row['Exon region end (bp)'])
                })
           
        # Creating an entry for each gene
        if gene_name not in result:
            result[gene_name] = {
                "chromosome": chromosome,
                "strand": strand,
                "5' UTR start": utr5_start,
                "5' UTR end": utr5_end,
                "3' UTR start": utr3_start,
                "3' UTR end": utr3_end,
                refseq_id: [
                    {"exons": exons}
                ]
            }

    return result

# Generate json
def generate_json(hg38_df, hg19_df, path="data"):
    """
    Converts hg38 and hg19 DataFrames into JSON and writes them to disk.

    Args:
        hg38_df (pd.DataFrame): hg38 dataset.
        hg19_df (pd.DataFrame): hg19 dataset.
        path (str, optional): Directory where JSON files will be saved (default = "data").

    Returns:
        None
    """

    os.makedirs(path, exist_ok=True)

    hg38_json = df_to_json(hg38_df)
    hg19_json = df_to_json(hg19_df)

    with open(os.path.join(path, f"hg38_{datetime.today().strftime('%Y-%m-%d')}.json"), "w") as file:
        json.dump(hg38_json, file, indent=2)

    with open(os.path.join(path, f"hg19_{datetime.today().strftime('%Y-%m-%d')}.json"), "w") as file:
        json.dump(hg19_json, file, indent=2)

    shared_log.logger.info("Exported JSON files for hg38 and hg19.")

    return

def build_genome(trim_utr=False, padding=0):
    """
    Runs the full genome-building pipeline: retrieves MANE/hg19 data, merges them, filters hg19 with hg38, and generates JSON outputs.

    Args:
        trim_utr (bool, optional): Whether to trim UTRs. Default = False.
        padding (int, optional): Padding to apply to exon boundaries. Default = 0.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: (hg38 dataset, filtered hg19 dataset).
    """

    # Mane Select DataFrame
    shared_log.logger.info("Creating Mane Select DataFrame.")
    refseq_select_df = refseq_select()
    mane_select_df = mane_select()

    # Mane Plus Clinical DataFrame
    shared_log.logger.info("Creating Mane Plus Clinical DataFrame.")
    refseq_pc_df = refseq_pc()
    mane_pc_df = mane_plus_clinical()

    # hg38 DataFrame
    shared_log.logger.info("Creating hg38 DataFrame.")
    hg38_df = hg38(refseq_select_df, refseq_pc_df, mane_select_df, mane_pc_df, trim_utr, padding)

    # hg19 DataFrame
    shared_log.logger.info("Creating hg19 DataFrame.")
    refseq_hg19_df = hg19_refseq()
    hg19_df = hg19()
    shared_log.logger.info("Filtering hg19 with hg38.")
    hg19_filtered_df = hg19_filter(refseq_hg19_df, hg19_df, hg38_df, trim_utr, padding)

    # Transforming DataFrames into Json to extract data
    generate_json(hg38_df, hg19_filtered_df)
   
    return hg38_df, hg19_filtered_df

def main ():

    parser = argparse.ArgumentParser(description='Build genome data files.')
    parser.add_argument('--genome', type=str, choices=['hg19', 'hg38'])
    parser.add_argument('--trim_utr', action='store_true', help="Add triming for 5' and 3' UTRs to include only coding regions")
    parser.add_argument('--padding', type=int, default=0, help="Applies padding to include splice site regions")
    args = parser.parse_args()

    build_genome(trim_utr=args.trim_utr, padding=args.padding)

# Calling main function
if __name__ == "__main__":
    main()