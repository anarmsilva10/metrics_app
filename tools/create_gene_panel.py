import argparse
import requests
import json
import os
import re
from log import shared_log


def fetch_gene_panels(ids_file: str, output_dir: str):
    # Load test IDs from file
    with open(ids_file, "r") as file:
        test_ids = [line.strip() for line in file.readlines() if line.strip()]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for test_id in test_ids:
        gtr_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gtr&id={test_id}&retmode=json"

        response = requests.get(gtr_url)

        if response.status_code == 200:
            data = response.json()

            try:
                panel_data = data['result'][test_id]
                panel_name = panel_data.get('testname')

                filename = re.sub(r'[^\w\-_\.]', '_', panel_name) + '.json'
                genes = []

                if 'analytes' in panel_data:
                    for analyte in panel_data.get('analytes', []):
                        gene_name = analyte.get('name')
                        if gene_name:
                            genes.append(gene_name)

                if genes:
                    json_file = os.path.join(output_dir, filename)
                    with open(json_file, "w") as file:
                        json.dump({'genes': genes}, file, indent=4)

                    shared_log.logger.info(f"Retrieved {len(genes)} genes for Gene Panel: {filename}")
                else:
                    shared_log.logger.error(f"No genes found for Gene Panel: {filename}.")

            except Exception as e:
                shared_log.logger.error(f"Failed to parse data for test ID: {test_id}. Error: {e}")

        else:
            shared_log.logger.error(f"HTTP error: {response.status_code} for test ID: {test_id}")


def main():
    parser = argparse.ArgumentParser(description="Fetch gene panels from GTR by test IDs.")
    parser.add_argument(
        "--ids-file",
        default="tools/data/gene_panels_ids.txt",
        help="Path to file containing GTR test IDs (one per line).",
    )
    parser.add_argument(
        "--output-dir",
        default="data/gene_panels",
        help="Directory to save gene panel JSON files.",
    )

    args = parser.parse_args()
    fetch_gene_panels(ids_file=args.ids_file, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
