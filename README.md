# Metrics App

## Description
Metrics App is a bioinformatics application developed in Python and Shiny to perform genetic analyses, with a focus on quality management of BAM/CRAM files. It supports the analysis of **somatic cases**, **individual genes**, and **gene panels** using BED files. The software calculates two main types of metrics:
 - **Statistical Metrics** – including Reads Mapped, Duplicated Reads, and Mapping Quality, which are computed using pysam.samtools.
 - **Coverage Metrics** – such as Breadth of Coverage and Depth of Coverage, calculated using three complementary methods that rely on pysam.coverage, pysam.depth or samtools depth.

To facilitate interpretation, the app provides graphical visualizations of coverage depth. Users can also download a comprehensive report summarizing all calculated metrics.

Additionally, the backend scripts responsible for metric computation can be executed independently from the command line, allowing integration into custom bioinformatics pipelines, with the command:
```python backend/coverage_metrics.py <cram_path> <bed_path> <output_path> --method <method> --stats --gene_selection <GENE_NAME> --exon_selection <EXON_NUMBER> --output_plot <output_plot_path>```

For more information, please check `metrics_documentation.ipynb`.

## Project Architecture
```
metrics_app
├── README.md
├── Dockerfile
├── environment.yml
├── app.py
├── backend
│   ├── coverage_metrics.py
│   ├── depth.py
│   ├── plot.py
│   ├── samtool_stats.py
│   ├── gene_panel.py
|   ├── helper.py
│   ├── retrieve_data.py
|   ├── s3.py
│   ├── single_gene.py
│   └── aws_config
|       ├── s3_buckets.py
│       └── s3_buckets.yaml
├── log
│   └── shared_log.py
├── data
|   ├── gene_panels
│   └── reference_genomes
└── tools
    ├── create_gene_panel.py
    └── create_gene_panel.ipynb
```

## Features
 - **Somatic Case Analysis** – Evaluate global sequencing quality across samples, including mapping, duplication, and coverage metrics.
 - **Single Gene Analysis** – Assess quality metrics for one gene, with detailed exon-level coverage visualization.
 - **Gene Panel Analysis** – Analyze multiple genes from a BED file and summarize per-panel, per-gene and per-exon quality metrics.
 - **Custom BED Analysis** – Upload custom BED regions for flexible, targeted quality evaluation.
 - **Visualization** – Generate interactive plots for coverage depth and summary tables of computed metrics.
 - **Reporting** – Export comprehensive reports (HTML, PDF, or CSV) including metrics and visualizations.
 - **S3 Integration** – Access BAM/CRAM files directly from Amazon S3, namely the 1000 Genomes Project, without local downloads for Single Gene, Gene Panel and Custom BED analysis. Note: to change S3 Bucket, please, read the aws_config README file.
 - **Command-line Support** – Run backend Python scripts independently for automation and pipeline integration.

 ## Requisites
 - Python 3.12.2
 - Docker
 - Conda (optional) - dependencies management

 ### Python Deperndencies
 The file environment.yml have all dependencies needed for this project:
```yaml
name: quality_app_env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.12
  - pip
  - samtools
  - pip:
      - boto3
      - pandas
      - plotly
      - pysam
      - pyyaml
      - shiny
      - shinywidgets
      - shinyswatch
      - requests
      - libsass
      - openpyxl
      - dotenv
 ```
 
## Installation and Setup

In order for this application to work, the log, depth_calculator, and aws_config repositories are required.

### Docker
  1. Clone the repository:
  ```
    git clone https://github.com/anarmsilva10/metrics_app
    cd metrics_app
  ```

  2. Build Docker Image:
  ```
    docker build -t metrics_app .
  ```

  3. Run Docker Container:
  ```
    docker run -p 8000:8000 metrics_app
  ```

The application will be accessible at http://localhost:8000.

### Conda
  1. Clone the repository:
  ```
    git clone https://github.com/anarmsilva10/metrics_app
    cd metrics_app
  ```

  2. Create Conda environment:
  ```
    conda env create -f environment.yml
    conda activate quality_app_env
  ```

  3. Run the application:
  ```
    shiny run app.py
  ```

## Required Files
The Metrics App requires two main input to run correctly:

  1. Gene Panels (JSON format):
      - Gene Panels are required for Gene Panel Analysis.
      - Gene panels must contain the HGNC symbol of the genes that belong to them. In the tools folder, there is a Jupiter notebook (```create_gene_panel.ipynb```) that helps with the construction of gene panels.
      - Gene panels must be stored in ```data/gene_panels/```.

  2. Reference Genome:
      - The reference genome is required for mapping and coverage analysis.
      - It must be placed under: ```data/reference_genome/```
      - Before running the app, it is necessary to populate the sequence cache. This prepares the reference genome for use by samtools and pysam:
          1. Recover the humane reference genome:
              Humane Reference Genome for [hg19](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/); Humane Reference Genome for [hg38](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)
          2. Populate the sequence chace:
            ```seq_cache_populate.pl -root=data/reference_genome/ -subdirs 2 *.fa ```
      - [More information](https://www.htslib.org/doc/reference_seqs.html)

## Aplication Pages
 - Somatics Page: Designed for **somatic case analysis**, allowing users to evaluate sequencing quality.
 - WES Page: Performs **whole-exome sequencing (WES)** quality analysis.
    - Single Gene SubPage: Focused on detailed analysis of a **single gene**, displays gene and exon-level coverage depth and quality metrics for the selected gene, provides coverage plot to monitor the quality.
    - Gene Panel SubPage: Performs quality evaluation of multiple genes defined in a **gene panel**, summarizes per-panel, per-gene and per-exon metrics to assess overall panel performance, rovides coverage visualizations for each gene and exon.
    - Custom BED SubPage: Enables users to upload or select a **custom BED file** to analyze any genomic region of interest, ideal for targeted or experimental regions outside standard panels.
 - About Page: Displays **application information**, version details, and credits.

## Contribution:
To contribute to this project:
 1. Fork the repository.
 2. Create a new branch: ```git checkout -b my-branch```.
 3. Commit your changes: ```git commit -m 'changes'```.
 4. Push to your branch: ```git push origin my-branch```.
 5. Open a pull request.

## Contact
Developer: Ana Rita Silva
Email: anar.m.silva10@gmail.com

## License
This project is licensed under the terms of the MIT License.