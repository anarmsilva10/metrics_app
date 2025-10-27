import boto3
from pathlib import Path
import tempfile
import os
from log import shared_log
from backend.aws_config.s3_buckets import S3BucketConfig

# Load bucket configurations from YAML file
cfg = getattr(S3BucketConfig(), "1000genomes")

# Create an S3 client using credentials from the config
s3_client = boto3.client("s3", **cfg.config)

bucket = cfg.phase3

# Get list of available samples to display for user selection
def list_samples() -> list[str]:
    """
    Get a list of available samples from the S3 bucket.

    Returns:
        list[str]: List of sample identifiers (folder names) available under the `phase3/data/` prefix.
    """

    shared_log.logger.info("Extracting Available Samples")

    base_prefix = "phase3/data/"

    response = s3_client.list_objects_v2(
        Bucket=bucket.Bucket,
        Prefix=base_prefix,
        Delimiter="/"
    )

    if "CommonPrefixes" not in response:
        return []

    samples = [
        prefix["Prefix"].replace(base_prefix, "").strip("/")
        for prefix in response["CommonPrefixes"]
    ]
    return samples

# Replace sample with the choosed ID sample
def get_sample(sample: str) -> str:
    """
    Format the S3 key prefix for a given sample.

    Args:
        sample (str): Sample identifier to replace in the configured prefix.

    Returns:
        str: S3 key prefix with the selected sample inserted.
    """

    return bucket.Prefix.format(sample=sample)

# Function to parse the files for shiny app
def list_files(sample: str) -> list[str]:
    """
    List available files (.bam or .cram) for a given sample.

    Args:
        sample (str): Sample identifier.

    Returns:
        list[str]: List of file keys (paths in S3) for the selected sample.
    """

    shared_log.logger.info("Extracting Available Files")

    prefix = get_sample(sample)

    # List objects in the bucket
    response = s3_client.list_objects_v2(Bucket=bucket.Bucket, Prefix=prefix)

    if "Contents" not in response:
        return []
   
    return [obj["Key"] for obj in response.get("Contents", [])
            if obj["Key"].endswith(".bam") or obj["Key"].endswith(".cram")]

download_cache = {}

# File and Index File temporary download
def download_files(key: str) -> tuple[str, str, str]:
    """
    Download a file (.bam or .cram) and its corresponding index file
    (.bai or .crai) from S3 to a local temporary directory.

    Ensures that if files already exist in a previously created temporary
    directory, they are reused instead of being downloaded again.

    Args:
        key (str): S3 key of the main file (.bam or .cram).

    Returns:
        tuple[str, str, str]: Paths to:
            - downloaded main file (.bam or .cram),
            - corresponding index file (.bai or .crai),
            - temporary directory containing both files.

    Raises:
        Exception: If the file type is unsupported (not .bam or .cram).
    """

    sample_name = Path(key).stem

    suffix = Path(key).suffix

    if suffix == '.cram':
        index_key = '.crai'
    elif suffix == '.bam':
        index_key = '.bai'
    else:
        shared_log.logger.error(f"Unsupported file type: {suffix}")
        raise ValueError(f"Unsupported file type: {suffix}")
    
    base_data_dir = Path('data')

    # Ensure that if the file exists already, it doens't have to download again
    for temp_dir in base_data_dir.iterdir():

        if temp_dir.is_dir():
            cram_path = temp_dir / f'{sample_name}{suffix}'
            index_path = temp_dir / f'{sample_name}{suffix}{index_key}'

            if cram_path.exists() and index_path.exists():
                shared_log.logger.info(f"Files already download: {sample_name}")
                download_cache[sample_name] = (str(cram_path), str(index_path), str(temp_dir))
                return str(cram_path), str(index_path), str(temp_dir)

    tempdir = tempfile.mkdtemp(dir=Path('data'))

    cram_path = os.path.join(tempdir, f'{sample_name}{suffix}')
    index_path = os.path.join(tempdir, f'{sample_name}{suffix}{index_key}')

    shared_log.logger.info(f"Downloading file: {cram_path}")
    with open(cram_path, 'wb') as file:
        s3_client.download_fileobj(bucket.Bucket, key, file)

    index = f'{key}{index_key}'
    shared_log.logger.info(f"Downloading file: {index_path}")
    with open(index_path, 'wb') as file:
        s3_client.download_fileobj(bucket.Bucket, index, file)

    download_cache[sample_name] = (str(cram_path), str(index_path), str(temp_dir))

    return str(cram_path), str(index_path), str(tempdir)