#!/usr/bin/env python3
"""
Takes a CSV file exported from inca or a variant store, aggregates data for
unique variants, and produces a VCF compatible with VEP.
"""

from glob import glob
import os
import subprocess
import sys

# if running in DNAnexus
if os.path.exists("/home/dnanexus"):
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "--no-index", "--no-deps"] + glob("packages/*")
    )

import argparse
import dxpy
from utils.parse_input import (
    clean_csv,
    filter_probeset,
    aggregate_uniq_vars
)
from utils.generate_output import (
    create_output_filename,
    initialise_vcf,
    write_vcf_header,
    write_rename_file,
    index_file,
    bcftools_annotate_vcf
)
from utils.dxpy_functions import (
    download_input_file,
    upload_output_file
)


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    Returns
    ----------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Create a VCF from a variant database CSV export file"
    )

    parser.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        choices=["inca", "variant_store"],
        help="Type of database (inca or variant_store) used to generate input data"
    )

    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        required=True,
        help="CSV file export of variant database"
    )

    parser.add_argument(
        "-o",
        "--output_filename",
        type=str,
        help="Output VCF filename"
    )

    parser.add_argument(
        "-g",
        "--genome_build",
        type=str,
        required=True,
        choices=["GRCh37", "GRCh38"],
        help="Reference genome build used for analysis"
    )

    parser.add_argument(
        "-set",
        "--probeset",
        type=str,
        choices=["germline", "somatic"],
        help="Probeset (germline or somatic) to filter inca data on"
    )

    parser.add_argument(
        "-c",
        "--capture",
        type=str,
        help="Assay and panel version used for capture e.g. MYE_v3 (variant_store only)"
    )

    args = parser.parse_args()

    return args


@dxpy.entry_point("main")
def main(database: str, input_file: str, output_filename: str,
         genome_build: str, probeset: str, capture : str):

    if os.path.exists("/home/dnanexus"):
        input_file = download_input_file(input_file)

    OUTPUT_FILENAME_ERROR = "Output filename must end with '.vcf'"
    if not output_filename:
        output_filename = create_output_filename(
            database, genome_build, probeset, capture)
    elif not output_filename.endswith(".vcf"):
        raise ValueError(OUTPUT_FILENAME_ERROR)

    valid_captures = ('CEN', 'WES', 'MYE', 'TSO500', 'HRD', 'PCAN')
    NO_CAPTURE_ERROR = "Capture (assay and panel version, e.g. MYE_v3) must be supplied for variant store data"
    INVALID_CAPTURE_ERROR = f"Capture must start with one of: {', '.join(valid_captures)}"
    if database == 'variant_store':
        if not capture:
            raise ValueError(NO_CAPTURE_ERROR)
        elif capture and not capture.startswith(valid_captures):
            raise ValueError(INVALID_CAPTURE_ERROR)

    aggregated_database = "aggregated_database.tsv"
    minimal_vcf = "minimal_vcf.vcf"
    header_filename = "header.vcf"
    temp_vcf = "initial_vcf.vcf"
    renaming_file = "renaming_annotation.txt"

    initial_df = clean_csv(database, input_file, genome_build)
    if database == 'inca':
        filtered_df = filter_probeset(initial_df, probeset, genome_build)
    elif database == 'variant_store':
        filtered_df = initial_df
    aggregated_df = aggregate_uniq_vars(database, filtered_df, aggregated_database)

    initialise_vcf(aggregated_df, minimal_vcf)
    write_vcf_header(database, genome_build, header_filename)
    if database == 'variant_store':
        write_rename_file(renaming_file, capture)
    index_file(aggregated_database)
    bcftools_annotate_vcf(
        database, aggregated_database, minimal_vcf, header_filename,
        temp_vcf, renaming_file, output_filename
    )

    if os.path.exists("/home/dnanexus"):
        output = {}
        output["output_vcf"] = upload_output_file(f"{output_filename}.gz")
        output["output_index"] = upload_output_file(f"{output_filename}.gz.tbi")

        return output


if os.path.exists("/home/dnanexus"):
    dxpy.run()
elif __name__ == "__main__":
    args = parse_args()
    main(args.database, args.input_file, args.output_filename,
         args.genome_build, args.probeset, args.capture)
