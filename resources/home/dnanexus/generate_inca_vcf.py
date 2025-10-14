from glob import glob
from datetime import datetime
import os
import config
import subprocess

if os.path.exists("/home/dnanexus"):
    # running in DNAnexus
    subprocess.check_call(
        ["pip", "install", "--no-index", "--no-deps"] + glob("packages/*")
    )

import pandas as pd
import argparse
import dxpy
import pysam
import pysam.bcftools
from utils.parse_input import (
    clean_csv,
    filter_probeset,
    aggregate_uniq_vars
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

    args = parser.parse_args()

    return args


def intialise_vcf(aggregated_df, minimal_vcf) -> None:
    """
    Initialise minimal VCF with CHROM, POS, ID, REF, ALT with minimal header

    Parameters
    ----------
    aggregated_df : pd.DataFrame
        Dataframe of aggregated data
    minimal_vcf : str
        Output filename for the minimal VCF
    """
    vcf_lines = []
    for _, row in aggregated_df.iterrows():
        vcf_line = (
            f"{row['CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t.\t.\t."
        )
        vcf_lines.append(vcf_line)

    with open(minimal_vcf, "w") as vcf_file:
        vcf_file.write(config.MINIMAL_VCF_HEADER)
        vcf_file.write("\n".join(vcf_lines) + "\n")


def write_vcf_header(db, genome_build, header_filename) -> None:
    """
    Write VCF header by populating INFO fields and specifying contigs

    Parameters
    ----------
    genome_build : str
        Genome build to specify contigs in header
    header_filename : str
        Output filename for the VCF header
    """
    if db == 'inca':
        config_field = config.INFO_FIELDS_INCA
    else:
        config_field = config.INFO_FIELDS_VARSTORE

    with open(header_filename, "w") as header_vcf:
        for field_info in config_field.values():
            info_line = f'##INFO=<ID={field_info["id"]},Number={field_info["number"]},Type={field_info["type"]},Description="{field_info["description"]}">\n'
            header_vcf.write(info_line)

        if genome_build == "GRCh37":
            header_vcf.write(config.GRCh37_CONTIG)
        else:
            header_vcf.write(config.GRCh38_CONTIG)


def index_file(file) -> None:
    """
    Bgzip and index file

    Parameters
    ----------
    file : str
        File to bgzip and index
    """
    pysam.tabix_compress(f"{file}", f"{file}.gz")
    pysam.tabix_index(f"{file}.gz", seq_col=0, start_col=1, end_col=1)


def bcftools_annotate_vcf(
    db, aggregated_database, minimal_vcf, header_filename, output_filename
) -> None:
    """
    Run bcftools annotate to annotate the minimal VCF with the aggregated info

    Parameters
    ----------
    aggregated_database : str
        Output filename of aggregated database
    minimal_vcf : str
        Output filename for the minimal VCF
    header_filename : str
        Output filename for the VCF header
    output_filename : str
        Output filename for annotated VCF
    """
    if db == 'inca':
        config_field = config.INFO_FIELDS_INCA
    else:
        config_field = config.INFO_FIELDS_VARSTORE

    info_fields = ",".join(item["id"] for item in config_field.values())

    # Run bcftools annotate with pysam
    annotate_output = pysam.bcftools.annotate(
        "-a",
        f"{aggregated_database}.gz",
        "-h",
        f"{header_filename}",
        "-c",
        f"CHROM,POS,REF,ALT,{info_fields}",
        f"{minimal_vcf}",
    )
    with open(output_filename, "w") as f:
        f.write(annotate_output)

    index_file(output_filename)


def download_input_file(remote_file) -> str:
    """
    Download given input file with same name as file in project
    Function from vcf_qc.py from eggd_vcf_qc

    Parameters
    ----------
    remote_file : dict
        DNAnexus input file

    Returns
    -------
    str
        name of locally downloaded file
    """
    local_name = dxpy.describe(remote_file).get("name")
    dxpy.bindings.dxfile_functions.download_dxfile(
        dxid=remote_file, filename=local_name
    )

    return local_name


def upload_output_file(outfile) -> None:
    """
    Upload output file to set folder in current project
    Function from vcf_qc.py from eggd_vcf_qc

    Parameters
    ----------
    outfile : str
        name of file to upload
    """
    output_project = os.environ.get("DX_PROJECT_CONTEXT_ID")
    output_folder = (
        dxpy.bindings.dxjob.DXJob(os.environ.get("DX_JOB_ID"))
        .describe()
        .get("folder", "/")
    )
    print(f"\nUploading {outfile} to {output_project}:{output_folder}")

    dxpy.set_workspace_id(output_project)
    dxpy.api.project_new_folder(
        output_project, input_params={"folder": output_folder, "parents": True}
    )

    url_file = dxpy.upload_local_file(
        filename=outfile,
        folder=output_folder,
        wait_on_close=True,
    )

    return dxpy.dxlink(url_file)


def create_output_filename(database, genome_build, probeset):
    """
    Generate an output filename if none is provided

    Parameters
    ----------
    database : str
        inca or variant_store
    genome_build : str
        GRCh37 or GRCh38
    probeset : str
        germline or somatic

    Returns
    -------
    str
        Name for output VCF
    """
    date = datetime.today().strftime("%y%m%d")
    output = f"{date}_{database}_{genome_build}"

    if probeset:
        output += f"_{probeset}"

    output += ".vcf"

    return output


@dxpy.entry_point("main")
def main(database: str, input_file: str, output_filename: str, genome_build: str, probeset: str):

    if os.path.exists("/home/dnanexus"):
        input_file = download_input_file(input_file)

    OUTPUT_FILENAME_ERROR = "Output filename must end with '.vcf'"
    if not output_filename:
        output_filename = create_output_filename(database, genome_build, probeset)
    elif not output_filename.endswith(".vcf"):
        raise ValueError(OUTPUT_FILENAME_ERROR)

    minimal_vcf = "minimal_vcf.vcf"
    header_filename = "header.vcf"
    aggregated_database = "aggregated_database.tsv"

    initial_df = clean_csv(database, input_file, genome_build)
    if database == 'inca':
        filtered_df = filter_probeset(initial_df, probeset, genome_build)
    else:
        filtered_df = initial_df
    aggregated_df = aggregate_uniq_vars(database, filtered_df, aggregated_database)

    intialise_vcf(aggregated_df, minimal_vcf)
    write_vcf_header(database, genome_build, header_filename)
    index_file(aggregated_database)
    bcftools_annotate_vcf(
        database, aggregated_database, minimal_vcf, header_filename, output_filename
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
    main(args.database, args.input_file, args.output_filename, args.genome_build, args.probeset)
