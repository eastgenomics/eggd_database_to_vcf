from glob import glob
from datetime import datetime
import os
import re
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


def clean_csv(database, input_file, genome_build) -> pd.DataFrame:
    """
    Clean up the input CSV by:
    - Convert to tab separated instead of comma
    - Rename to CHROM, POS, REF, ALT
    - Move CHROM, POS, REF, ALT to be first 4 columns
    - Remove all new lines or tabs within cells

    Parameters
    ----------
    database : str
        inca or variant_store, affects column names
    input_file : str
        Filepath to input CSV
    genome_build : str
        Genome build to specify columns

    Returns
    -------
    pd.DataFrame
        Dataframe with cleaned up data
    """
    df = pd.read_csv(
        input_file,
        delimiter=",",
        low_memory=False,
    )

    if database == 'inca':
        df["date_last_evaluated"] = pd.to_datetime(df["date_last_evaluated"])
        df.loc[:, "germline_classification"] = df[
            "germline_classification"].str.replace(" ", "_")
        df.loc[:, "oncogenicity_classification"] = df[
            "oncogenicity_classification"].str.replace(" ", "_")

        columns = {
            "chromosome": "CHROM",
            "start": "POS",
            "reference_allele": "REF",
            "alternate_allele": "ALT",
        }
        if genome_build == "GRCh38":
            columns["start_38"] = columns.pop("start")

    else:
        df['alternatealleles'] = df['alternatealleles'].map(lambda x: x.lstrip('[').rstrip(']'))
        columns = {
            "contigname": "CHROM",
            "start": "POS",
            "referenceallele": "REF",
            "alternatealleles": "ALT",
            "filters": "FILTER"
        }

    df.rename(columns=columns, inplace=True)
    df = df[
        list(columns.values())
        + [col for col in df.columns if col not in list(columns.values())]
    ]
    df = df.applymap(
        lambda x: x.replace("\n", " ").strip() if isinstance(x, str) else x
    )

    return df


def filter_probeset(cleaned_csv, probeset, genome_build) -> pd.DataFrame:
    """
    Filter cleaned data to interpreted variants for specified germline/somatic probesets

    Parameters
    ----------
    cleaned_csv : pd.DataFrame
        Dataframe with cleaned up data
    probeset : str
        Germline or somatic choice
    genome_build : str
        Genome build to filter

    Returns
    -------
    pd.DataFrame
        Dataframe filtered by probeset
    """
    interpreted_df = cleaned_csv[cleaned_csv["interpreted"].str.lower() == "yes"]
    if (
        interpreted_df["germline_classification"].isnull()
        & interpreted_df["oncogenicity_classification"].isnull()
    ).any():
        raise ValueError(
            "Both germline and oncogenicity classification are null in at least one row."
        )

    if genome_build == "GRCh37":
        prefiltered_df = interpreted_df[
            interpreted_df["ref_genome"].str.contains("grch37", na=False, case=False)
        ]
    else:
        prefiltered_df = interpreted_df[
            interpreted_df["ref_genome_38"].str.contains("grch38", na=False, case=False)
        ]

    if probeset:
        filtered_df = prefiltered_df.loc[prefiltered_df["allele_origin"] == probeset]
    else:
        filtered_df = prefiltered_df

    filtered_df = filtered_df.drop_duplicates()
    filtered_df = filtered_df.dropna(subset=["date_last_evaluated"])

    return filtered_df


def get_latest_entry(sub_df) -> pd.Series:
    """
    Get latest entry by date

    Parameters
    ----------
    sub_df : pd.DataFrame
        Dataframe per group of CHROM, POS, REF, ALT

    Returns
    -------
    pd.Series
        Latest entry per group by date
    """
    latest_idx = sub_df["date_last_evaluated"].idxmax()
    latest_entry = sub_df.loc[latest_idx]
    return latest_entry


def aggregate_hgvs(hgvs_series) -> str:
    """
    Given the 'attributes' column for one variant group from variant store
    data, extract all HGVSc from each row and combine into a single list.

    Parameters
    ----------
    hgvs_series : pd.Series
        'attributes' column from variant store data

    Returns
    -------
    str
        All HGVSc per variant joined
    """
    all_hgvs = []
    for attr_string in hgvs_series:
        # extract the NM_*:c.* value from a string of pipe-separated values
        all_hgvs += re.findall(r"(?<=\|)(NM_.[^\|]*:c\..*?)(?=\|)", attr_string)
    uniq_hgvs = list(set([x for x in all_hgvs if x]))

    return "|".join(uniq_hgvs)


def format_total_classifications(classifications) -> str:
    """
    Counts all classifications, including the latest classification.
    Returns classifications in the format: classification(count)|classification(count)

    Parameters
    ----------
    classifications : pd.Series
        Germline or somatic classifications per variant

    Returns
    -------
    str
        All classifications per variant joined
    """
    counts = classifications.value_counts()
    formatted_counts = [
        f"{classification}({count})" for classification, count in counts.items()
    ]
    return "|".join(formatted_counts)


def sort_aggregated_data(aggregated_df) -> pd.DataFrame:
    """
    Sort aggregate data

    Parameters
    ----------
    aggregated_df : pd.DataFrame
        Dataframe of aggregated data

    Returns
    -------
    pd.DataFrame
        Dataframe sorted by CHROM and POS
    """
    # Define chromosome order: numeric first, then X and Y and sort
    chromosome_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    aggregated_df["CHROM"] = pd.Categorical(
        aggregated_df["CHROM"], categories=chromosome_order, ordered=True
    )
    aggregated_df = aggregated_df.sort_values(by=["CHROM", "POS", "REF", "ALT"])

    return aggregated_df


def aggregate_uniq_vars(db, probeset_df, aggregated_database) -> pd.DataFrame:
    """
    Aggregate data for each unique variant
    Similaritites to create_vcf_from_inca_csv.py by Raymond Miles

    Parameters
    ----------
    probeset_df : pd.DataFrame
        Dataframe filtered by probeset
    probeset : str
        Germline or somatic choice
    aggregated_database : str
        Output filename for aggregated data

    Returns
    -------
    pd.DataFrame
        Dataframe of aggregated data
    """

    aggregated_data = []
    grouped = probeset_df.groupby(["CHROM", "POS", "REF", "ALT"])

    if db == 'variant_store':
        uniq_sample_count = len(probeset_df["sampleid"].dropna().unique())

    for _, group in grouped:
        if db == 'inca':
            latest_entry = get_latest_entry(group)
            latest_germline = latest_entry["germline_classification"]
            latest_oncogenicity = latest_entry["oncogenicity_classification"]
            latest_date = latest_entry["date_last_evaluated"]
            latest_sample_id = latest_entry["specimen_id"]
            hgvs = "|".join(group["hgvsc"].dropna().unique())
            total_germline = format_total_classifications(
                group["germline_classification"]
            )
            total_oncogenicity = format_total_classifications(
                group["oncogenicity_classification"]
            )

            aggregated_data.append(
                {
                    "CHROM": latest_entry["CHROM"],
                    "POS": latest_entry["POS"],
                    "REF": latest_entry["REF"],
                    "ALT": latest_entry["ALT"],
                    "latest_germline": latest_germline,
                    "latest_oncogenicity": latest_oncogenicity,
                    "latest_date": latest_date,
                    "latest_sample_id": latest_sample_id,
                    "total_germline": total_germline,
                    "total_oncogenicity": total_oncogenicity,
                    "aggregated_hgvs": hgvs,
                }
            )

        else:
            hgvs = aggregate_hgvs(group['attributes'])

            aggregated_data.append(
                {
                    "CHROM": group['CHROM'].unique()[0],
                    "POS": group['POS'].unique()[0],
                    "REF": group['REF'].unique()[0],
                    "ALT": group['ALT'].unique()[0],
                    "aggregated_hgvs": hgvs,
                    "variant_sample_count": len(group['sampleid'].dropna().unique()),
                    "total_samples": uniq_sample_count,
                }
            )

    aggregated_df = pd.DataFrame(aggregated_data)

    if aggregated_df.empty:
        raise AssertionError(
            "There are no variants to process. Please check inputs and filters."
        )

    aggregated_df["POS"] = aggregated_df["POS"].astype("Int64")
    aggregated_df = sort_aggregated_data(aggregated_df)
    aggregated_df.to_csv(aggregated_database, sep="\t", index=False, header=False)
    return aggregated_df


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
    """ Generate an output filename if none is provided

    Args:
        database (str): inca or variant_store
        genome_build (str): GRCh37 or GRCh38
        probeset (str): germline or somatic

    Returns:
        str: Name for output VCF
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

    if not output_filename:
        output_filename = create_output_filename(database, genome_build, probeset)
    elif not output_filename.endswith(".vcf"):
        raise ValueError("Output filename must end with '.vcf'")

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
