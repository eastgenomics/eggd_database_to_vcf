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
    Parse and validate command-line arguments for the VCF generation workflow.
    
    Parameters:
        None
    
    Returns:
        args (argparse.Namespace): Parsed command-line arguments with attributes:
            - database (str): Source database type, either "inca" or "variant_store".
            - input_file (str): Path to the CSV export of the variant database.
            - output_filename (str or None): Desired output VCF filename, if provided.
            - genome_build (str): Reference genome build, either "GRCh37" or "GRCh38".
            - probeset (str or None): Optional probeset filter for inca data, either "germline" or "somatic".
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
    Normalize and clean the input CSV into a standardized DataFrame suitable for downstream VCF generation.
    
    Performs database-specific normalizations, renames and reorders key columns so `CHROM`, `POS`, `REF`, `ALT` are the first four columns, and cleans string cells by removing internal newlines and trimming whitespace.
    
    Parameters:
        database (str): Source database type; must be "inca" or "variant_store". Determines column mappings and database-specific normalizations.
        input_file (str): Path to the input CSV file.
        genome_build (str): Genome build indicator ("GRCh37" or "GRCh38"); for "GRCh38" uses the build-specific start column when available.
    
    Returns:
        pd.DataFrame: Cleaned DataFrame with `CHROM`, `POS`, `REF`, `ALT` as the first four columns, other original columns preserved; string cells have internal newlines removed and are trimmed. For "inca", `date_last_evaluated` is parsed to datetimes and classification fields have spaces replaced with underscores. For "variant_store", square brackets are stripped from the `alternatealleles` values.
    """
    df = pd.read_csv(
        input_file,
        delimiter=",",
        low_memory=False,
    )

    if database == 'inca':
        df["date_last_evaluated"] = pd.to_datetime(
            df["date_last_evaluated"], errors="coerce")
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
        df['alternatealleles'] = df['alternatealleles'].map(
            lambda x: x.lstrip('[').rstrip(']')
            if isinstance(x, str) else x)
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
    Filter a cleaned DataFrame to interpreted variants matching the requested genome build and optional probeset.
    
    Parameters:
        cleaned_csv (pd.DataFrame): Cleaned variant table containing at least the columns
            "interpreted", "germline_classification", "oncogenicity_classification",
            "ref_genome", "ref_genome_38", "allele_origin", and "date_last_evaluated".
        probeset (str | None): If provided, filter rows where `allele_origin` equals this value (case-insensitive).
        genome_build (str): "GRCh37" or "GRCh38" determining which reference column to use
            ("ref_genome" for GRCh37, "ref_genome_38" for GRCh38).
    
    Returns:
        pd.DataFrame: Rows where `interpreted` is "yes", matching the specified genome build
        and probeset (if given), with duplicates removed and only rows that have a
        non-null `date_last_evaluated`.
    
    Raises:
        ValueError: If any row has both `germline_classification` and `oncogenicity_classification` null.
    """
    interpreted_df = cleaned_csv[cleaned_csv["interpreted"].str.lower() == "yes"]
    CLASSIFICATION_ERROR = "Both germline and oncogenicity classification are null in at least one row."
    if (
        interpreted_df["germline_classification"].isnull()
        & interpreted_df["oncogenicity_classification"].isnull()
    ).any():
        raise ValueError(CLASSIFICATION_ERROR)

    if genome_build == "GRCh37":
        prefiltered_df = interpreted_df[
            interpreted_df["ref_genome"].str.contains("grch37", na=False, case=False)
        ]
    else:
        prefiltered_df = interpreted_df[
            interpreted_df["ref_genome_38"].str.contains("grch38", na=False, case=False)
        ]

    if probeset:
        filtered_df = prefiltered_df.loc[
            prefiltered_df["allele_origin"].str.lower() == probeset]
    else:
        filtered_df = prefiltered_df

    filtered_df = filtered_df.drop_duplicates()
    filtered_df = filtered_df.dropna(subset=["date_last_evaluated"])

    return filtered_df


def get_latest_entry(sub_df) -> pd.Series:
    """
    Select the row with the most recent `date_last_evaluated` from a grouped DataFrame.
    
    Parameters:
        sub_df (pd.DataFrame): DataFrame representing a single group of variants (CHROM, POS, REF, ALT). Must contain a `date_last_evaluated` column.
    
    Returns:
        pd.Series: The row with the maximum `date_last_evaluated`.
    """
    latest_idx = sub_df["date_last_evaluated"].idxmax()
    latest_entry = sub_df.loc[latest_idx]
    return latest_entry


def aggregate_hgvs(hgvs_series) -> str:
    """
    Aggregate NM_*:c.* HGVSc entries from a series of variant-store attribute strings into a single pipe-separated string.
    
    Parameters:
        hgvs_series (pd.Series): Series of attribute strings from variant_store rows; non-string and missing values are ignored.
    
    Returns:
        str: Unique HGVSc values joined with "|" (empty string if none found).
    """
    all_hgvs = []
    # extract NM_*:c.* value from a string of pipe-separated values
    pattern = re.compile(r"(?<=\|)(NM_.[^\|]*:c\..*?)(?=\|)")

    for attr_string in hgvs_series.dropna():
        if isinstance(attr_string, str):
            all_hgvs += pattern.findall(attr_string)
    uniq_hgvs = list(set([x for x in all_hgvs if x]))

    return "|".join(uniq_hgvs)


def format_total_classifications(classifications) -> str:
    """
    Format classification counts as a pipe-separated string.
    
    Counts occurrences of each classification in the provided pandas Series and returns them as `classification(count)` entries joined by `|`.
    
    Parameters:
        classifications (pd.Series): Series of classification labels (e.g., germline or oncogenicity).
    
    Returns:
        str: Formatted string of classifications with counts, e.g. "Pathogenic(3)|Likely_pathogenic(1)".
    """
    counts = classifications.value_counts()
    formatted_counts = [
        f"{classification}({count})" for classification, count in counts.items()
    ]
    return "|".join(formatted_counts)


def sort_aggregated_data(aggregated_df) -> pd.DataFrame:
    """
    Order chromosomes as 1–22, X, Y and return the DataFrame sorted by CHROM, POS, REF, and ALT.
    
    Parameters:
        aggregated_df (pd.DataFrame): Aggregated variant data containing at least the columns `CHROM`, `POS`, `REF`, and `ALT`.
    
    Returns:
        pd.DataFrame: The input DataFrame with `CHROM` coerced to the ordered chromosome sequence 1–22, X, Y and rows sorted by `CHROM`, `POS`, `REF`, then `ALT`.
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
    Aggregate per-variant records into a single-row summary table and write the result to a tab-separated file.
    
    For db == "inca", each variant row contains the latest germline and oncogenicity classifications, the latest evaluation date and specimen id, counts of classification occurrences, and aggregated HGVSc strings. For db == "variant_store", each variant row contains aggregated HGVSc extracted from the attributes field, the number of samples with that variant in the group, and the total number of unique samples across the input.
    
    Parameters:
        db (str): Source database type; expected values are "inca" or "variant_store" and determine the aggregation schema.
        probeset_df (pd.DataFrame): DataFrame of variant records already filtered by probeset/genome build.
        aggregated_database (str): Path to write the aggregated TSV output (no header, tab-separated).
    
    Returns:
        pd.DataFrame: Aggregated DataFrame with one row per unique (CHROM, POS, REF, ALT) and columns appropriate to the selected `db`.
    
    Raises:
        AssertionError: If no variants are produced after aggregation.
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
    VARIANT_ERROR = "There are no variants to process. Please check inputs and filters."
    if aggregated_df.empty:
        raise AssertionError(VARIANT_ERROR)

    aggregated_df["POS"] = aggregated_df["POS"].astype("Int64")
    aggregated_df = sort_aggregated_data(aggregated_df)
    aggregated_df.to_csv(aggregated_database, sep="\t", index=False, header=False)
    return aggregated_df


def intialise_vcf(aggregated_df, minimal_vcf) -> None:
    """
    Create a minimal VCF file containing one record per row from the aggregated DataFrame.
    
    Each output record contains CHROM, POS, ID ('.'), REF, ALT, QUAL ('.'), FILTER ('.'), and INFO ('.'). The file is written to `minimal_vcf` and prefixed with the configured minimal VCF header.
    
    Parameters:
        aggregated_df (pd.DataFrame): Aggregated variant rows; must include columns 'CHROM', 'POS', 'REF', and 'ALT'.
        minimal_vcf (str): Path to the output minimal VCF file.
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
    Write a VCF header file containing INFO field definitions and contig records for the specified database and genome build.
    
    Parameters:
        db (str): Database identifier; expected values are 'inca' or 'variant_store' and determine which INFO field definitions are used.
        genome_build (str): Genome build identifier, either 'GRCh37' or 'GRCh38', used to select the contig records written to the header.
        header_filename (str): Path to the output header file to create.
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
    Annotate a minimal VCF with aggregated INFO fields using bcftools and write the annotated VCF.
    
    Parameters:
        db (str): Database type determining which INFO field definitions to use; expected values are "inca" or "variant_store".
        aggregated_database (str): Path prefix to the aggregated tab-separated file (the function will use `{aggregated_database}.gz` for annotation).
        minimal_vcf (str): Path to the minimal VCF to be annotated.
        header_filename (str): Path to a file containing VCF header lines (INFO and contig definitions) to provide to bcftools.
        output_filename (str): Path where the annotated VCF will be written; the file will be created/overwritten and then indexed.
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
    Upload a local file to the current DNAnexus project's job output folder and return a DNAnexus link to the uploaded file.
    
    Parameters:
        outfile (str): Path to the local file to upload.
    
    Returns:
        str: A DNAnexus link (dx://...) pointing to the uploaded file.
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

    """
    Orchestrates the full CSV-to-VCF conversion workflow and produces an annotated VCF.
    
    Processes the provided input CSV according to the selected `database` and `genome_build`, optionally filters by `probeset`, aggregates unique variants, constructs a minimal VCF, writes a VCF header, annotates the VCF with the aggregated data, and indexes outputs. Validates or generates `output_filename` and, when running in DNAnexus, downloads the input and uploads the final VCF and its index.
    
    Parameters:
        database (str): Source database type; expected values are "inca" or "variant_store".
        input_file (str): Path or remote identifier for the input CSV.
        output_filename (str): Desired output VCF filename; if empty a default is created. Must end with ".vcf".
        genome_build (str): Genome build identifier, e.g. "GRCh37" or "GRCh38".
        probeset (str): Optional probeset filter, e.g. "germline" or "somatic"; used when `database` is "inca".
    
    Returns:
        dict: When running inside DNAnexus, returns a mapping with keys "output_vcf" and "output_index" containing DNAnexus links to the uploaded VCF and its index. Otherwise returns `None`.
    
    Raises:
        ValueError: If `output_filename` is provided but does not end with ".vcf".
    """
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