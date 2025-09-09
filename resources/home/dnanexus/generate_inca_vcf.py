from glob import glob
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

def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    Returns
    ----------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(description="Create VCF from previous interpretations")

    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        required=True,
        help=(
            "CSV file exported from the previous interpretations "
            "database"
        ),
    )

    parser.add_argument(
        "-o",
        "--output_filename",
        type=str,
        help=(
            "Output VCF filename"
        ),
    )

    parser.add_argument(
        "-g",
        "--genome_build",
        type=str,
        required=True,
        choices=["GRCh37", "GRCh38"],
        help="Genome build the samples were run in",
    )

    parser.add_argument(
        "-set",
        "--probeset",
        type=str,
        help=(
            "probeset_id or allele_origin to filter. Comma-separated if more than one."
        ),
    )

    args = parser.parse_args()

    return args

def clean_csv(input_file, genome_build) -> pd.DataFrame:
    '''
    Clean up the Inca database CSV by: 
    - Convert to tab separated instead of comma
    - Rename to CHROM, POS, REF, ALT
    - Move CHROM, POS, REF, ALT to be first 4 columns
    - Remove all new lines or tabs within cells

    Parameters
    ----------
    input_file : str
        Filepath to Inca CSV
    genome_build : str
        Genome build to specify columns

    Returns
    -------
    pd.DataFrame
        Dataframe with cleaned up data
    '''
    df = pd.read_csv(
        input_file,
        delimiter=",",
        parse_dates=['date_last_evaluated'],
        low_memory=False,
    )
    columns = {"chromosome": "CHROM",
                "start": "POS",
                "reference_allele": "REF",
                "alternate_allele": "ALT"}
    if genome_build == "GRCh38":
        columns["start_38"] = columns.pop("start")
    df.rename(columns=columns, inplace=True)
    df = df[["CHROM", "POS", "REF", "ALT"] + [ col for col in df.columns if col not in ["CHROM", "POS", "REF", "ALT"]]]
    df = df.applymap(lambda x: x.replace("\n", " ").strip() if isinstance(x, str) else x)

    return df

def filter_probeset(cleaned_csv, probeset, genome_build) -> pd.DataFrame:
    '''
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
    '''
    interpreted_df = cleaned_csv[cleaned_csv['interpreted'].str.lower() == "yes"]
    if (interpreted_df['germline_classification'].isnull() & interpreted_df['oncogenicity_classification'].isnull()).any():
        raise ValueError("Both germline and oncogenicity classification are null in at least one row.")
    
    if genome_build == "GRCh37":
        prefiltered_df = interpreted_df[interpreted_df['ref_genome'].str.contains("grch37", na=False, case=False)]
    else:
        prefiltered_df = interpreted_df[interpreted_df['ref_genome_38'].str.contains("grch38", na=False, case=False)]

    all_dfs = []
    for origin in probeset:
        if origin in ["99347387", "96527893"]:
            column = "probeset_id"
        elif origin.lower() in ["germline", "somatic"]:
            column = "allele_origin"
        else:
            raise ValueError(f"Invalid argument: '{origin}'. Expected one of: germline, somatic, 99347387, or 96527893.")
        filtered_df = prefiltered_df.loc[prefiltered_df[column] == origin]
        all_dfs.append(filtered_df)
    
    probeset_df = pd.concat(all_dfs, ignore_index=True)
    probeset_df = probeset_df.drop_duplicates()

    return probeset_df

def get_latest_entry(sub_df) -> pd.Series:
    '''
    Get latest entry by date

    Parameters
    ----------
    sub_df : pd.DataFrame
        Dataframe per group of CHROM, POS, REF, ALT

    Returns
    -------
    pd.Series
        Latest entry per group by date
    '''
    latest_idx = sub_df['date_last_evaluated'].idxmax()
    latest_entry = sub_df.loc[latest_idx]
    return latest_entry

def aggregate_hgvs(hgvs_series) -> str:
    '''
    Aggregates all unique HGVS

    Parameters
    ----------
    hgvs_series : pd.Series
        HGVSc per variant

    Returns
    -------
    str
        All HGVSc per variant joined
    '''
    unique_hgvs = hgvs_series.dropna().unique()
    return "|".join(unique_hgvs)

def format_total_classifications(classifications) -> str:
    '''
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
    '''
    counts = classifications.value_counts()
    formatted_counts = [f"{classification}({count})" for classification, count in counts.items()]
    return "|".join(formatted_counts)

def sort_aggregated_data(aggregated_df) -> pd.DataFrame:
    '''
    Sort aggregate data
    
    Parameters
    ----------
    aggregated_df : pd.DataFrame
        Dataframe of aggregated data

    Returns
    -------
    pd.DataFrame
        Dataframe sorted by CHROM and POS
    '''
    # Define chromosome order: numeric first, then X and Y and sort
    chromosome_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    aggregated_df['CHROM'] = pd.Categorical(aggregated_df['CHROM'], categories=chromosome_order, ordered=True)
    aggregated_df = aggregated_df.sort_values(by=['CHROM', 'POS'])

    return aggregated_df

def aggregate_uniq_vars(probeset_df, probeset, aggregated_database) -> pd.DataFrame:
    '''
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
    '''
    probeset_df.loc[:, 'germline_classification'] = probeset_df['germline_classification'].str.replace(' ', '_')
    probeset_df.loc[:, 'oncogenicity_classification'] = probeset_df['oncogenicity_classification'].str.replace(' ', '_')
    probeset_df.loc[:, 'CHROM'] = probeset_df['CHROM'].str.replace(' ', '')
    probeset_df = probeset_df.dropna(subset=['date_last_evaluated'])

    aggregated_data = []
    for origin in probeset:
        if origin in ["germline", "99347387"]:
            classification = "germline"
        else:
            classification = "oncogenicity"

        grouped = probeset_df.groupby(['CHROM', 'POS', 'REF', 'ALT', f'{classification}_classification'])

        for _, group in grouped:
            latest_entry = get_latest_entry(group)
            latest_germline = latest_entry['germline_classification']
            latest_oncogenicity = latest_entry['oncogenicity_classification']
            latest_date = latest_entry['date_last_evaluated']
            latest_sample_id = latest_entry['specimen_id']
            hgvs = aggregate_hgvs(group['hgvsc'])
            total_germline = format_total_classifications(group['germline_classification'])
            total_oncogenicity = format_total_classifications(group['oncogenicity_classification'])

            aggregated_data.append({
                'CHROM': latest_entry['CHROM'],
                'POS': latest_entry['POS'],
                'REF': latest_entry['REF'],
                'ALT': latest_entry['ALT'],
                'latest_germline': latest_germline,
                'latest_oncogenicity': latest_oncogenicity,
                'latest_date': latest_date,
                'latest_sample_id': latest_sample_id,
                'total_germline': total_germline,
                'total_oncogenicity': total_oncogenicity,
                'aggregated_hgvs': hgvs
        })

    
    aggregated_df = pd.DataFrame(aggregated_data)

    if aggregated_df.empty:
        raise AssertionError("There are no variants to process. Please check inputs and filters.")
    
    aggregated_df['POS'] = aggregated_df['POS'].astype('Int64')
    aggregated_df = sort_aggregated_data(aggregated_df)
    aggregated_df.to_csv(aggregated_database, sep="\t", index=False, header=False)
    return aggregated_df
    

def intialise_vcf(aggregated_df, minimal_vcf) -> None:
    ''' 
    Initialise minimal VCF with CHROM, POS, ID, REF, ALT with minimal header

    Parameters
    ----------
    aggregated_df : pd.DataFrame
        Dataframe of aggregated data
    minimal_vcf : str
        Output filename for the minimal VCF
    '''
    vcf_lines = []
    for _, row in aggregated_df.iterrows():
        vcf_line = f"{row['CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t.\t.\t."
        vcf_lines.append(vcf_line)

    with open(minimal_vcf, "w") as vcf_file:
        vcf_file.write(config.MINIMAL_VCF_HEADER)
        vcf_file.write("\n".join(vcf_lines) + "\n")

def write_vcf_header(genome_build, header_filename) -> None:
    '''
    Write VCF header by populating INFO fields and specifying contigs

    Parameters
    ----------
    genome_build : str
        Genome build to specify contigs in header
    '''
    with open(header_filename, "w") as header_vcf:
        for field_info in config.INFO_FIELDS.values():
            info_line = f'##INFO=<ID={field_info["id"]},Number={field_info["number"]},Type={field_info["type"]},Description="{field_info["description"]}">\n'
            header_vcf.write(info_line)
        
        if genome_build == "GRCh37":
            header_vcf.write(config.GRCh37_CONTIG)
        else:
            header_vcf.write(config.GRCh38_CONTIG)

def index_file(file) -> None:
    '''
    Bgzip and index file
    
    Parameters
    ----------
    file : str
        File to bgzip and index
    '''
    pysam.tabix_compress(f"{file}", f"{file}.gz")
    pysam.tabix_index(f"{file}.gz", seq_col=0, start_col=1, end_col=1)

def bcftools_annotate_vcf(aggregated_database, minimal_vcf, header_filename, output_filename) -> None:
    '''
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
    '''
    # Run bcftools annotate with pysam
    info_fields = ",".join(item["id"] for item in config.INFO_FIELDS.values())
    annotate_output = pysam.bcftools.annotate(
        "-a", f"{aggregated_database}.gz", 
        "-h", f"{header_filename}", 
        "-c", f"CHROM,POS,REF,ALT,{info_fields}", 
        f"{minimal_vcf}")
    with open(output_filename, 'w') as f:
        f.write(annotate_output)

    index_file(output_filename)

def download_input_file(remote_file) -> str:
    '''
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
    '''
    local_name = dxpy.describe(remote_file).get("name")
    dxpy.bindings.dxfile_functions.download_dxfile(
        dxid=remote_file, filename=local_name
    )

    return local_name

def upload_output_file(outfile) -> None:
    '''
    Upload output file to set folder in current project
    Function from vcf_qc.py from eggd_vcf_qc

    Parameters
    ----------
    outfile : str
        name of file to upload
    '''
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

@dxpy.entry_point("main")
def main(input_file: str,
         output_filename: str,
         genome_build: str,
         probeset: str):
    if os.path.exists("/home/dnanexus"):
        input_file = download_input_file(input_file)

    probeset = [x.strip().lower() for x in probeset.split(",")]
    minimal_vcf = "minimal_vcf.vcf"
    header_filename = "header.vcf"
    aggregated_database = f"{'_'.join(probeset)}_aggregated_database.tsv"

    cleaned_csv = clean_csv(input_file, genome_build)
    probeset_df = filter_probeset(cleaned_csv, probeset, genome_build)
    aggregated_df = aggregate_uniq_vars(probeset_df, probeset, aggregated_database)

    intialise_vcf(aggregated_df, minimal_vcf)
    write_vcf_header(genome_build, header_filename)
    index_file(aggregated_database)
    bcftools_annotate_vcf(aggregated_database, minimal_vcf, header_filename, output_filename)

    if os.path.exists("/home/dnanexus"):
        output = {}
        output["output_vcf"] = upload_output_file(f"{output_filename}.gz")
        output["output_index"] = upload_output_file(f"{output_filename}.gz.tbi")

        return output


if os.path.exists("/home/dnanexus"):
    dxpy.run()
elif __name__ == "__main__":
    args = parse_args()
    main(args.input_file, args.output_filename, args.genome_build, args.probeset)
