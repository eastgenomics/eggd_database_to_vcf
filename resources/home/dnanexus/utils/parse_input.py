"""
Functions to read, clean, filter and parse data from a database CSV export
using pandas dataframes.
"""

import re
import pandas as pd


def clean_csv(database, input_file, genome_build) -> pd.DataFrame:
    """
    Clean up the input CSV by:
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

    EMPTY_DF_ERROR = "Imported dataframe is empty"
    if df.empty:
        raise ValueError(EMPTY_DF_ERROR)

    if database == 'inca':
        df["date_last_evaluated"] = pd.to_datetime(
            df["date_last_evaluated"], errors="coerce")
        df.loc[:, "germline_classification"] = df[
            "germline_classification"].str.replace(" ", "_", regex=False)
        df.loc[:, "oncogenicity_classification"] = df[
            "oncogenicity_classification"].str.replace(" ", "_", regex=False)

        columns = {
            "chromosome": "CHROM",
            "start": "POS",
            "reference_allele": "REF",
            "alternate_allele": "ALT",
        }
        if genome_build == "GRCh38":
            if 'start_38' not in df.columns:
                raise KeyError('Missing "start_38" column for GRCh38 INCA export')
            columns["start_38"] = columns.pop("start")

    elif database == 'variant_store':
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
        lambda x: x.replace("\n", " ").replace("\t", " ").strip()
        if isinstance(x, str) else x
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
    # extract NM_*:c.* value from a string of pipe-separated values
    pattern = re.compile(r"(?<=\|)(NM_.[^\|]*:c\..*?)(?=\|)")

    for attr_string in hgvs_series.dropna():
        if isinstance(attr_string, str):
            all_hgvs += pattern.findall(attr_string)
    uniq_hgvs = sorted({x for x in all_hgvs if x})

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

    # Define chromosome order: numeric first, then X and Y
    chromosome_order = [str(i) for i in range(1, 23)] + ["X", "Y"]

    # convert CHROM field to string before sorting
    aggregated_df["CHROM"] = (aggregated_df["CHROM"].astype(str))
    aggregated_df["CHROM"] = pd.Categorical(
        aggregated_df["CHROM"], categories=chromosome_order, ordered=True
    )

    aggregated_df = aggregated_df.sort_values(by=["CHROM", "POS", "REF", "ALT"])

    return aggregated_df


def aggregate_uniq_vars(db, capture, threshold_af, probeset_df, aggregated_database) -> pd.DataFrame:
    """
    Aggregate data for each unique variant
    Similaritites to create_vcf_from_inca_csv.py by Raymond Miles

    Parameters
    ----------
    db : str
        Type of database (inca or variant_store)
    capture : str
        Capture (assay and version) used to generate variant store data, e.g MYE_v3
    threshold_af : float | None
        If set, include SAMPLE_IDS when variant_proportion < threshold_af (variant_store only)
    probeset_df : pd.DataFrame
        Dataframe filtered by probeset
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
        NO_SAMPLES_ERROR = "sampleid column is empty"
        if not uniq_sample_count:
            raise ValueError(NO_SAMPLES_ERROR)

    for _, group in grouped:
        if db == 'inca':
            latest_entry = get_latest_entry(group)
            latest_germline = latest_entry["germline_classification"]
            latest_oncogenicity = latest_entry["oncogenicity_classification"]
            latest_date = latest_entry["date_last_evaluated"]
            latest_sample_id = latest_entry["specimen_id"]
            hgvs = "|".join(sorted(group["hgvsc"].dropna().unique()))
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

        elif db == 'variant_store':
            chrom = group['CHROM'].unique()[0]
            variant_count = len(group['sampleid'].dropna().unique())
            variant_proportion = variant_count / uniq_sample_count

            # calculate AC, AN and AF for non-X/Y germline variants
            ac_het = ac_hom = an = af = ''
            if capture.startswith(('CEN', 'WES')) and chrom not in ['X', 'Y']:
                ac_het = len(group[group['calls'] == '[0, 1]'])
                ac_hom = 2 * len(group[group['calls'] == '[1, 1]'])
                an = 2 * uniq_sample_count
                af = (ac_het + ac_hom) / an

            # include sample ids if a threshold AF is specified
            if threshold_af and (variant_proportion < threshold_af):
                sample_ids = "|".join(sorted(group["sampleid"].dropna().unique()))
            else:
                sample_ids = ""

            aggregated_data.append(
                {
                    "CHROM": chrom,
                    "POS": group['POS'].unique()[0],
                    "REF": group['REF'].unique()[0],
                    "ALT": group['ALT'].unique()[0],
                    "variant_proportion": variant_proportion,
                    "variant_count": variant_count,
                    "total_samples": uniq_sample_count,
                    "ac_het": ac_het,
                    "ac_hom": ac_hom,
                    "an": an,
                    "af": af,
                    "sample_ids": sample_ids,
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
