"""
Functions to create an output VCF file compatible with VEP
"""

from datetime import datetime
import config
import pandas as pd
import pysam
import pysam.bcftools


def create_output_filename(database, genome_build, probeset) -> str:
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


def initialise_vcf(df, minimal_vcf) -> None:
    """
    Initialise minimal VCF with CHROM, POS, REF, ALT with minimal header

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe of aggregated data
    minimal_vcf : str
        Output filename for the minimal VCF
    """
    vcf_lines = [
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t."
        for chrom,pos,ref,alt
        in zip(df['CHROM'], df['POS'], df['REF'], df['ALT'])
        ]

    with open(minimal_vcf, "w") as vcf_file:
        vcf_file.write(config.MINIMAL_VCF_HEADER)
        vcf_file.write("\n".join(vcf_lines) + "\n")


def write_vcf_header(db, genome_build, header_filename) -> None:
    """
    Write VCF header by populating INFO fields and specifying contigs

    Parameters
    ----------
    db : str
        Type of database (inca or variant_store)
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
    pysam.tabix_compress(file, f"{file}.gz", force=True)

    # handle vcf and tsv files differently
    if file.endswith(".vcf"):
        pysam.tabix_index(f"{file}.gz", preset="vcf", force=True)
    else:
        pysam.tabix_index(f"{file}.gz", seq_col=0, start_col=1, end_col=1, force=True)


def bcftools_annotate_vcf(
    db, aggregated_database, minimal_vcf, header_filename, output_filename
) -> None:
    """
    Run bcftools annotate to annotate the minimal VCF with the aggregated info

    Parameters
    ----------
    db : str
        Type of database (inca or variant_store)
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
