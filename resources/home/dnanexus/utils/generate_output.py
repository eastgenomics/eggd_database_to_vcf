"""
Functions to create an output VCF file compatible with VEP
"""

from datetime import datetime
import config
import pandas as pd
import pysam
import pysam.bcftools


def create_output_filename(db, genome, probeset, capture) -> str:
    """
    Generate an output filename if none is provided

    Parameters
    ----------
    db : str
        inca or variant_store
    genome : str
        GRCh37 or GRCh38
    probeset : str
        germline or somatic
    capture : str
        Assay and panel version, e.g. MYE_v3

    Returns
    -------
    str
        Name for output VCF
    """
    date = datetime.today().strftime("%y%m%d")
    output = f"{date}_{db}_{genome}"

    if db == 'inca' and probeset:
        output += f"_{probeset}"
    elif db == 'variant_store' and capture:
        output += f"_{capture}"

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
    elif db == "variant_store":
        config_field = config.INFO_FIELDS_VARSTORE
    else:
        raise ValueError(f"Unsupported database: {db}")

    with open(header_filename, "w") as header_vcf:
        for field_info in config_field.values():
            info_line = f'##INFO=<ID={field_info["id"]},Number={field_info["number"]},Type={field_info["type"]},Description="{field_info["description"]}">\n'
            header_vcf.write(info_line)

        if genome_build == "GRCh37":
            header_vcf.write(config.GRCh37_CONTIG)
        else:
            header_vcf.write(config.GRCh38_CONTIG)


def write_rename_file(renaming_file, capture):
    """Generate a file specifying how INFO fields should be renamed to include
    capture assay and version.

    Args:
        renaming_file : str
            Name of text file to hold INFO fields to change
        capture : str
            Assay and version for variant store data, e.g. MYE_v3
    """

    # each line has the form "<field to rename> <replacement name>"
    rename_lines = [
        f"INFO/VARIANT_PROPORTION {capture}_PROPORTION",
        f"INFO/VARIANT_COUNT {capture}_COUNT",
        f"INFO/TOTAL_SAMPLES {capture}_TOTAL",
        f"INFO/AC_HET {capture}_AC_HET",
        f"INFO/AC_HOM {capture}_AC_HOM",
        f"INFO/AN {capture}_AN",
        f"INFO/AF {capture}_AF",
        f"INFO/SAMPLE_IDS {capture}_IDS",
        ]

    with open(renaming_file, 'w') as writer:
        writer.write("\n".join(rename_lines) + "\n")


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
    db, aggregated_db, minimal_vcf, header_filename, temp_vcf, renaming_file, output_filename
) -> None:
    """
    Run bcftools annotate to annotate the minimal VCF with the aggregated info.
    For variant store data, also use bcftools annotate to rename INFO fields.

    Parameters
    ----------
    db : str
        Type of database (inca or variant_store)
    aggregated_db : str
        Output filename of aggregated database
    minimal_vcf : str
        Output filename for the minimal VCF
    header_filename : str
        Output filename for the VCF header
    temp_vcf : str
        Intermediate file for annotated variant store data prior to renaming
    renaming_file : str
        Defines how INFO fields should be renamed for variant store data
    output_filename : str
        Output filename for annotated VCF
    """
    if db == 'inca':
        config_fields = config.INFO_FIELDS_INCA
        annotation_output = output_filename
    elif db == 'variant_store':
        config_fields = config.INFO_FIELDS_VARSTORE
        annotation_output = temp_vcf

    info_fields = ",".join(f'INFO/{item["id"]}' for item in config_fields.values())

    # Run bcftools annotate with pysam
    pysam.bcftools.annotate(
        "-a", f"{aggregated_db}.gz",
        "-h", f"{header_filename}",
        "-c", f"CHROM,POS,REF,ALT,{info_fields}",
        "-O", "v",
        "-o", f"{annotation_output}",
        f"{minimal_vcf}",
        catch_stdout=False
    )

    # Rename INFO fields to contain assay and version
    if db == 'variant_store':
        pysam.bcftools.annotate(
            "--rename-annots", renaming_file,
            "-o", f"{output_filename}",
            f"{temp_vcf}",
            catch_stdout=False
        )

    index_file(output_filename)
