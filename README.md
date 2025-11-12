# eggd_database_to_vcf

## What does this app do?

Given a database export CSV from either inca or a variant store, aggregates variant information and produces a VEP-compatible VCF.

## What data are required for this app to run?

**Packages**

* Python packages as seen in resources/home/dnanexus/packages.

**Inputs**

* `database [str]` (required): Source of exported CSV file. Options: inca, variant_store
* `genome_build [str]` (required): Genome build, required to populate the VCF header with contigs. Options: GRCh37, GRCh38
* `input_file [str]` (required): Name/path of input CSV file
* `capture [str]` (required for variant store): Assay and panel version used to generate the data. Example: MYE_v3
* `output_filename [str]` (optional): Output filename for annotated VCF
* `probeset [str]` (optional, inca only): Probeset to filter data on. Options: germline, somatic
* `threshold_af [float]` (optional, variant_store only): If provided, includes sample IDs for variants with AF < threshold_af. Range: 0-1

## How to run

```bash
# in DNAnexus
dx run $app_id \
-idatabase=$database \
-igenome_build=$genome_build \
-iinput_file=$input_file \
[ -icapture=$capture ] \
[ -ioutput_filename=$output_filename ] \
[ -iprobeset=$probeset ] \
[ -ithreshold_af=$threshold_af ] \
-y

# locally
python3 generate_vcf.py \
--database $database \
--genome_build $genome_build \
--input_file $input_file \
[ --capture $capture ] \
[ --output_filename $output_filename ] \
[ --probeset $probeset ] \
[ --threshold_af $threshold_af ]
```

## What does this app output?

This app outputs an annotated VCF of aggregated variant information, compatible with VEP.
