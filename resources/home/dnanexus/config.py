INFO_FIELDS_INCA = {
    "latest_germline": {
        "id": "LATEST_GERMLINE",
        "number": 1,
        "type": "String",
        "description": "Latest germline classification"},
    "latest_oncogenicity": {
        "id": "LATEST_ONCOGENICITY",
        "number": 1,
        "type": "String",
        "description": "Latest somatic classification"},
    "latest_date": {
        "id": "LATEST_DATE",
        "number": 1,
        "type": "String",
        "description": "Latest date of classification"},
    "latest_sample_id": {
        "id": "LATEST_SAMPLE_ID",
        "number": 1,
        "type": "String",
        "description": "Latest sample ID with classification"},
    "total_germline": {
        "id": "TOTAL_GERMLINE",
        "number": 1,
        "type": "String",
        "description": "Total number of germline classifications"},
    "total_oncogenicity": {
        "id": "TOTAL_ONCOGENICITY",
        "number": 1,
        "type": "String",
        "description": "Total number of oncogenicity classifications"},
    "aggregated_hgvs": {
        "id": "AGGREGATED_HGVS",
        "number": 1,
        "type": "String",
        "description": "Aggregated HGVSc"}
    }

INFO_FIELDS_VARSTORE = {
    "het_count": {
        "id": "HET_COUNT",
        "number": 1,
        "type": "Integer",
        "description": "Number of samples heterozygous for the variant"},
    "hom_count": {
        "id": "HOM_COUNT",
        "number": 1,
        "type": "Integer",
        "description": "Number of samples homozygous for the variant"},
    "allele_count": {
        "id": "AC",
        "number": 1,
        "type": "Integer",
        "description": "Total number of ALT alleles seen"},
    "allele_number": {
        "id": "AN",
        "number": 1,
        "type": "Integer",
        "description": "Total number of alleles seen"},
    "allele_frequency": {
        "id": "AF",
        "number": 1,
        "type": "Float",
        "description": "Allele frequency within the relevant variant store"},
}

GRCh37_CONTIG = '''\
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>\
'''

GRCh38_CONTIG = '''\
##contig=<ID=1,length=248956422,assembly=hg38>
##contig=<ID=2,length=242193529,assembly=hg38>
##contig=<ID=3,length=198295559,assembly=hg38>
##contig=<ID=4,length=190214555,assembly=hg38>
##contig=<ID=5,length=181538259,assembly=hg38>
##contig=<ID=6,length=170805979,assembly=hg38>
##contig=<ID=7,length=159345973,assembly=hg38>
##contig=<ID=8,length=145138636,assembly=hg38>
##contig=<ID=9,length=138394717,assembly=hg38>
##contig=<ID=10,length=133797422,assembly=hg38>
##contig=<ID=11,length=135086622,assembly=hg38>
##contig=<ID=12,length=133275309,assembly=hg38>
##contig=<ID=13,length=114364328,assembly=hg38>
##contig=<ID=14,length=107043718,assembly=hg38>
##contig=<ID=15,length=101991189,assembly=hg38>
##contig=<ID=16,length=90338345,assembly=hg38>
##contig=<ID=17,length=83257441,assembly=hg38>
##contig=<ID=18,length=80373285,assembly=hg38>
##contig=<ID=19,length=58617616,assembly=hg38>
##contig=<ID=20,length=64444167,assembly=hg38>
##contig=<ID=21,length=46709983,assembly=hg38>
##contig=<ID=22,length=50818468,assembly=hg38>
##contig=<ID=X,length=156040895,assembly=hg38>
##contig=<ID=Y,length=57227415,assembly=hg38>
##contig=<ID=M,length=16569,assembly=hg38>\
'''

MINIMAL_VCF_HEADER = '''\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
'''
