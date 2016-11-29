"""
Rules for downloading sequence data to analyze.
"""

import os


#######################
### Local Variables ###
#######################

### Reference ###

# Reference file location
STREP_REF = 'local/strep/reference/NC_003028.fasta'
STREP_PBP_BED = 'data/strep/NC_003028.pbp.bed'

# Indexes
STREP_REF_FAI = STREP_REF + '.fai'
STREP_BWA_INDEX_LIST=['{}.{}'.format(STREP_REF, ext) for ext in ('amb', 'ann', 'bwt', 'pac', 'sa')]
STREP_PICARD_INDEX=STREP_REF.replace('.fasta', '.dict')
STREP_RTG_INDEX_FLAG='local/strep/reference/rtg/done'


#################
### Functions ###
#################

def _strep_get_all_sample_file_names(wildcards):

    input_file_names = list()

    for accession in STREP_ACCESSIONS:
        list.append('local/strep/samples/{accession}/{accession}_1.fastq.gz'.format(accession))
        list.append('local/strep/samples/{accession}/{accession}_2.fastq.gz'.format(accession))


#############
### Rules ###
#############

### RTG on Reference ###

# strep_rtg_index_reference
#
# Index reference for RTG.
rule strep_rtg_index_reference:
    input:
        ref=STREP_REF
    output:
        rtg_flag=STREP_RTG_INDEX_FLAG
    log:
        'local/strep/reference/log/strep_rtg_index_reference.log'
    shell:
        """rm -rf $(dirname {output.rtg_flag}); """  # RTG fails if the directory exists, and it is automatically created by snakemake
        """bin/rtg format -o $(dirname {output.rtg_flag}) {input.ref} >{log} 2>&1"""

### Reference ###

# strep_reference_index_picard
#
# Index reference for Picard.
rule strep_reference_index_picard:
    input:
        ref=STREP_REF
    output:
        index=STREP_PICARD_INDEX
    shell:
        """java -jar lib/picard.jar CreateSequenceDictionary R={input.ref} O={output.index}"""

# strep_reference_index_bwa
#
# Index reference for BWA.
rule strep_reference_index_bwa:
    input:
        ref=STREP_REF
    output:
        index_files=STREP_BWA_INDEX_LIST
    shell:
        """bwa index {input.ref}"""

# strep_reference_index_samtools
#
# Index reference for samtools
rule strep_reference_index_samtools:
    input:
        ref=STREP_REF
    output:
        fai=STREP_REF_FAI
    shell:
        """samtools faidx {input.ref}"""

# strep_cp_reference
#
# Copy reference from data directory to the working directory.
rule strep_cp_reference:
    input:
        ref=config['strep']['reference']
    output:
        ref=STREP_REF
    shell:
        """cp -f {input.ref} {output.ref}"""


### Samples ###

# strep_get_sample
#
# Dump paired-end FASTQ files from SRA file.
rule strep_get_sample_fastq:
    input:
        sra='local/strep/temp/samples/{accession}.sra'
    output:
        fastq_1='local/strep/samples/{accession}/{accession}_1.fastq.gz',
        fastq_2='local/strep/samples/{accession}/{accession}_2.fastq.gz'
    log:
        'local/strep/samples/log/{accession}/strep_get_sample_fastq.log'
    run:
        # Download
        out_dir = os.path.dirname(output.fastq_1)
        shell("""fastq-dump --split-files --gzip {input.sra} --outdir {out_dir} >{log} 2>&1""")

# strep_get_sample_sra
#
# Download SRA file for sample.
rule strep_get_sample_sra:
    output:
        sra=temp('local/strep/temp/samples/{accession}.sra')
    log:
        'local/strep/samples/log/{accession}/strep_get_sample_sra.log'
    run:
        sra_url = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{}/{}/{}.sra'.format(
            wildcards.accession[0:6], wildcards.accession, wildcards.accession
        )

        dest_dir = os.path.dirname(output.sra)

        shell('wget -qc {sra_url} -P {dest_dir} -o {log}')
