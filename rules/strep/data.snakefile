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

# Indexes
STREP_REF_FAI = STREP_REF + '.fai'
STREP_BWA_INDEX_LIST=['{}.{}'.format(STREP_REF, ext) for ext in ('amb', 'ann', 'bwt', 'pac', 'sa')]
STREP_PICARD_INDEX=STREP_REF.replace('.fasta', '.dict')


#############
### Rules ###
#############

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

# strep_link_reference
#
# Link reference from data directory to the working directory.
rule strep_link_reference:
    input:
        ref=config['strep']['reference']
    output:
        ref=STREP_REF
    shell:
        """ln -sf ../../../{input.ref} {output.ref}"""


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
    run:
        # Download
        out_dir = os.path.dirname(output.fastq_1)
        shell("""fastq-dump --split-files --gzip {input.sra} --outdir {out_dir}""")

# strep_get_sample_sra
#
# Download SRA file for sample.
rule strep_get_sample_sra:
    output:
        sra=temp('local/strep/temp/samples/{accession}.sra')
    log:
        'local/strep/log/strep_get_sample_sra/{accession}.log'
    run:
        sra_url = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{}/{}/{}.sra'.format(
            wildcards.accession[0:6], wildcards.accession, wildcards.accession
        )

        dest_dir = os.path.dirname(output.sra)

        shell('wget -qc {sra_url} -P {dest_dir} -o {log}')
