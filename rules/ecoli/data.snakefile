"""
Data management and preparation.
"""

import pandas as pd


###################
### Definitions ###
###################

# Reference
ECOLI_REF = 'local/ecoli/reference/NC_000913.fasta'

# Indexes
ECOLI_REF_FAI = ECOLI_REF + '.fai'
ECOLI_BWA_INDEX_LIST=['{}.{}'.format(ECOLI_REF, ext) for ext in ('amb', 'ann', 'bwt', 'pac', 'sa')]
ECOLI_RTG_INDEX_FLAG='local/ecoli/reference/rtg/done'
ECOLI_REF_BED = ECOLI_REF + '.bed'


#############
### Rules ###
#############

#
# RTG on Reference
#

# ecoli_rtg_index_reference
#
# Index reference for RTG.
rule ecoli_rtg_index_reference:
    input:
        ref=ECOLI_REF
    output:
        rtg_flag=touch(ECOLI_RTG_INDEX_FLAG)
    log:
        'local/ecoli/reference/log/ecoli_rtg_index_reference.log'
    shell:
        """rm -rf $(dirname {output.rtg_flag}); """  # RTG fails if the directory exists, and it is automatically created by snakemake
        """bin/rtg format -o $(dirname {output.rtg_flag}) {input.ref} >{log} 2>&1"""


#
# Reference
#

# ecoli_reference_make_bed
#
# Make a BED file of the complete reference loci.
rule ecoli_reference_make_bed:
    input:
        fai=ECOLI_REF_FAI
    output:
        bed=ECOLI_REF_BED
    run:

        # Read and get BED columns
        df = pd.read_table(input.fai, header=None)
        df['#CHROM'] = df[0]
        df['POS'] = 0
        df['END'] = df[1]

        # Write
        df.ix[:, ('#CHROM', 'POS', 'END')].to_csv(output.bed, sep='\t', index=False)


# ecoli_reference_index_bwa
#
# Index reference for BWA.
rule ecoli_reference_index_bwa:
    input:
        ref=ECOLI_REF
    output:
        index_files=ECOLI_BWA_INDEX_LIST
    log:
        'local/ecoli/reference/log/bwa_index.log'
    shell:
        """bin/bwa index {input.ref} >{log} 2>&1"""

# ecoli_reference_index_samtools
#
# Index reference for samtools
rule ecoli_reference_index_samtools:
    input:
        ref=ECOLI_REF
    output:
        fai=ECOLI_REF_FAI
    log:
        'local/ecoli/reference/log/samtools_faidx.log'
    shell:
        """bin/samtools faidx {input.ref} >{log} 2>&1"""

# ecoli_cp_reference
#
# Copy reference from data directory to the working directory.
rule ecoli_cp_reference:
    input:
        ref=config['ecoli']['reference']
    output:
        ref=ECOLI_REF
    shell:
        """gunzip -c {input.ref} > {output.ref}"""


#
# Contigs
#

# ecoli_data_extract_contigs
#
# Compress contigs
rule ecoli_data_extract_contigs:
    input:
        tar_gz='local/strep/temp/data/Dataset_S7.tar.gz'
    output:
        fa_gz='local/ecoli/samples/{accession}.fa.gz'
    run:
        shell("""tar -Ozxf {input.tar_gz} Supplement_S7/{wildcards.accession}.scaffolds_min500bp.fa | gzip > {output.fa_gz}""")

# ecoli_data_dl_contigs
#
# Get contigs.
rule ecoli_data_dl_contigs:
    output:
        tar_gz=temp('local/strep/temp/data/Dataset_S7.tar.gz')
    shell:
        """wget -O {output.tar_gz} http://genome.cshlp.org/content/suppl/2014/11/11/gr.180190.114.DC1/Dataset_S7.tar.gz"""
