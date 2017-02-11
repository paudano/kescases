"""
Rules for downloading sequence data to analyze.
"""

import collections
import gzip
import io
import os

import pandas as pd
import numpy as np

from Bio import SeqIO

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
        input_file_names.append('local/strep/samples/{}/{}_1.fastq.gz'.format(accession, accession))
        input_file_names.append('local/strep/samples/{}/{}_2.fastq.gz'.format(accession, accession))

    return input_file_names


#############
### Rules ###
#############

#
# Sequnece data information
#


rule strep_count_reads_and_bases:
    input:
        fq=_strep_get_all_sample_file_names
    output:
        tab='local/strep/samples/sample_seq_summary.tab'
    run:

        # Create an empty data frame
        df = pd.DataFrame([], columns=['reads', 'bases'], dtype=(np.int64, np.int64))

        # Read sequence data for each accession
        for accession in STREP_ACCESSIONS:
            record_count = 0
            base_count = 0

            # Read base and record count
            for in_file_name in ['local/strep/samples/{}/{}_{}.fastq.gz'.format(accession, accession, n) for n in range(1, 3)]:
                with io.TextIOWrapper(gzip.open(in_file_name)) as fq_file:
                    for seq_record in SeqIO.parse(fq_file, 'fastq'):
                        record_count += 1
                        base_count += len(seq_record.seq)

            # Add to table
            ser = pd.Series({'reads': record_count, 'bases': base_count}, dtype=(np.int64, np.int64))
            ser.name = accession

            df = df.append(ser)

        # Make integers
        df = df.astype(np.int64)

        # Write table
        df.to_csv(output.tab, sep='\t', index=True, index_label='accession')



#
# RTG on Reference
#

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
        """bin/rtg format -o $(dirname {output.rtg_flag}) {input.ref} >{log} 2>&1; """
        """touch {output.rtg_flag}"""

#
# Reference
#

# strep_reference_index_picard
#
# Index reference for Picard.
rule strep_reference_index_picard:
    input:
        ref=STREP_REF
    output:
        index=STREP_PICARD_INDEX
    shell:
        """java -jar {tools.picard} """
            """CreateSequenceDictionary """
            """R={input.ref} """
            """O={output.index}"""

# strep_reference_index_bwa
#
# Index reference for BWA.
rule strep_reference_index_bwa:
    input:
        ref=STREP_REF
    output:
        index_files=STREP_BWA_INDEX_LIST
    shell:
        """bin/bwa index {input.ref}"""

# strep_reference_index_samtools
#
# Index reference for samtools
rule strep_reference_index_samtools:
    input:
        ref=STREP_REF
    output:
        fai=STREP_REF_FAI
    shell:
        """bin/samtools faidx {input.ref}"""

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
