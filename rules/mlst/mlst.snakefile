"""
Rules and variables common to the MLST experiment.
"""

from kescaseslib import interval

import pandas as pd


###################
### Definitions ###
###################

### Accessions ###
MLST_ACCESSIONS = list()

with open(config['mlst']['accessions'], 'r') as in_file:
    for line in in_file:
        line = line.strip()

        if not line:
            continue

        MLST_ACCESSIONS.append(line)

def _mlst_get_match_series(accession):
    """
    Get a series of booleans for each gene where 'True' indicates if the MLST and KesMLST
    approaches agreed.
    """

    match_series = pd.read_table(
        'local/mlst/results/{}/mlst_merged_calls.tab'.format(accession),
        header=0,
        index_col=0
    ).apply(lambda row: row['MLST'] == row['KesMLST'], axis=1)

    match_series.name = accession

    return match_series


#############
### Rules ###
#############

#
# Merge MLST results from both callers
#

rule mlst_summary_table:
    input:
        tab=expand('local/mlst/results/{accession}/mlst_merged_calls.tab', accession=MLST_ACCESSIONS)
    output:
        tab='local/mlst/summary/mlst_summary.tab'
    run:

        pd.concat(
            [_mlst_get_match_series(accession) for accession in MLST_ACCESSIONS],
            axis=1
        ).to_csv(output.tab, sep='\t', index_label='Gene')



# mlst_merge_results
#
# Merge results from both pipelines for each sample.
rule mlst_merge_results:
    input:
        asm='local/mlst/results/{accession}/assemble/mlst_calls.tab',
        kes='local/mlst/results/{accession}/kesmlst/mlst_calls.tab'
    output:
        tab='local/mlst/results/{accession}/mlst_merged_calls.tab'
    run:

        # Read dataframes
        df_asm = pd.read_table(input.asm, names=('Gene', 'MLST'), index_col=0)
        df_kes = pd.read_table(input.kes, skiprows=1, names=('Gene', 'KesMLST', 'Distance'), index_col=0)

        df_asm = df_asm.ix[df_asm.index != 'ST', :]

        # Concatenate
        df_merged = pd.concat([df_asm['MLST'], df_kes['KesMLST'], df_kes['Distance']], axis=1)

        # Write
        df_merged.to_csv(output.tab, sep='\t', index_label='Gene')
#
# Run MLST with K-mers
#

# mlst_kesmlst_run_mlst
#
# Run MLST on a k-mer database.
rule mlst_kesmlst_run_mlst:
    input:
        ikc='local/mlst/results/{accession}/kesmlst/kmertable.ikc'
    output:
        tab='local/mlst/results/{accession}/kesmlst/mlst_calls.tab',
        time='local/mlst/results/{accession}/kesmlst/bm/kesmlst.time',
        trace='local/mlst/results/{accession}/kesmlst/bm/kesmlst.trace'
    log:
        'local/mlst/results/{accession}/kesmlst/log/kesmlst.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -Xmx2G -Dlogback.configurationFile=lib/kesmlst/logback.xml """
            """-jar {tools.kesmlst} """
            """data/mlst/db/mlst.fa """
            """{input.ikc} """
            """>{output.tab} """
            """2>{log}"""

# mlst_kesmlst_make_ikc
#
# Read FASTQ files and k-merize them to an IKC (Indexed K-mer Count) file. Kestrel
# will read k-mer data from this file.
rule mlst_kesmlst_make_ikc:
    input:
        fq_1='local/mlst/samples/{accession}/{accession}_1.fastq.gz',
        fq_2='local/mlst/samples/{accession}/{accession}_2.fastq.gz'
    output:
        ikc=protected('local/mlst/results/{accession}/kesmlst/kmertable.ikc'),
        time='local/mlst/results/{accession}/kesmlst/bm/kmertable.time',
        trace='local/mlst/results/{accession}/kesmlst/bm/kmertable.trace'
    log:
        'local/mlst/results/{accession}/kesmlst/log/kmertable.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -Xmx3G -jar {tools.kanalyze} """
            """count """
            """-k 31 """
            """--countfilter=kmercount:5 """
            """--quality=10 """
            """-m ikc """
            """--minsize 15 """
            """-o {output.ikc} """
            """{input.fq_1} {input.fq_2} """
            """> {log}"""


#
# Run MLST from assembly
#

# mlst_assemble_parse_results
#
# Summarize MLST calls.
rule mlst_assemble_parse_results:
    input:
        csv='local/mlst/results/{accession}/assemble/mlst/MLST.csv'
    output:
        tab='local/mlst/results/{accession}/assemble/mlst_calls.tab'
    run:
        df = pd.read_csv(input.csv)

        if df.shape[0] != 1:
            call_series = pd.Series({'pgm': 'NA', 'fumC': 'NA', 'aroE': 'NA', 'gdh': 'NA', 'pdhC': 'NA', 'adk': 'NA', 'abcZ': 'NA', 'ST': 'NA'})
        else:
            call_series = df.ix[0, ['pgm', 'fumC', 'aroE', 'gdh', 'pdhC', 'adk', 'abcZ', 'ST']]

        call_series.to_csv(output.tab, sep='\t')

# mlst_assemble_run_mlst
#
# Run MLST on assembled scaffolds.
rule mlst_assemble_run_mlst:
    input:
        scaffolds='local/mlst/temp/{accession}/assemble/genome/scaffolds.fasta'
    output:
        csv=protected('local/mlst/results/{accession}/assemble/mlst/MLST.csv'),
        time='local/mlst/results/{accession}/assemble/bm/mlst_asm.time',
        trace='local/mlst/results/{accession}/assemble/bm/mlst_asm.trace'
    log:
        'local/mlst/results/{accession}/assemble/log/mlst.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/python2 scripts/mlst/run_MLST.py """
            """-o $(dirname {output.csv}) """
            """-i data/mlst/db """
            """-g $(dirname {input.scaffolds}) """
            """-p data/mlst/db/mlst_profiles.txt """
            """--blast_exe bin/blastn """
            """--force """
            """> {log} 2>&1"""

# mlst_assemble_link_scaffolds
#
# The MLST script reads all FASTA Files from a directory as input. Link the scaffolds fasta to its own temp directory
# for input into the MLST script (for rule mlst_run_from_assembly).
rule mlst_assemble_link_scaffolds:
    input:
        scaffolds='local/mlst/results/{accession}/assemble/spades/scaffolds.fasta'
    output:
        scaffolds=temp('local/mlst/temp/{accession}/assemble/genome/scaffolds.fasta')
    shell:
        """ln -sf $(readlink -f {input.scaffolds}) {output.scaffolds}"""

# mlst_assemble_run_assembly
#
# Assemble contigs for a sample.
rule mlst_assemble_run_assembly:
    input:
        fq_1='local/mlst/samples/{accession}/{accession}_1.fastq.gz',
        fq_2='local/mlst/samples/{accession}/{accession}_2.fastq.gz'
    output:
        scaffolds=protected('local/mlst/results/{accession}/assemble/spades/scaffolds.fasta'),
        time='local/mlst/results/{accession}/assemble/bm/spades.time',
        trace='local/mlst/results/{accession}/assemble/bm/spades.trace'
    log:
        'local/mlst/results/{accession}/assemble/log/spades.log'
    shell:
        """rm -rf $(dirname {output.scaffolds}); """
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/spades.py """
            """-k 27,49,71,93,115,127 """
            """--careful """
            """-1 {input.fq_1} -2 {input.fq_2} """
            """-o $(dirname {output.scaffolds}) """
            """> {log}"""
