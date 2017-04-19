"""
Master snakefile for the E. coli pipeline.
"""

localrules: ecoli_fetch, ecoli_tables, ecoli_plot


###################
### Definitions ###
###################

# Read accessions
ECOLI_ACCESSIONS = list()

with open(config['ecoli']['accessions'], 'r') as in_file:
    for line in in_file:
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        ECOLI_ACCESSIONS.append(line)


###############
### Include ###
###############

include: 'data.snakefile'
include: 'varcall/assembly.snakefile'
include: 'varcall/kestrel.snakefile'
include: 'varcall/quality.snakefile'
include: 'summary.snakefile'


#############
### Rules ###
#############

# ecoli_figures
#
# Make E. coli summary plots.
rule ecoli_figures:
    input:
        'local/ecoli/summary/plots/con_size_hist.pdf'

# ecoli_tables
#
# Make E. coli summary tables.
rule ecoli_tables:
    input:
        'local/ecoli/summary/table/consensus_len.tab',
        'local/ecoli/summary/table/haplotype_len.tab',
        'local/ecoli/summary/table/summary_len.tab',
        'local/ecoli/summary/table/call_stats_con.tab',
        'local/ecoli/summary/table/call_stats_hap.tab',
        'local/ecoli/summary/table/call_stats_summary.tab'

# ecoli_fetch
#
# Download E. coli data.
rule ecoli_fetch:
    input:
        expand('local/ecoli/samples/{accession}.fa.gz', accession=ECOLI_ACCESSIONS)
