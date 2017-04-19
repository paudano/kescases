"""
Rules and variables common to the strep experiment.
"""

from kescaseslib import interval


###################
### Definitions ###
###################

### Accessions ###
STREP_ACCESSIONS = list()

with open(config['strep']['accessions'], 'r') as in_file:
    for line in in_file:
        line = line.strip()

        if not line:
            continue

        STREP_ACCESSIONS.append(line)

STREP_INTERVAL_CONTAINER = interval.IntervalContainer()
STREP_INTERVAL_CONTAINER.add_bed(config['strep']['pbp_bed'])


###############
### Include ###
###############

include: 'data.snakefile'
include: 'varcall/assembly.snakefile'
include: 'varcall/gatk.snakefile'
include: 'varcall/kestrel.snakefile'
include: 'varcall/quality.snakefile'
include: 'varcall/benchmarks.snakefile'
include: 'summary.snakefile'
include: 'ani.snakefile'


#############
### Rules ###
#############

# strep_figures
#
# Make strep figures.
rule strep_figures:
    input:
        'local/strep/summary/plots/bm/rt/runtime_cpu_by_reads.pdf',
        'local/strep/summary/plots/bm/rt/runtime_cpu_noasm_by_reads.pdf',
        'local/strep/summary/plots/bm/phylo/phylo_sero.pdf',
        'local/strep/summary/plots/bm/phylo/phylo_variant.pdf',
        'local/strep/summary/plots/bm/size/size_bam_vs_ikc.pdf',
        'local/strep/summary/plots/bm/size/size_ikc_bam_by_reads.pdf',
        'local/strep/summary/plots/bm/varjitter/gatk_variant_jitter_all.pdf',
        'local/strep/summary/plots/bm/varjitter/kestrel_variant_jitter_all.pdf'

# strep_tables
#
# Get Strep results tables.
rule strep_tables:
    input:
        'local/strep/summary/kestrel/quality/summary_stats.tab',
        'local/strep/summary/kestrel/variants.tab',
        'local/strep/summary/kestrel/variants_removed.tab',
        'local/strep/summary/gatk/quality/summary_stats.tab',
        'local/strep/summary/gatk/variants.tab',
        'local/strep/summary/gatk/variants_removed.tab',
        'local/strep/summary/benchmarks/assemble_runtime.tab',
        'local/strep/summary/benchmarks/assemble_runtime_summary_byreads.tab',
        'local/strep/summary/benchmarks/assemble_trace.tab',
        'local/strep/summary/benchmarks/gatk_runtime.tab',
        'local/strep/summary/benchmarks/gatk_runtime_summary.tab',
        'local/strep/summary/benchmarks/gatk_runtime_summary_byreads.tab',
        'local/strep/summary/benchmarks/gatk_trace.tab',
        'local/strep/summary/benchmarks/gatk_trace_summary.tab',
        'local/strep/summary/benchmarks/kestrel_runtime_summary.tab',
        'local/strep/summary/benchmarks/kestrel_runtime_summary_byreads.tab',
        'local/strep/summary/benchmarks/kestrel_trace.tab',
        'local/strep/summary/benchmarks/kestrel_trace_summary.tab',
        'local/strep/summary/benchmarks/seq_file_size.tab',
        'local/strep/summary/benchmarks/seq_file_size_summary.tab'


# strep_fetch
#
# Fetch Strep data.
rule strep_fetch:
    input:
        expand('local/strep/samples/{accession}/{accession}_1.fastq.gz', accession=STREP_ACCESSIONS),
        expand('local/strep/samples/{accession}/{accession}_2.fastq.gz', accession=STREP_ACCESSIONS)

