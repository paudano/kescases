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

STREP_ALL_NOFIGURES = [
    'local/strep/summary/benchmarks/ikc_segment_size.tab',
    'local/strep/summary/benchmarks/kestrel_runtime.tab',
    'local/strep/summary/benchmarks/gatk_runtime.tab',
    'local/strep/summary/benchmarks/assemble_runtime.tab',
    'local/strep/summary/benchmarks/kestrel_trace.tab',
    'local/strep/summary/benchmarks/gatk_trace.tab',
    'local/strep/summary/benchmarks/assemble_trace.tab'
]

STREP_ALL_FIGURES = [
    'local/strep/summary/plots/bm/memory_trace_bm.pdf',
    'local/strep/summary/plots/bm/runtime_cpu.pdf',
    'local/strep/summary/plots/bm/runtime_cpu_noasm.pdf',
    'local/strep/summary/plots/phylo/phylo_variant.pdf',
    'local/strep/summary/plots/phylo/phylo_sero.pdf'
]

STREP_ALL = STREP_ALL_NOFIGURES + STREP_ALL_FIGURES

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

rule strep_all_nofigures:
    input: STREP_ALL_NOFIGURES

rule strep_all_figures:
    input: STREP_ALL_FIGURES

rule strep_all:
    input: STREP_ALL
