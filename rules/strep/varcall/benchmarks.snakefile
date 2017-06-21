"""
Merge benchmark data.
"""

import pandas as pd
from kescaseslib import bm


#############
### Rules ###
#############

#
# Memory usage
#

# strep_bm_get_asm_mem
#
# Merge runtimes for each step of the assembly pipeline into one table.
rule strep_bm_get_asm_mem:
    input:
        bm_asm='local/strep/results/{accession}/assemble/bm/spades.trace',
        bm_align='local/strep/results/{accession}/assemble/bm/bwa.trace',
        bm_sort='local/strep/results/{accession}/assemble/bm/sortsam.trace',
        bm_index='local/strep/results/{accession}/assemble/bm/indexbam.trace',
        bm_pileup='local/strep/results/{accession}/assemble/bm/pileup.trace'
    output:
        tab='local/strep/results/{accession}/assemble/bm/trace.tab'
    run:

        pd.concat(
            [
                bm.MaxMemUsage(input.bm_asm).to_series('assemble'),
                bm.MaxMemUsage(input.bm_align).to_series('align'),
                bm.MaxMemUsage(input.bm_sort).to_series('sort'),
                bm.MaxMemUsage(input.bm_index).to_series('index'),
                bm.MaxMemUsage(input.bm_pileup).to_series('pileup')
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')

# strep_bm_get_gatk_mem
#
# Merge max memory usage for each step of the GATK into one table.
rule strep_bm_get_gatk_mem:
    input:
        bm_align='local/strep/results/{accession}/gatk/bm/bwa.trace',
        bm_sort='local/strep/results/{accession}/gatk/bm/sort_sam.trace',
        bm_mark='local/strep/results/{accession}/gatk/bm/mark_dup.trace',
        bm_index='local/strep/results/{accession}/gatk/bm/index_mark_dup.trace',
        bm_target='local/strep/results/{accession}/gatk/bm/realign_target.trace',
        bm_realign='local/strep/results/{accession}/gatk/bm/realign.trace',
        bm_var='local/strep/results/{accession}/gatk/bm/variants.trace'
    output:
        tab='local/strep/results/{accession}/gatk/bm/trace.tab'
    run:

        pd.concat(
            [
                bm.MaxMemUsage(input.bm_align).to_series('align'),
                bm.MaxMemUsage(input.bm_sort).to_series('sort'),
                bm.MaxMemUsage(input.bm_mark).to_series('mark'),
                bm.MaxMemUsage(input.bm_index).to_series('index'),
                bm.MaxMemUsage(input.bm_target).to_series('target'),
                bm.MaxMemUsage(input.bm_realign).to_series('realign'),
                bm.MaxMemUsage(input.bm_var).to_series('variant')
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')

# strep_bm_get_kestrel_mem
#
# Merge max memory usage for each step in the Kestrel pipeline.
rule strep_bm_get_kestrel_mem:
    input:
        bm_table='local/strep/results/{accession}/kestrel/bm/kmertable.trace',
        bm_var='local/strep/results/{accession}/kestrel/bm/variants.trace'
    output:
        tab='local/strep/results/{accession}/kestrel/bm/trace.tab'
    run:

        pd.concat(
            [
                bm.MaxMemUsage(input.bm_table).to_series('table'),
                bm.MaxMemUsage(input.bm_var).to_series('variant')
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')


#
# Runtimes
#

# strep_bm_get_asm_runtime
#
# Merge runtimes for each step of the assembly pipeline into one table.
rule strep_bm_get_asm_runtime:
    input:
        bm_asm='local/strep/results/{accession}/assemble/bm/spades.time',
        bm_align='local/strep/results/{accession}/assemble/bm/bwa.time',
        bm_sort='local/strep/results/{accession}/assemble/bm/sortsam.time',
        bm_index='local/strep/results/{accession}/assemble/bm/indexbam.time',
        bm_pileup='local/strep/results/{accession}/assemble/bm/pileup.time'
    output:
        tab='local/strep/results/{accession}/assemble/bm/time.tab'
    run:

        pd.concat(
            [
                bm.RunTime(input.bm_asm).to_series('assemble'),
                bm.RunTime(input.bm_align).to_series('align'),
                bm.RunTime(input.bm_sort).to_series('sort'),
                bm.RunTime(input.bm_index).to_series('index'),
                bm.RunTime(input.bm_pileup).to_series('pileup')
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')

# strep_bm_get_gatk_runtime
#
# Merge runtimes for each step of the GATK into one table.
rule strep_bm_get_gatk_runtime:
    input:
        bm_align='local/strep/results/{accession}/gatk/bm/bwa.time',
        bm_sort='local/strep/results/{accession}/gatk/bm/sort_sam.time',
        bm_mark='local/strep/results/{accession}/gatk/bm/mark_dup.time',
        bm_index='local/strep/results/{accession}/gatk/bm/index_mark_dup.time',
        bm_target='local/strep/results/{accession}/gatk/bm/realign_target.time',
        bm_realign='local/strep/results/{accession}/gatk/bm/realign.time',
        bm_var='local/strep/results/{accession}/gatk/bm/variants.time'
    output:
        tab='local/strep/results/{accession}/gatk/bm/time.tab'
    run:

        pd.concat(
            [
                bm.RunTime(input.bm_align).to_series('align'),
                bm.RunTime(input.bm_sort).to_series('sort'),
                bm.RunTime(input.bm_mark).to_series('mark'),
                bm.RunTime(input.bm_index).to_series('index'),
                bm.RunTime(input.bm_target).to_series('target'),
                bm.RunTime(input.bm_realign).to_series('realign'),
                bm.RunTime(input.bm_var).to_series('variant')
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')

# strep_bm_get_kestrel_runtime
#
# Merge runtimes for each step of the Kestrel into one table.
rule strep_bm_get_kestrel_runtime:
    input:
        bm_table='local/strep/results/{accession}/kestrel/bm/kmertable.time',
        bm_var='local/strep/results/{accession}/kestrel/bm/variants.time'
    output:
        tab='local/strep/results/{accession}/kestrel/bm/time.tab'
    run:

        pd.concat(
            [
                bm.RunTime(input.bm_table).to_series('table'),
                bm.RunTime(input.bm_var).to_series('variant')
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')
