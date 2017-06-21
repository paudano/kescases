"""
Merge benchmarking information for MLST.
"""

import pandas as pd
from kescaseslib import bm

#############
### Rules ###
#############

#
# Memory usage
#

# mlst_bm_get_mlst_mem
#
# Merge runtimes for each step of the assembly pipeline into one table.
rule mlst_bm_get_mlst_mem:
    input:
        bm_asm='local/mlst/results/{accession}/assemble/bm/spades.trace',
        bm_mlst='local/mlst/results/{accession}/assemble/bm/mlst_asm.trace'
    output:
        tab='local/mlst/results/{accession}/assemble/bm/trace.tab'
    run:

        pd.concat(
            [
                bm.MaxMemUsage(input.bm_asm).to_series('assemble'),
                bm.MaxMemUsage(input.bm_mlst).to_series('mlst')
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')

# mlst_bm_get_kesmlst_mem
#
# Merge runtimes for each step of the kesmlst pipeline into one table.
rule mlst_bm_get_kesmlst_mem:
    input:
        bm_kmertable='local/mlst/results/{accession}/kesmlst/bm/kmertable.trace',
        bm_kesmlst='local/mlst/results/{accession}/kesmlst/bm/kesmlst.trace'
    output:
        tab='local/mlst/results/{accession}/kesmlst/bm/trace.tab'
    run:

        pd.concat(
            [
                bm.MaxMemUsage(input.bm_kmertable).to_series('kmertable'),
                bm.MaxMemUsage(input.bm_kesmlst).to_series('kesmlst')
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')


#
# Runtimes
#

# mlst_bm_get_mlst_runtime
#
# Merge runtimes for each step of the assembly pipeline into one table.
rule mlst_bm_get_mlst_runtime:
    input:
        bm_asm='local/mlst/results/{accession}/assemble/bm/spades.time',
        bm_mlst='local/mlst/results/{accession}/assemble/bm/mlst_asm.time'
    output:
        tab='local/mlst/results/{accession}/assemble/bm/time.tab'
    run:

        pd.concat(
            [
                bm.RunTime(input.bm_asm).to_series('assemble'),
                bm.RunTime(input.bm_mlst).to_series('mlst'),
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')

# mlst_bm_get_kesmlst_runtime
#
# Merge runtimes for each step of the kesmlst pipeline into one table.
rule mlst_bm_get_kesmlst_runtime:
    input:
        bm_kmertable='local/mlst/results/{accession}/kesmlst/bm/kmertable.time',
        bm_kesmlst='local/mlst/results/{accession}/kesmlst/bm/kesmlst.time'
    output:
        tab='local/mlst/results/{accession}/kesmlst/bm/time.tab'
    run:

        pd.concat(
            [
                bm.RunTime(input.bm_kmertable).to_series('kmertable'),
                bm.RunTime(input.bm_kesmlst).to_series('kesmlst'),
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='time')
