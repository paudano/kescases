"""
Rules to summarize data.
"""

import pandas as pd
import numpy as np

# ecoli_summary_plot_consensus_size
#
# Make consensus size distribution histogram.
rule ecoli_summary_plot_consensus_size:
    input:
        tab_con='local/ecoli/summary/table/consensus_len.tab',
        tab_hap='local/ecoli/summary/table/haplotype_len.tab'
    output:
        pdf='local/ecoli/summary/plots/con_size_hist.pdf'
    shell:
        """Rscript scripts/plots/con_size_hist.R {input.tab_con} {input.tab_hap} {output.pdf}"""


#
# Summary tables
#

rule ecoli_summary_summary_size_table:
    input:
        tab_con='local/ecoli/summary/table/consensus_len.tab',
        tab_hap='local/ecoli/summary/table/haplotype_len.tab'
    output:
        tab='local/ecoli/summary/table/summary_len.tab'
    run:

        # Read
        con = pd.read_table(input.tab_con, header=0, index_col=0, squeeze=True)
        con.name = 'Consensus'

        hap = pd.read_table(input.tab_hap, header=0, index_col=0, squeeze=True)
        hap.name = 'Haplotype'

        # Merge
        df = pd.concat(
            [con, hap],
            axis=1
        ).transpose()

        # Mean
        mean_len = df.apply(np.mean, axis=1)

        # Write
        mean_len.to_csv(output.tab, sep='\t', index=True, index_label='Type')



# ecoli_summary_make_summary_table
#
# Generate summary statistics for all calls (merge from samples).
rule ecoli_summary_make_summary_table:
    input:
        tab_con='local/ecoli/summary/table/call_stats_con.tab',
        tab_hap='local/ecoli/summary/table/call_stats_hap.tab'
    output:
        tab='local/ecoli/summary/table/call_stats_summary.tab'
    run:
        
        # Get tables
        df_con = pd.read_table(input.tab_con, header=0, index_col=0).ix[:, ('TP', 'FP', 'FN')]
        df_hap = pd.read_table(input.tab_hap, header=0, index_col=0).ix[:, ('TP', 'FP', 'FN')]

        # Sum counts
        sum_con = df_con.apply(sum, axis=0)
        sum_con.name = 'CONSENSUS'

        sum_hap = df_hap.apply(sum, axis=0)
        sum_hap.name = 'HAPLOTYPE'

        # Merge
        df = pd.concat(
            [sum_con, sum_hap],
            axis=1
        ).transpose()

        # Summary stats
        df['TPR'] = df.apply(lambda row: row['TP'] / (row['TP'] + row['FN']), axis=1)
        df['FDR'] = df.apply(lambda row: row['FP'] / (row['FP'] + row['TP']), axis=1)

        # Write
        df.to_csv(output.tab, sep='\t', index=True, index_label='SAMPLE')

# ecoli_summary_get_table
#
# Merge call stats with consensus regions.
rule ecoli_summary_get_table:
    input:
        call='local/ecoli/temp/table/call_stats_{varset}.tab',
        con_size='local/ecoli/summary/table/consensus_len.tab'
    output:
        tab='local/ecoli/summary/table/call_stats_{varset,con|hap}.tab'
    run:

        # Merge
        df = pd.concat(
            [
                pd.read_table(input.con_size, header=0, index_col=0),
                pd.read_table(input.call, header=0, index_col=0)
            ],
            axis=1
        )

        # Write
        df.to_csv(output.tab, sep='\t', header=True)


#
# Full tables
#

# ecoli_summary_get_call_stats
#
# Get call statistics for all samples.
rule ecoli_summary_get_call_stats:
    input:
        tab=expand('local/ecoli/results/{accession}/kestrel/variants_{{varset}}.tab', accession=ECOLI_ACCESSIONS)
    output:
        tab=temp('local/ecoli/summary/table/call_stats_{varset,con|hap}.tab')
    run:

        # Initialize counts
        n_tp = pd.Series(dtype=np.int32)
        n_tp.name = 'TP'

        n_fp = pd.Series(dtype=np.int32)
        n_fp.name = 'FP'

        n_fn = pd.Series(dtype=np.int32)
        n_fn.name = 'FN'

        # Get counts for each sample
        for accession in ECOLI_ACCESSIONS:
            call = pd.read_table(
                'local/ecoli/results/{}/kestrel/variants_{}.tab'.format(accession, wildcards.varset),
                header=0,
                usecols=['CALL'],
                squeeze=True
            )

            n_tp[accession] = sum(call == 'TP')
            n_fp[accession] = sum(call == 'FP')
            n_fn[accession] = sum(call == 'FN')

        # Get summary statistics
        tpr = n_tp / (n_tp + n_fn)
        tpr.name = 'TPR'

        fdr = n_fp / (n_fp + n_tp)
        fdr.name = 'FDR'

        # Merge and write
        pd.concat(
            [n_tp, n_fp, n_fn, tpr, fdr],
            axis=1
        ).to_csv(output.tab, sep='\t', header=True, index=True, index_label='SAMPLE')

# ecoli_summary_get_haplotype_size
#
# Get the size of all consensus regions for each sample.
rule ecoli_summary_get_haplotype_size:
    input:
        bed=expand('local/ecoli/results/{accession}/kestrel/haplotypes.bed', accession=ECOLI_ACCESSIONS)
    output:
        tab='local/ecoli/summary/table/haplotype_len.tab'
    run:

        # Initialize size
        con_size = pd.Series(dtype=np.int32)
        con_size.name = 'SIZE'

        # Get consensus size for each sample
        for accession in ECOLI_ACCESSIONS:
            df = pd.read_table('local/ecoli/results/{}/kestrel/haplotypes.bed'.format(accession), header=None)

            con_size[accession] = sum(df[2] - df[1])

        # Write
        con_size.to_csv(output.tab, sep='\t', header=True, index_label='SAMPLE')

# ecoli_summary_get_consensus_size
#
# Get the size of all consensus regions for each sample.
rule ecoli_summary_get_consensus_size:
    input:
        bed=expand('local/ecoli/results/{accession}/kestrel/consensus_regions.bed', accession=ECOLI_ACCESSIONS)
    output:
        tab='local/ecoli/summary/table/consensus_len.tab'
    run:

        # Initialize size
        con_size = pd.Series(dtype=np.int32)
        con_size.name = 'SIZE'

        # Get consensus size for each sample
        for accession in ECOLI_ACCESSIONS:
            df = pd.read_table('local/ecoli/results/{}/kestrel/consensus_regions.bed'.format(accession), header=None)

            con_size[accession] = sum(df[2] - df[1])

        # Write
        con_size.to_csv(output.tab, sep='\t', header=True, index_label='SAMPLE')
