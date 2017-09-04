"""
Rules to summarize data.
"""

import pandas as pd
import numpy as np


#
# Plot
#

# ecoli_summary_plot_consensus_size
#
# Make consensus size distribution histogram.
rule ecoli_summary_plot_consensus_size:
    input:
        tab_con='local/ecoli/summary/kestrel/consensus_len.tab',
        tab_hap='local/ecoli/summary/kestrel/haplotype_len.tab'
    output:
        pdf='local/ecoli/summary/plots/con_size_hist.pdf'
    shell:
        """Rscript scripts/plots/con_size_hist.R {input.tab_con} {input.tab_hap} {output.pdf}"""


#
# Summary tables
#

# ecoli_summary_summary_size_table
#
# Get average consensus and haplotype sizes.
rule ecoli_summary_summary_size_table:
    input:
        tab_con='local/ecoli/summary/kestrel/consensus_len.tab',
        tab_hap='local/ecoli/summary/kestrel/haplotype_len.tab'
    output:
        tab='local/ecoli/summary/kestrel/summary_len.tab'
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
        tab_con='local/ecoli/summary/{pipeline}/call_stats_con.tab',
        tab_hap='local/ecoli/summary/{pipeline}/call_stats_hap.tab'
    output:
        tab='local/ecoli/summary/{pipeline,kestrel|gatk}/call_stats_summary.tab'
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
        call='local/ecoli/temp/{pipeline}/call_stats_{varset}.tab',
        con_size='local/ecoli/summary/kestrel/consensus_len.tab'
    output:
        tab='local/ecoli/summary/{pipeline,kestrel|gatk}/call_stats_{varset,con|hap}.tab'
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
        tab=expand('local/ecoli/results/{accession}/{{pipeline}}/variants_{{varset}}.tab', accession=ECOLI_ACCESSIONS)
    output:
        tab=temp('local/ecoli/temp/{pipeline,kestrel|gatk}/call_stats_{varset,con|hap}.tab')
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
                'local/ecoli/results/{}/{}/variants_{}.tab'.format(accession, wildcards.pipeline, wildcards.varset),
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
        tab='local/ecoli/summary/kestrel/haplotype_len.tab'
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
        tab='local/ecoli/summary/kestrel/consensus_len.tab'
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


#
# Kestrel variant call depth
#

# ecoli_summary_tp_call_graph
#
# Kestrel-only BED graph.
rule ecoli_summary_tp_call_graph:
    input:
        bed='local/ecoli/temp/summary/bed/variant_depth_{pipeline}.bed',
        fai=ECOLI_REF_FAI
    output:
        bed='local/ecoli/summary/bed/variant_depth_{pipeline,kestrel|gatk}.bed'
    shell:
        """echo -e "#CHROM\tPOS\tEND\tDEPTH" >{output.bed}; """
        """bedtools genomecov -bg -i {input.bed} -g {input.fai} """
        """>>{output.bed}"""

# ecoli_summary_tp_call_depth
#
# Merge BED records from all samples into one file.
rule ecoli_summary_tp_call_depth:
    input:
        bed=expand('local/ecoli/results/{accession}/{{pipeline}}/bed/variants_con_tp.bed', accession=ECOLI_ACCESSIONS)
    output:
        bed=temp('local/ecoli/temp/summary/bed/variant_depth_{pipeline,kestrel|gatk}.bed')
    shell:
        """sort -k1,1 -k2,2n -m {input.bed} """
        """>{output.bed}"""

#
# Merged variant calls
#

# ecoli_summary_merge_variants_in_region
#
# Merge variants in one region.
rule ecoli_summary_merge_variants_in_region:
    input:
        bed=expand('local/ecoli/results/{accession}/{{pipeline}}/bed/variants_con_{{call_type}}.bed', accession=ECOLI_ACCESSIONS)
    output:
        bed='local/ecoli/summary/bed/regions/region_{pipeline}_{call_type}_{pos}_{end}.bed'
    run:

        loc_pos = int(wildcards.pos)
        loc_end = int(wildcards.end)

        # Read each sample, subset, and add to list
        df_list = list()

        for accession in ECOLI_ACCESSIONS:
            df = pd.read_table(
                'local/ecoli/results/{}/{}/bed/variants_con_tp.bed'.format(
                    accession, wildcards.pipeline, wildcards.call_type
                ),
                header=0
            )

            df['ACCESSION'] = accession

            df = df.loc[df.apply(lambda row: row['END'] >= loc_pos and row['POS'] <= loc_end, axis=1)]

            df_list.append(df)

        # Merge
        df = pd.concat(df_list, axis=0)

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        df.to_csv(output.bed, sep='\t', index=False)


#
# Kestrel-only coverage BED
#

# ecoli_summary_kes_only_bed_graph
#
# Kestrel-only BED graph.
rule ecoli_summary_kes_only_bed_graph:
    input:
        bed='local/ecoli/temp/summary/bed/kes_only_merged.bed',
        fai=ECOLI_REF_FAI
    output:
        bed='local/ecoli/summary/bed/kes_only_merged.bed'
    shell:
        """echo -e "#CHROM\tPOS\tEND\tLENGTH\tDEPTH" >{output.bed}; """
        """bedtools genomecov -bg -i {input.bed} -g {input.fai} | """
        """awk -vOFS="\t" '{{print $1, $2, $3, $3 - $2, $4}}' """
        """>>{output.bed}"""

# ecoli_summary_bed_merge_kes_only
#
# Merge BED records from all samples into one file.
rule ecoli_summary_bed_merge_kes_only:
    input:
        bed=expand('local/ecoli/results/{accession}/kestrel/bed/regions_no_gatk.bed', accession=ECOLI_ACCESSIONS)
    output:
        bed=temp('local/ecoli/temp/summary/bed/kes_only_merged.bed')
    shell:
        """sort -k1,1 -k2,2n -m {input.bed} """
        """>{output.bed}"""
