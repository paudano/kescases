"""
Rules to generate plots and summaries for the Strep project
"""

import os
import pandas


###################
### Definitions ###
###################

def _summarize_runtime_table(tab_file_name, accession):
    """
    Generate a series by adding all the runtimes for each type (real, sys, user).

    :param tab_file_name: Name of the table file with runtimes for each step.
    :param accession: Sample accession.
    """

    summary_series = pd.read_table(tab_file_name, header=0, index_col=0).apply(sum, axis=1)

    summary_series.name = accession

    return summary_series

def _summarize_trace_table(tab_file_name, accession):
    """
    Get a series describing the maximum vss and rss from a table of rss and vss for each step.

    :param tab_file_name: Table of rss and vss for each step.
    :param accession: Sample accession.
    """

    summary_series = pd.read_table(tab_file_name, header=0, index_col=0).apply(max, axis=1)

    summary_series.name = accession

    return summary_series


#############
### Rules ###
#############

#
# Phylogeny circos plot
#

# strep_phylogeny_variant
#
# Generate a circos plot of variants against phylogeny.
rule strep_phylogeny_variant:
    input:
        ani_ident='local/strep/ani/ani_tables/ANIm_percentage_identity.tab',
        variants='local/strep/summary/kestrel/variants.tab'
    output:
        pdf='local/strep/summary/plots/phylo/phylo_variant.pdf'
    shell:
        """Rscript scripts/plots/strep_phylo_variant.R {input.ani_ident} {input.variants} {output.pdf}"""

# strep_phylogeny_serotype
#
# Generate a circos plot of variants against serotypes.
rule strep_phylogeny_serotype:
    input:
        ani_ident='local/strep/ani/ani_tables/ANIm_percentage_identity.tab',
        sra='data/strep/SraRunTable.txt'
    output:
        pdf='local/strep/summary/plots/phylo/phylo_sero.pdf'
    shell:
        """Rscript scripts/plots/strep_phylo_sero.R {input.ani_ident} {input.sra} {output.pdf}"""


#
# Benchmark plots
#

# strep_summary_bm_runtime_box_plot
#
# Make a boxplot of run times.
rule strep_summary_bm_runtime_box_plot:
    input:
        tab_kes='local/strep/summary/benchmarks/kestrel_runtime.tab',
        tab_gatk='local/strep/summary/benchmarks/gatk_runtime.tab',
        tab_asm='local/strep/summary/benchmarks/assemble_runtime.tab'
    output:
        pdf_cpu_all='local/strep/summary/plots/bm/runtime_cpu.pdf',
        pdf_real_all='local/strep/summary/plots/bm/runtime_real.pdf',
        pdf_cpu_noasm='local/strep/summary/plots/bm/runtime_cpu_noasm.pdf',
        pdf_real_noasm='local/strep/summary/plots/bm/runtime_real_noasm.pdf'
    shell:
        """Rscript scripts/plots/runtime_box_plot.R {input.tab_kes} {input.tab_gatk} {input.tab_asm} local/strep/summary/plots/bm/runtime"""

# strep_summary_bm_mem_box_plot
#
# Make a boxplot of memory usage.
rule strep_summary_bm_mem_box_plot:
    input:
        tab_kes='local/strep/summary/benchmarks/kestrel_trace.tab',
        tab_gatk='local/strep/summary/benchmarks/gatk_trace.tab',
        tab_asm='local/strep/summary/benchmarks/assemble_trace.tab'
    output:
        pdf='local/strep/summary/plots/bm/memory_trace_bm.pdf',
        eps='local/strep/summary/plots/bm/memory_trace_bm.eps',
        rdata='local/strep/summary/plots/bm/memory_trace_bm.RData'
    shell:
        """Rscript scripts/plots/memory_box_plot.R {input.tab_kes} {input.tab_gatk} {input.tab_asm} local/strep/summary/plots/bm/memory_trace_bm"""

rule strep_summary_bam_vs_ikc_size:
    input:
        size_tab='local/strep/summary/benchmarks/seq_file_size.tab',
        blacklist='data/strep/NC_003028.blacklist.tab'
    output:
        pdf='local/strep/summary/plots/size/size_bam_vs_ikc.pdf',
        eps='local/strep/summary/plots/size/size_bam_vs_ikc.eps',
        rdata='local/strep/summary/plots/size/size_bam_vs_ikc.RData'
    shell:
        """Rscript scripts/plots/bam_vs_ikc_file_size.R {input.size_tab} {input.blacklist} local/strep/summary/plots/size/size_bam_vs_ikc"""

#
# Benchmark tables
#

# strep_summary_merge_trace_tables
#
# Merge runtime benchmarks for each sample.
rule strep_summary_merge_trace_tables:
    input:
        tab=expand('local/strep/results/{accession}/{{pipeline}}/bm/trace.tab', accession=STREP_ACCESSIONS)
    output:
        tab='local/strep/summary/benchmarks/{pipeline}_trace.tab'
    run:

        pd.concat(
            [
                _summarize_trace_table(
                    'local/strep/results/{}/{}/bm/trace.tab'.format(accession, wildcards.pipeline), accession
                ) for accession in STREP_ACCESSIONS
            ],
            axis=1
        ).transpose().to_csv(output.tab, sep='\t', index=True, index_label='accession', float_format='%.2f')

# strep_summary_merge_runtime_tables
#
# Merge runtime benchmarks for each sample.
rule strep_summary_merge_runtime_tables:
    input:
        tab=expand('local/strep/results/{accession}/{{pipeline}}/bm/time.tab', accession=STREP_ACCESSIONS)
    output:
        tab='local/strep/summary/benchmarks/{pipeline}_runtime.tab'
    run:

        pd.concat(
            [
                _summarize_runtime_table(
                    'local/strep/results/{}/{}/bm/time.tab'.format(accession, wildcards.pipeline), accession
                ) for accession in STREP_ACCESSIONS
            ],
            axis=1
        ).transpose().to_csv(output.tab, sep='\t', index=True, index_label='accession', float_format='%.2f')

# strep_summary_merge_ikc_disk_usage
#
# Get IKC segment file count and disk usage for each accession. Sizes are in KB.
rule strep_summary_merge_ikc_disk_usage:
    input:
        seg_size=expand('local/strep/results/{accession}/kestrel/bm/seg_size', accession=STREP_ACCESSIONS),
        seg_count=expand('local/strep/results/{accession}/kestrel/bm/seg_count', accession=STREP_ACCESSIONS),
        ikc=expand('local/strep/results/{accession}/kestrel/kmertable.ikc', accession=STREP_ACCESSIONS),
        bam=expand('local/strep/results/{accession}/gatk/sample.bam', accession=STREP_ACCESSIONS)
    output:
        tab='local/strep/summary/benchmarks/seq_file_size.tab'
    run:

        ikc_stat_list = list()

        for accession in STREP_ACCESSIONS:

            # Read sizes and counts
            with open('local/strep/results/{}/kestrel/bm/seg_size'.format(accession), 'r') as in_file:
                seg_size = int(next(in_file).strip())

            with open('local/strep/results/{}/kestrel/bm/seg_count'.format(accession), 'r') as in_file:
                seg_n = int(next(in_file).strip())

            ikc_size = int(os.path.getsize('local/strep/results/{}/kestrel/kmertable.ikc'.format(accession)) / 1024)
            bam_size = int(os.path.getsize('local/strep/results/{}/gatk/sample.bam'.format(accession)) / 1024)

            # Add series to the stat list
            ikc_stat_list.append(pd.Series(
                {'seg_size': seg_size, 'seg_n': seg_n, 'ikc_size': ikc_size, 'bam_size': bam_size}
            ))

            ikc_stat_list[-1].name = accession

        # Merge and write
        df = pd.concat(ikc_stat_list, axis=1).transpose()
        df = df.ix[:, ('ikc_size', 'bam_size', 'seg_size', 'seg_n')]

        df.to_csv(output.tab, sep='\t', index=True, index_label='accession')


#
# Jitter plot
#

# strep_summary_make_variant_jitter_plots
#
# Make jitter plot of the depth of all variants.
rule strep_summary_make_variant_jitter_plots:
    input:
        tab='local/strep/summary/{pipeline}/variants.tab'
    output:
        all_pdf='local/strep/summary/plots/varjitter/{pipeline}_variant_jitter_all.pdf',
        all_eps='local/strep/summary/plots/varjitter/{pipeline}_variant_jitter_all.eps',
        lc_pdf='local/strep/summary/plots/varjitter/{pipeline}_variant_jitter_lc.pdf',
        lc_eps='local/strep/summary/plots/varjitter/{pipeline}_variant_jitter_lc.eps',
        rdata='local/strep/summary/plots/varjitter/{pipeline}_variant_jitter.RData'
    shell:
        """Rscript scripts/plots/variant_jitter.R {wildcards.pipeline} {input.tab} local/strep/summary/plots/varjitter/{wildcards.pipeline}_variant_jitter"""

#
# Merge variant calls from all samples
#

# strep_get_variant_call_stats
#
# Get variant call summary statistics.
rule strep_get_variant_call_stats:
    input:
        tab='local/strep/summary/{pipeline}/variants.tab'
    output:
        tab='local/strep/summary/{pipeline}/quality/summary_stats.tab'
    run:

        # Read table and group by call evaluation
        df_merged = pd.read_table(input.tab, header=0, na_filter=False)
        gr_call_stat = df_merged.groupby('CALL')

        # Count classes
        tp = sum(df_merged['CALL'] == 'TP')
        fp = sum(df_merged['CALL'] == 'FP')
        fn = sum(df_merged['CALL'] == 'FN')

        # Create dataframe
        df_summary = pd.DataFrame(
            {
                'TP': (tp, ),
                'FP': (fp, ),
                'FN': (fn, ),
                'TPR': (tp / (tp + fn), ),
                'FDR': (fp / (tp + fp), ),
                'PPV': (tp / (tp + fp), )
            }
        )

        df_summary = df_summary.ix[:, ('TP', 'FP', 'FN', 'TPR', 'FDR', 'PPV')]

        # Write
        df_summary.to_csv(output.tab, sep='\t', index=False, float_format='%.2f')

# strep_merge_variant_calls
#
# Merge variant calls from all samples into one table.
rule strep_merge_variant_calls:
    input:
        tab=expand('local/strep/results/{accession}/{{pipeline}}/variants.tab', accession=STREP_ACCESSIONS)
    output:
        tab_pass='local/strep/summary/{pipeline}/variants.tab',
        tab_fail='local/strep/summary/{pipeline}/variants_removed.tab'
    run:

        df_merged = None

        # Read input files
        for input_file_name in input.tab:
            df_sample = pd.read_table(input_file_name, header=0, na_filter=False)

            # Ignore empty tables
            if df_sample.shape[0] == 0:
                continue

            if df_merged is not None:
                df_merged = df_merged.append(df_sample, ignore_index=True)
            else:
                df_merged = df_sample

        # Separate by filter
        df_pass = df_merged.ix[df_merged['FILTER'] == 'PASS', :]
        df_fail = df_merged.ix[df_merged['FILTER'] != 'PASS', :]

        # Write
        df_pass.to_csv(output.tab_pass, sep='\t', index=False)
        df_fail.to_csv(output.tab_fail, sep='\t', index=False)
