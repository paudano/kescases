"""
Rules to generate plots and summaries for the Strep project
"""

import os
import pandas as pd
import numpy as np

from kescaseslib import bm


#############
### Rules ###
#############

# strep_summary_bam_vs_ikc_file_size
#
# Summarize file sizes by number of reads.
rule strep_summary_bam_vs_ikc_file_size:
    input:
        tab='local/strep/summary/benchmarks/seq_file_size.tab'
    output:
        pdf='local/strep/summary/plots/size/size_ikc_bam_by_reads.pdf',
        eps='local/strep/summary/plots/size/size_ikc_bam_by_reads.eps',
        rdata='local/strep/summary/plots/size/size_ikc_bam_by_reads.RData'
    shell:
        """Rscript scripts/plots/strep_file_reads_size_plot.R {input.tab} local/strep/summary/plots/size/size_ikc_bam_by_reads"""

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
# Variant Call Depth
#

# strep_variant_call_depth
#
# Plot depth of each variant as annotated by the variant caller vs the alignment depth at the variant locus.
rule strep_variant_call_depth:
    input:
        var_kes='local/strep/summary/{pipeline}/variants.tab'
    output:
        pdf_kes='local/strep/summary/plots/variant/variant_depth_{pipeline}.pdf'
    shell:
        """Rscript scripts/plots/var_vs_align_depth.R {input.var_kes} {output.pdf_kes} {wildcards.pipeline}"""

#
# Benchmark plots
#

# strep_summary_bm_runtime_box_plot_by_reads
#
# Make runtime boxplots with CPU time divided by millions of reads.
rule strep_summary_bm_runtime_box_plot_by_reads:
    input:
        tab_kes='local/strep/temp/benchmarks/kestrel_runtime_withreads.tab',
        tab_gatk='local/strep/temp/benchmarks/gatk_runtime_withreads.tab',
        tab_asm='local/strep/temp/benchmarks/assemble_runtime_withreads.tab'
    output:
        pdf_all='local/strep/summary/plots/bm/rt/runtime_cpu_by_reads.pdf',
        pdf_noasm='local/strep/summary/plots/bm/rt/runtime_cpu_noasm_by_reads.pdf',
        pdf_all_multi='local/strep/summary/plots/bm/rt/multipart/runtime_cpu_by_reads.pdf',
        pdf_noasm_multi='local/strep/summary/plots/bm/rt/multipart/runtime_cpu_noasm_by_reads.pdf'
    shell:
        """Rscript scripts/plots/runtime_box_plot_by_reads.R {input.tab_kes} {input.tab_gatk} {input.tab_asm} {output.pdf_all} {output.pdf_noasm} {output.pdf_all_multi} {output.pdf_noasm_multi}"""

# strep_summary_runtime_with_reads
#
# Merge runtime and file size tables.
rule strep_summary_runtime_with_reads:
    input:
        tab_rt='local/strep/summary/benchmarks/{pipeline}_runtime.tab',
        tab_file='local/strep/summary/benchmarks/seq_file_size.tab'
    output:
        tab='local/strep/temp/benchmarks/{pipeline}_runtime_withreads.tab'
    run:

        # Read
        pd.concat(
            [
                pd.read_table(input.tab_rt, header=0, index_col=0),
                pd.read_table(input.tab_file, header=0, index_col=0)
            ],
            axis=1
        ).to_csv(output.tab, sep='\t', index=True, index_label='accession')

# strep_summary_bm_runtime_box_plot
#
# Make a boxplot of run times.
rule strep_summary_bm_runtime_box_plot:
    input:
        tab_kes='local/strep/summary/benchmarks/kestrel_runtime.tab',
        tab_gatk='local/strep/summary/benchmarks/gatk_runtime.tab',
        tab_asm='local/strep/summary/benchmarks/assemble_runtime.tab'
    output:
        pdf_cpu_all='local/strep/summary/plots/bm/rt/runtime_cpu.pdf',
        pdf_cpu_noasm='local/strep/summary/plots/bm/rt/runtime_cpu_noasm.pdf',
    shell:
        """Rscript scripts/plots/runtime_box_plot.R {input.tab_kes} {input.tab_gatk} {input.tab_asm} {output.pdf_cpu_all} {output.pdf_cpu_noasm}"""

# strep_summary_bm_mem_box_plot
#
# Make a boxplot of memory usage.
rule strep_summary_bm_mem_box_plot:
    input:
        tab_kes='local/strep/summary/benchmarks/kestrel_trace.tab',
        tab_gatk='local/strep/summary/benchmarks/gatk_trace.tab',
        tab_asm='local/strep/summary/benchmarks/assemble_trace.tab'
    output:
        pdf='local/strep/summary/plots/bm/mem/memory_trace_bm.pdf',
        pdf_multi='local/strep/summary/plots/bm/mem/multipart/memory_trace_bm.pdf'
    shell:
        """Rscript scripts/plots/memory_box_plot.R {input.tab_kes} {input.tab_gatk} {input.tab_asm} {output.pdf} {output.pdf_multi}"""

# strep_summary_bam_vs_ikc_size
#
# BAM vs IKC size plot.
rule strep_summary_bam_vs_ikc_size:
    input:
        size_tab='local/strep/summary/benchmarks/seq_file_size.tab',
        blacklist='data/strep/NC_003028.blacklist.tab'
    output:
        pdf='local/strep/summary/plots/size/size_bam_vs_ikc.pdf',
        pdf_multi='local/strep/summary/plots/size/multipart/size_bam_vs_ikc.pdf',
    shell:
        """Rscript scripts/plots/bam_vs_ikc_file_size.R {input.size_tab} {input.blacklist} {output.pdf} {output.pdf_multi}"""


#
# Benchmark tables
#

# strep_summary_filesize
#
# Get file size summary.
rule strep_summary_filesize:
    input:
        tab='local/strep/summary/benchmarks/seq_file_size.tab'
    output:
        tab='local/strep/summary/benchmarks/seq_file_size_summary.tab'
    run:

        # Read table
        df = pd.read_table(input.tab, header=0, index_col=0)

        # Get stats
        stat_mean = df.apply(np.mean, axis=0)
        stat_mean.name = 'MEAN'

        stat_median = df.apply(np.mean, axis=0)
        stat_median.name = 'MED'

        stat_sd = df.apply(np.std, axis=0)
        stat_sd.name = 'SD'

        stat_min = df.apply(np.min, axis=0)
        stat_min.name = 'MIN'

        stat_max = df.apply(np.max, axis=0)
        stat_max.name = 'MAX'

        # Merge stats
        df_stat = pd.concat(
            [stat_mean, stat_median, stat_sd, stat_min, stat_max],
            axis=1
        )

        # Write
        df_stat.to_csv(output.tab, sep='\t', index=True, index_label='STAT', float_format='%.2f')

# strep_summary_runtime_by_reads
#
# Get runtime performance by number of reads.
rule strep_summary_runtime_by_reads:
    input:
        tab_rt='local/strep/summary/benchmarks/{pipeline}_runtime.tab',
        tab_reads='local/strep/samples/sample_seq_summary.tab'
    output:
        tab='local/strep/summary/benchmarks/{pipeline}_runtime_summary_byreads.tab'
    run:

        # Read
        rt = pd.read_table(input.tab_rt, header=0, index_col=0)['user']
        size = pd.read_table(input.tab_reads, header=0, index_col=0)

        # Get stats
        rt_reads = rt / size['reads'] * 1e6
        rt_reads.name = 'byread'

        rt_bases = rt / size['bases'] * 1e6
        rt_bases.name = 'bybase'

        # Merge
        df_sample = pd.concat(
            [rt_reads, rt_bases],
            axis=1
        )

        # Get stats
        stat_mean = df_sample.apply(np.mean, axis=0)
        stat_mean.name = 'MEAN'

        stat_median = df_sample.apply(np.median, axis=0)
        stat_median.name = 'MED'

        stat_sd = df_sample.apply(np.std, axis=0)
        stat_sd.name = 'SD'

        # Merge stats
        df_stat = pd.concat(
            [stat_mean, stat_median, stat_sd],
            axis=1
        )

        # Write
        df_stat.to_csv(output.tab, sep='\t', index=True, index_label='STAT', float_format='%.2f')

# strep_summary_summarize_bm
#
# Summarize benchmark statistics (trace and runtime).
rule strep_summary_summarize_runtime:
    input:
        tab='local/strep/summary/benchmarks/{pipeline}_{stat}.tab'
    output:
        tab='local/strep/summary/benchmarks/{pipeline}_{stat}_summary.tab'
    run:

        # Read
        df = pd.read_table(input.tab, header=0, index_col=0)

        # Get stats
        stat_mean = df.apply(np.mean, axis=0)
        stat_mean.name = 'MEAN'

        stat_median = df.apply(np.median, axis=0)
        stat_median.name = 'MED'

        stat_sd = df.apply(np.std, axis=0)
        stat_sd.name = 'SD'

        stat_min = df.apply(np.min, axis=0)
        stat_min.name = 'MIN'

        stat_max = df.apply(np.max, axis=0)
        stat_max.name = 'MAX'

        # Merge
        df_stat = pd.concat(
            [stat_mean, stat_median, stat_sd, stat_min, stat_max],
            axis=1
        )

        # Write
        df_stat.to_csv(output.tab, sep='\t', index=True, index_label='STAT', float_format='%.2f')


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
                bm.summarize_trace_table(
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
                bm.summarize_runtime_table(
                    'local/strep/results/{}/{}/bm/time.tab'.format(accession, wildcards.pipeline), accession
                ) for accession in STREP_ACCESSIONS
            ],
            axis=1
        ).transpose().to_csv(output.tab, sep='\t', index=True, index_label='accession', float_format='%.2f')

# strep_summary_merge_ikc_disk_usage
#
# Get IKC segment file count and disk usage (sizes in KB) along with the number
# of sequence reads and bases for each sample.
rule strep_summary_merge_ikc_disk_usage:
    input:
        seg_size=expand('local/strep/results/{accession}/kestrel/bm/seg_size', accession=STREP_ACCESSIONS),
        seg_count=expand('local/strep/results/{accession}/kestrel/bm/seg_count', accession=STREP_ACCESSIONS),
        ikc=expand('local/strep/results/{accession}/kestrel/kmertable.ikc', accession=STREP_ACCESSIONS),
        bam=expand('local/strep/results/{accession}/gatk/sample.bam', accession=STREP_ACCESSIONS),
        seq_tab='local/strep/samples/sample_seq_summary.tab'
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

        # Read sequence data
        df_seq = pd.read_table(input.seq_tab, header=0, index_col=0)

        # Merge and write
        df = pd.concat(ikc_stat_list, axis=1).transpose()
        df = pd.concat([df, df_seq], axis=1)
        df = df.ix[:, ('ikc_size', 'bam_size', 'seg_size', 'seg_n', 'reads', 'bases')]

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
        pdf='local/strep/summary/plots/varjitter/{pipeline}_variant_jitter_all.pdf',
        pdf_multi='local/strep/summary/plots/varjitter/multipart/{pipeline}_variant_jitter_all.pdf'
    shell:
        """Rscript scripts/plots/variant_jitter.R {wildcards.pipeline} {input.tab} {output.pdf} {output.pdf_multi} {wildcards.pipeline}"""

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
