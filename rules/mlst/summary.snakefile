"""
Summary statistics and tables for MLST.
"""

###################
### Definitions ###
###################

def _mlst_get_all_sample_file_names(wildcards):

    input_file_names = list()

    for accession in MLST_ACCESSIONS:
        input_file_names.append('local/mlst/samples/{}/{}_1.fastq.gz'.format(accession, accession))
        input_file_names.append('local/mlst/samples/{}/{}_2.fastq.gz'.format(accession, accession))

    return input_file_names


#############
### Rules ###
#############

# mlst_summary_runtime_by_reads
#
# Get runtime performance by number of reads.
rule mlst_summary_runtime_by_reads:
    input:
        tab_rt='local/mlst/summary/benchmarks/{pipeline}_runtime.tab',
        tab_reads='local/mlst/samples/sample_seq_summary.tab'
    output:
        tab='local/mlst/summary/benchmarks/{pipeline}_runtime_summary_byreads.tab'
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

# mlst_summary_merge_runtime_tables
#
# Merge runtime benchmarks for each sample.
rule mlst_summary_merge_runtime_tables:
    input:
        tab=expand('local/mlst/results/{accession}/{{pipeline}}/bm/time.tab', accession=MLST_ACCESSIONS)
    output:
        tab='local/mlst/summary/benchmarks/{pipeline}_runtime.tab'
    run:

        pd.concat(
            [
                bm.summarize_runtime_table(
                    'local/mlst/results/{}/{}/bm/time.tab'.format(accession, wildcards.pipeline), accession
                ) for accession in MLST_ACCESSIONS
            ],
            axis=1
        ).transpose().to_csv(output.tab, sep='\t', index=True, index_label='accession', float_format='%.2f')

# mlst_summary_summarize_bm
#
# Summarize benchmark statistics (trace and runtime).
rule mlst_summary_summarize_bm:
    input:
        tab='local/mlst/summary/benchmarks/{pipeline}_{stat}.tab'
    output:
        tab='local/mlst/summary/benchmarks/{pipeline}_{stat}_summary.tab'
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

# mlst_count_reads_and_bases
#
# Count the number of reads and bases in each sample.
rule mlst_count_reads_and_bases:
    input:
        fq=_mlst_get_all_sample_file_names
    output:
        tab='local/mlst/samples/sample_seq_summary.tab'
    run:

        # Create an empty data frame
        df = pd.DataFrame([], columns=['reads', 'bases'], dtype=(np.int64, np.int64))

        # Read sequence data for each accession
        for accession in MLST_ACCESSIONS:
            record_count = 0
            base_count = 0

            # Read base and record count
            for in_file_name in ['local/mlst/samples/{}/{}_{}.fastq.gz'.format(accession, accession, n) for n in range(1, 3)]:
                with io.TextIOWrapper(gzip.open(in_file_name)) as fq_file:
                    for seq_record in SeqIO.parse(fq_file, 'fastq'):
                        record_count += 1
                        base_count += len(seq_record.seq)

            # Add to table
            ser = pd.Series({'reads': record_count, 'bases': base_count}, dtype=(np.int64, np.int64))
            ser.name = accession

            df = df.append(ser)

        # Make integers
        df = df.astype(np.int64)

        # Write table
        df.to_csv(output.tab, sep='\t', index=True, index_label='accession')
