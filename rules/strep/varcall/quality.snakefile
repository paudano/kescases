"""
Variant call quality, metrics, and summary statistics
"""

import pandas as pd

from kescaseslib import interval


#################
### Functions ###
#################

def _merge_variant_data_frames(input):
    """
    :param input: An input object with `input.tp`, `input.fp`, and `input.fn` pointing to the
        table files extracted from the VCF files (rule `strep_variant_vcf_eval_to_table`), and
        `input.pbp` pointing to the PBP BED file where genes are listed.

    :return: A merged and sorted dataframe with annotated regions.
    """

    # Get first table (save column order)
    df_tp = pd.read_table(input.tp, header=0)

    df = df_tp

    # Merge FP
    df_fp = pd.read_table(input.fp, header=0)

    if df_fp.shape[0] > 0:
        df = pd.concat([df, df_fp])

    # Merge FN
    df_fn = pd.read_table(input.fn, header=0)

    if df_fn.shape[0] > 0:
        df = pd.concat([df, df_fn])

    # Rearrange columns
    df = df.ix[:, df_tp.columns.tolist()]

    # Sort
    df = df.sort_values(['CHROM', 'POS'])

    # Annotate by PBP region
    pbp_interval = interval.IntervalContainer()
    pbp_interval.add_bed(input.pbp)

    df['REGION'] = df.apply(lambda row: str(pbp_interval.get_interval(row['CHROM'], row['POS'])), axis=1)

    # Return dataframe
    return df

def _get_depth_column_from_pileup(df, pileup_file_name):
    """
    Get a columns for `df` that contains the alignment depth at each locus.

    :param df: Merged dataframe.
    :param pileup_file_name: Name of the pileup file to extract read depths from.

    :return: A Series that can be added to `df`.
    """
    df_pileup = pd.read_table(pileup_file_name, header=None)

    depth_df = pd.Series(
        list(df_pileup[3]),
        index=pd.MultiIndex.from_tuples(list(zip(df_pileup[0], df_pileup[1])), names=('chrom', 'pos'))
    )

    return df.apply(lambda row: depth_df[(row['CHROM'], row['POS'])] if (row['CHROM'], row['POS']) in depth_df else 0, axis=1)


#############
### Rules ###
#############

# strep_variant_merge_kestrel_vcf_eval
#
# Merge Kestrel TP, FP, and FN calls annotated by vcfeval.
rule strep_variant_merge_kestrel_vcf_eval:
    input:
        tp='local/strep/temp/{accession}/kestrel/tp.tab',
        fp='local/strep/temp/{accession}/kestrel/fp.tab',
        fn='local/strep/temp/{accession}/kestrel/fn.tab',
        pbp=STREP_PBP_BED
    output:
        tab='local/strep/results/{accession}/kestrel/variants.tab'
    run:

        # Get merged variants
        df = _merge_variant_data_frames(input)

        # Write table
        df.to_csv(output.tab, sep='\t', index=False)

# strep_variant_merge_gatk_vcf_eval
#
# Merge GATK TP, FP, and FN calls annotated by vcfeval.
rule strep_variant_merge_gatk_vcf_eval:
    input:
        tp='local/strep/temp/{accession}/gatk/tp.tab',
        fp='local/strep/temp/{accession}/gatk/fp.tab',
        fn='local/strep/temp/{accession}/gatk/fn.tab',
        pileup='local/strep/results/{accession}/gatk/pileup.tab',
        pbp=STREP_PBP_BED
    output:
        tab='local/strep/results/{accession}/gatk/variants.tab'
    run:

        # Get merged variants
        df = _merge_variant_data_frames(input)

        # Annotate alignment depth
        df['DEPTH'] = _get_depth_column_from_pileup(df, input.pileup)

        # Write table
        df.to_csv(output.tab, sep='\t', index=False)


# strep_variant_vcf_eval_to_table
#
# Convert compared variants to a table file
rule strep_variant_vcf_eval_to_table:
    input:
        vcf='local/strep/results/{accession}/{pipeline}/vcfeval/{call}.vcf.gz',
    output:
        tab=temp('local/strep/temp/{accession}/{pipeline}/{call}.tab')
    run:

        call_type = wildcards.call.upper()

        shell(
            """zcat {input.vcf} | """
            """fgrep -v '##INFO' | """
            """awk -v OFS="\t" '{{$8 = "."; print}}' | """
            """bin/vcf2tsv -g | """
            """awk -vOFS="\t" -vCALLTYPE={call_type} '(NR == 1) {{$(NF + 1) = "CALL"; print}} (NR > 1) {{$(NF + 1) = CALLTYPE; print}}' """
#            """xargs -I LINE echo -e "LINE\t{call_type}" """
            """> {output.tab}"""
        )

# strep_variant_vcfeval
#
# Compare variants in GATK or Kestrel calls to the assembly calls.
rule strep_variant_vcfeval:
    input:
        vcf='local/strep/results/{accession}/{pipeline,kestrel|gatk}/variants.vcf.gz',
        base_vcf='local/strep/results/{accession}/assemble/variants.vcf.gz',
        ref=STREP_REF,
        pbp=STREP_PBP_BED,
        rtg_flag=STREP_RTG_INDEX_FLAG
    output:
        tp='local/strep/results/{accession}/{pipeline}/vcfeval/tp.vcf.gz',
        fp='local/strep/results/{accession}/{pipeline}/vcfeval/fp.vcf.gz',
        fn='local/strep/results/{accession}/{pipeline}/vcfeval/fn.vcf.gz',
        baseline='local/strep/results/{accession}/{pipeline}/vcfeval/tp-baseline.vcf.gz',
    log:
        'local/strep/results/{accession}/{pipeline}/log/strep_variant_vcfeval.log'
    shell:
        """rm -rf $(dirname {output.tp}); """
        """bin/rtg vcfeval """
            """-t $(dirname {input.rtg_flag}) """
            """-b {input.base_vcf} """
            """-c {input.vcf} """
            """-o $(dirname {output.tp}) """
            """--bed-regions {input.pbp} """
            """>{log} 2>&1"""
