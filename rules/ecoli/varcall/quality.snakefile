"""
Variant call quality, metrics, and summary statistics
"""

import pandas as pd
import numpy as np

from kescaseslib import interval
from kescaseslib import kesutil
from kescaseslib import variant


#################
### Functions ###
#################

def _merge_variant_data_frames_ecoli(input):
    """
    :param input: An input object with `input.tp`, `input.fp`, and `input.fn` pointing to the
        table files extracted from the VCF files (rule `ecoli_variant_vcf_eval_to_table`), and
        `input.ref_bed` pointing to the BED file with a record covering the reference genome.

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

    # Return dataframe
    return df

def _get_depth_column_from_pileup_ecoli(df, pileup_file_name):
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

def _set_filter_and_id_ecoli(df, filter_container):
    """
    Set the FILTER and ID fields in a dataframe of variants.

    :param df: Dataframe of all variants. Must contain fields "CHROM", "POS", "REF", and "ALT".
    :param filter_container: Interval container with filtered regions loaded. The `tag` field of each
        interval in the container is used to annotated the FILTER column for variants within it. Filters
        are processed in the order they are read.

    :return: Updated data frame.
    """

    if df.shape[0] == 0:
        return df

    # Get list of variants
    var_list = df.apply(
        lambda row: variant.Variant(row['CHROM'], row['POS'], row['REF'], row['ALT'], None, row['SAMPLE']),
        axis=1
    )

    # Set ID
    df['ID'] = [str(variant) for variant in var_list]

    # Set FILTER
    df['FILTER'] = 'PASS'
#    df['FILTER'] = [
#        filter_region.tag if filter_region is not None else 'PASS' for filter_region in [
#            filter_container.get_interval(variant.chrom, variant.start, variant.get_end()) for variant in var_list
#        ]
#    ]

    return df


#############
### Rules ###
#############


#
# Kestrel-only coverage BED
#

# ecoli_bed_intersect_variant_regions
#
# Get variants in Kestrel-only regions.
rule ecoli_bed_intersect_variant_regions:
    input:
        variants='local/ecoli/results/{accession}/kestrel/bed/variants_con_tp.bed',
        regions='local/ecoli/results/{accession}/kestrel/bed/regions_no_gatk.bed'
    output:
        bed='local/ecoli/results/{accession}/kestrel/bed/variants_no_gatk.bed',
    shell:
        """bedtools intersect -a {input.variants} -b {input.regions} -wa -u """
        """>{output.bed}"""

# ecoli_bed_intersect_kes_no_gatk
#
# Find variant in regions GATK did not call.
rule ecoli_bed_regions_kes_no_gatk:
    input:
        kes='local/ecoli/temp/{accession}/kestrel/bed/variants_con_tp_50.bed',
        gatk='local/ecoli/temp/{accession}/gatk/bed/variants_con_tp_50.bed',
        fai=ECOLI_REF_FAI
    output:
        kes='local/ecoli/results/{accession}/kestrel/bed/regions_no_gatk.bed'
    shell:
        """bedtools intersect -a {input.kes} -b {input.gatk} -wa -v -header -sorted | """
        """bedtools slop -i stdin -g {input.fai} -b -50 | """
        """awk -vOFS="\\t" '$3 > $2' """
        """>{output.kes}"""

# ecoli_bed_pad_regions
#
# Pad BED regions.
rule ecoli_bed_pad_regions:
    input:
        bed='local/ecoli/results/{accession}/{pipeline}/bed/variants_{varset}_tp.bed',
        fai=ECOLI_REF_FAI
    output:
        bed=temp('local/ecoli/temp/{accession}/{pipeline,kestrel|gatk}/bed/variants_{varset}_tp_{pad}.bed')
    shell:
        """bedtools slop -i {input.bed} -g {input.fai} -b {wildcards.pad} -header | """
        """bedtools merge """
        """>{output.bed}"""

# ecoli_bed_variant_to_bed
#
# Make BED of variant calls.
rule ecoli_bed_variant_to_bed:
    input:
        tab='local/ecoli/results/{accession}/{pipeline}/variants_{varset}.tab'
    output:
        bed='local/ecoli/results/{accession}/{pipeline,kestrel|gatk}/bed/variants_{varset}_{call_type}.bed'
    run:

        call_type = wildcards.call_type.upper()

        line_count = 0

        with open(input.tab, 'r') as in_file:
            with open(output.bed, 'w') as out_file:

                # Discard first line
                line = next(in_file)

                # Write header
                out_file.write('#CHROM\tPOS\tEND\tID\tTYPE\tLEN\tCALL\n')

                # Read each line
                for line in in_file:
                    line_count += 1

                    line = line.strip()

                    if not line:
                        continue

                    # Tokenize
                    tok = line.split('\t')

                    if len(tok) not in {12, 14}:
                        raise ValueError('Expected 12 or 14 fields, received {} (line={})'.format(len(tok), line_count))

                    # Only TP calls
                    if tok[-1] != call_type and call_type != 'ALL':
                        continue

                    # Get attributes
                    var_chrom = tok[0]
                    var_pos = int(tok[1]) - 1
                    var_id = tok[2]
                    var_len_ref = len(tok[3])
                    var_len_alt = len(tok[4])
                    var_call = tok[-1]

                    # Check
                    if var_len_ref != 1 and var_len_alt != 1:
                        raise ValueError('REF or ALT length must be 1: REF={}, ALT={} (line={})'.format(tok[3], tok[4], line_count))

                    # Get attributes
                    if var_len_ref == 1 and var_len_alt == 1:

                        # SNP
                        var_type = 'SNP'
                        var_len = 1
                        var_end = var_pos + 1

                    elif var_len_ref > 1:

                        # DEL
                        var_type = 'DEL'
                        var_len = var_len_ref
                        var_end = var_pos + var_len_ref

                    else:

                        # INS
                        var_type = 'INS'
                        var_len = var_len_alt
                        var_end = var_pos + 1

                    # Write
                    out_file.write('\t'.join((var_chrom, str(var_pos), str(var_end), var_id, var_type, str(var_len), var_call)))
                    out_file.write('\n')

#
# Annotate calls
#

# ecoli_annotate_variant_calls
#
# Add ID and FILTER columns.
rule ecoli_annotate_variant_calls:
    input:
        tab='local/ecoli/temp/{accession}/{pipeline}/vcfeval/{varset}/variants_unannotated.tab',
        nc_bed='local/ecoli/results/{accession}/assemble/no_consensus.bed'
    output:
        tab='local/ecoli/results/{accession}/{pipeline,kestrel|gatk}/variants_{varset}.tab'
    run:

        # Read variants and the blacklist
        df = pd.read_table(input.tab, header=0)

        # Initialize filters
        filter_container = interval.IntervalContainer()

        #filter_container.add_bed(input.nc_bed, 'NO_CONSENSUS')
        #filter_container.add_blacklist(input.bl_tab, wildcards.accession)

        # Annotate filter
        df = _set_filter_and_id_ecoli(df, filter_container)

        # Write
        df.to_csv(output.tab, sep='\t', index=False)


#
# VCFEVAL
#

# ecoli_variant_merge_kestrel_vcf_eval
#
# Merge Kestrel TP, FP, and FN calls annotated by vcfeval. Add a filed for the alignment depth.
rule ecoli_variant_merge_kestrel_vcf_eval:
    input:
        tp='local/ecoli/temp/{accession}/{pipeline}/vcfeval/{varset}/tp.tab',
        fp='local/ecoli/temp/{accession}/{pipeline}/vcfeval/{varset}/fp.tab',
        fn='local/ecoli/temp/{accession}/{pipeline}/vcfeval/{varset}/fn.tab',
        ref_bed=ECOLI_REF_BED
    output:
        tab=temp('local/ecoli/temp/{accession}/{pipeline,kestrel|gatk}/vcfeval/{varset}/variants_unannotated.tab')
    run:

        # Get merged variants
        df = _merge_variant_data_frames_ecoli(input)

#        # Annotate alignment depth
#        if df.shape[0] > 0:
#            df['DEPTH'] = _get_depth_column_from_pileup_ecoli(df, input.pileup)
#        else:
#            df['DEPTH'] = pd.Series({}, dtype=np.int32)

        # Write table
        df.to_csv(output.tab, sep='\t', index=False)


# ecoli_variant_vcfeval_to_table
#
# Convert compared variants to a table file
rule ecoli_variant_vcfeval_to_table:
    input:
        vcf='local/ecoli/results/{accession}/{pipeline}/vcfeval/{varset}/{call}.vcf.gz',
    output:
        tab=temp('local/ecoli/temp/{accession}/{pipeline,kestrel|gatk}/vcfeval/{varset,con|hap}/{call,tp|fp|fn}.tab')
    run:

        call_type = wildcards.call.upper()

        shell(
            """zcat {input.vcf} | """
            """fgrep -v '##INFO' | """
            """awk -v OFS="\t" '{{$8 = "."; print}}' | """
            """bin/vcf2tsv -g | """
            """awk -vOFS="\t" -vCALLTYPE={call_type} '(NR == 1) {{$(NF + 1) = "CALL"; print}} (NR > 1) {{$(NF + 1) = CALLTYPE; print}}' """
            """> {output.tab}"""
        )

# ecoli_variant_vcfeval
#
# Compare variants in GATK or Kestrel calls to the assembly calls.
rule ecoli_variant_vcfeval:
    input:
        vcf='local/ecoli/results/{accession}/{pipeline}/variants.vcf.gz',
        base_vcf='local/ecoli/results/{accession}/assemble/variants.vcf.gz',
        ref=ECOLI_REF,
        rtg_flag=ECOLI_RTG_INDEX_FLAG,
        con_bed='local/ecoli/results/{accession}/kestrel/consensus_regions.bed',
        hap_bed='local/ecoli/results/{accession}/kestrel/haplotypes.bed'
    output:
        tp='local/ecoli/results/{accession}/{pipeline,kestrel|gatk}/vcfeval/{varset,con|hap}/tp.vcf.gz',
        fp='local/ecoli/results/{accession}/{pipeline,kestrel|gatk}/vcfeval/{varset,con|hap}/fp.vcf.gz',
        fn='local/ecoli/results/{accession}/{pipeline,kestrel|gatk}/vcfeval/{varset,con|hap}/fn.vcf.gz',
        bl='local/ecoli/results/{accession}/{pipeline,kestrel|gatk}/vcfeval/{varset,con|hap}/tp-baseline.vcf.gz'
    log:
        'local/ecoli/results/{accession}/{pipeline}/log/{pipeline}_variant_vcfeval_{varset}.log'
    run:

        # Get regions based varset wildcard
        if wildcards.varset == 'con':
            bed_regions = input.con_bed

        elif wildcards.varset == 'hap':
            bed_regions = input.hap_bed

        else:
            raise RuntimeError('Unrecognized variant set for vcfeval: {}'.format(wildcards.varset))

        # Run vcfeval
        if kesutil.has_uncommented_lines(input.base_vcf):
            shell(
                """rm -rf $(dirname {output.tp}); """
                """bin/rtg vcfeval """
                    """-t $(dirname {input.rtg_flag}) """
                    """-b {input.base_vcf} """
                    """-c {input.vcf} """
                    """--bed-regions={bed_regions} """
                    """-o $(dirname {output.tp}) """
                    """>{log} 2>&1"""
            )
        else:
            shell(
                """cp {input.base_vcf} {output.fn}; """
                """cp {input.base_vcf} {output.bl}; """
                """cp {input.vcf} {output.tp}; """
                """cp {input.vcf} {output.fp}"""
            )


#
# Get haplotype regions
#

# ecoli_variant_get_haplotype_regions
#
# Get regions where haplotypes were generated.
rule ecoli_variant_get_haplotype_regions:
    input:
        bam='local/ecoli/results/{accession}/kestrel/haplotypes.bam'
    output:
        bed='local/ecoli/results/{accession}/kestrel/haplotypes.bed'
    shell:
        """bamToBed -i {input.bam} | """
        """bedtools merge > {output.bed}"""
