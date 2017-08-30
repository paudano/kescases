### GATK Pipeline ###

# ecoli_kestrel_call_all
#
# Get calls for all samples.
rule ecoli_gatk_call_all:
    input:
        tab_con=expand('local/ecoli/results/{accession}/gatk/variants_con.tab', accession=ECOLI_ACCESSIONS),
        tab_hap=expand('local/ecoli/results/{accession}/gatk/variants_hap.tab', accession=ECOLI_ACCESSIONS)
    run:
        pass

# ecoli_gatk_call_variants
#
# Call variants.
rule ecoli_gatk_call_variants:
    input:
        bam='local/ecoli/results/{accession}/gatk/sample.bam',
        consensus_bed='local/ecoli/results/{accession}/kestrel/consensus_regions.bed'
    output:
        vcf=temp('local/ecoli/temp/{accession}/gatk/variants.vcf'),
        time='local/ecoli/results/{accession}/gatk/bm/variants.time',
        trace='local/ecoli/results/{accession}/gatk/bm/variants.trace'
    log:
        'local/ecoli/results/{accession}/gatk/log/ecoli_gatk_call_variants.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.gatk} """
            """-T HaplotypeCaller """
            """-L {input.consensus_bed} """
            """-R {ECOLI_REF} """
            """-I {input.bam} """
            """-o {output.vcf} """
            """-ploidy 1 """
            """>{log} 2>&1"""

# ecoli_gatk_indel_realign
#
# Perform indel-based realignment.
rule ecoli_gatk_indel_realign:
    input:
        bam='local/ecoli/temp/{accession}/gatk/align_marked.bam',
        bai='local/ecoli/temp/{accession}/gatk/align_marked.bai',
        intervals='local/ecoli/results/{accession}/gatk/realigner.intervals',
        ref=ECOLI_REF
    output:
        bam='local/ecoli/results/{accession}/gatk/sample.bam',
        time='local/ecoli/results/{accession}/gatk/bm/realign.time',
        trace='local/ecoli/results/{accession}/gatk/bm/realign.trace'
    log:
        'local/ecoli/results/{accession}/gatk/log/ecoli_gatk_indel_realign.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.gatk} """
            """-T IndelRealigner """
            """-R {input.ref} """
            """-targetIntervals {input.intervals} """
            """-I {input.bam} """
            """-o {output.bam}"""
            """>{log} 2>&1"""

# ecoli_gatk_realign_target_creator
#
# Indel realignment target creator.
rule ecoli_gatk_realign_target_creator:
    input:
        bam='local/ecoli/temp/{accession}/gatk/align_marked.bam',
        bai='local/ecoli/temp/{accession}/gatk/align_marked.bai',
        ref=ECOLI_REF,
        index=ECOLI_PICARD_INDEX
    output:
        intervals='local/ecoli/results/{accession}/gatk/realigner.intervals',
        time='local/ecoli/results/{accession}/gatk/bm/realign_target.time',
        trace='local/ecoli/results/{accession}/gatk/bm/realign_target.trace'
    log:
        'local/ecoli/results/{accession}/gatk/log/ecoli_gatk_realign_target_creator.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.gatk} """
            """-T RealignerTargetCreator """
            """-R {input.ref} """
            """-I {input.bam} -o {output.intervals} """
            """>{log} 2>&1"""

# ecoli_gatk_index_marked_bam
#
# Index BAM.
rule ecoli_gatk_index_marked_bam:
    input:
        bam='local/ecoli/temp/{accession}/gatk/align_marked.bam'
    output:
        bai=temp('local/ecoli/temp/{accession}/gatk/align_marked.bai'),
        time='local/ecoli/results/{accession}/gatk/bm/index_mark_dup.time',
        trace='local/ecoli/results/{accession}/gatk/bm/index_mark_dup.trace'
    log:
        'local/ecoli/results/{accession}/gatk/log/ecoli_gatk_index_marked_bam.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.picard} """
            """BuildBamIndex """
            """I={input.bam} """
            """O={output.bai} """
            """>{log} 2>&1"""


# ecoli_gatk_mark_duplicates
#
# Mark duplicate reads.
rule ecoli_gatk_mark_duplicates:
    input:
        bam='local/ecoli/temp/{accession}/gatk/align_init.bam'
    output:
        bam=temp('local/ecoli/temp/{accession}/gatk/align_marked.bam'),
        metrics='local/ecoli/results/{accession}/gatk/mark_dup.metrics',
        time='local/ecoli/results/{accession}/gatk/bm/mark_dup.time',
        trace='local/ecoli/results/{accession}/gatk/bm/mark_dup.trace'
    log:
        'local/ecoli/results/{accession}/gatk/log/ecoli_gatk_mark_duplicates.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.picard} """
            """MarkDuplicatesWithMateCigar """
            """MINIMUM_DISTANCE=500 """
            """INPUT={input.bam} OUTPUT={output.bam} """
            """METRICS_FILE={output.metrics} """
            """>{log} 2>&1"""

# ecoli_gatk_sort_align
#
# Create a sorted BAM from the aligned reads.
rule ecoli_gatk_sort_align:
    input:
        sam='local/ecoli/temp/{accession}/gatk/align_init_nosort.sam'
    output:
        bam=temp('local/ecoli/temp/{accession}/gatk/align_init.bam'),
        time='local/ecoli/results/{accession}/gatk/bm/sort_sam.time',
        trace='local/ecoli/results/{accession}/gatk/bm/sort_sam.trace'
    log:
        'local/ecoli/results/{accession}/gatk/log/ecoli_gatk_sort_align.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.picard} SortSam """
            """I={input.sam} O={output.bam} """
            """SO=coordinate CREATE_INDEX=true """
            """>{log} 2>&1"""

# ecoli_gatk_align
#
# Align sequences against the reference.
rule ecoli_gatk_align:
    input:
        fq_1='local/ecoli/results/{accession}/kestrel/reads/simulated_reads1.fq.gz',
        fq_2='local/ecoli/results/{accession}/kestrel/reads/simulated_reads2.fq.gz',
        ref=ECOLI_REF,
        index=ECOLI_BWA_INDEX_LIST
    output:
        sam=temp('local/ecoli/temp/{accession}/gatk/align_init_nosort.sam'),
        time='local/ecoli/results/{accession}/gatk/bm/bwa.time',
        trace='local/ecoli/results/{accession}/gatk/bm/bwa.trace'
    log:
        'local/ecoli/results/{accession}/gatk/log/ecoli_gatk_align.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/bwa mem -R "@RG\\tID:Ecoli{wildcards.accession}\\tSM:{wildcards.accession}" {input.ref} """
            """{input.fq_1} {input.fq_2} """
            """>{output.sam} 2>{log}"""
