"""
Rules for calling variants by aligning sequence reads and running GATK HaplotypeCaller.
"""

#############
### Rules ###
#############


### Post-GATK Analysis ###

# strep_gatk_get_alignment_pileup
#
# Get alignment pileup
rule strep_gatk_get_alignment_pileup:
    input:
        bam='local/strep/results/{accession}/gatk/sample.bam',
        bai='local/strep/results/{accession}/gatk/sample.bam.bai',
        interval=config['strep']['pbp_bed']
    output:
        pileup='local/strep/results/{accession}/gatk/pileup.tab'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_samtools_index_bam.log'
    shell:
        """bin/samtools mpileup """
            """-f {STREP_REF} """
            """-l {input.interval} """
            """--max-depth 750 """
            """{input.bam} """
            """-o {output.pileup}"""
            """>{log} 2>&1"""

# strep_gatk_samtools_index_bam
#
# Generate a samtools index of the BAM file.
rule strep_gatk_samtools_index_bam:
    input:
        bam='local/strep/results/{accession}/gatk/sample.bam'
    output:
        bai='local/strep/results/{accession}/gatk/sample.bam.bai'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_samtools_index_bam.log'
    shell:
        """bin/samtools index {input} >{log} 2>&1"""


### GATK Pipeline ###

# strep_gatk_call_variants
#
# Call variants.
rule strep_gatk_call_variants:
    input:
        bam='local/strep/results/{accession}/gatk/sample.bam',
        interval=config['strep']['pbp_bed']
    output:
        vcf=temp('local/strep/temp/{accession}/gatk/variants.vcf'),
        time='local/strep/results/{accession}/gatk/bm/variants.time',
        trace='local/strep/results/{accession}/gatk/bm/variants.trace'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_call_variants.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.gatk} """
            """-T HaplotypeCaller """
            """-L {input.interval} """
            """-R {STREP_REF} """
            """-I {input.bam} """
            """-o {output.vcf} """
            """-ploidy 1 """
            """>{log} 2>&1"""

# strep_gatk_indel_realign
#
# Perform indel-based realignment.
rule strep_gatk_indel_realign:
    input:
        bam='local/strep/temp/{accession}/gatk/align_marked.bam',
        bai='local/strep/temp/{accession}/gatk/align_marked.bai',
        intervals='local/strep/results/{accession}/gatk/realigner.intervals',
        ref=STREP_REF
    output:
        bam='local/strep/results/{accession}/gatk/sample.bam',
        time='local/strep/results/{accession}/gatk/bm/realign.time',
        trace='local/strep/results/{accession}/gatk/bm/realign.trace'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_indel_realign.log'
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

# strep_gatk_realign_target_creator
#
# Indel realignment target creator.
rule strep_gatk_realign_target_creator:
    input:
        bam='local/strep/temp/{accession}/gatk/align_marked.bam',
        bai='local/strep/temp/{accession}/gatk/align_marked.bai',
        ref=STREP_REF,
        index=STREP_PICARD_INDEX
    output:
        intervals='local/strep/results/{accession}/gatk/realigner.intervals',
        time='local/strep/results/{accession}/gatk/bm/realign_target.time',
        trace='local/strep/results/{accession}/gatk/bm/realign_target.trace'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_realign_target_creator.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.gatk} """
            """-T RealignerTargetCreator """
            """-R {input.ref} """
            """-I {input.bam} -o {output.intervals} """
            """>{log} 2>&1"""

rule strep_gatk_index_marked_bam:
    input:
        bam='local/strep/temp/{accession}/gatk/align_marked.bam'
    output:
        bai=temp('local/strep/temp/{accession}/gatk/align_marked.bai'),
        time='local/strep/results/{accession}/gatk/bm/index_mark_dup.time',
        trace='local/strep/results/{accession}/gatk/bm/index_mark_dup.trace'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_index_marked_bam.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.picard} """
            """BuildBamIndex """
            """I={input.bam} """
            """O={output.bai} """
            """>{log} 2>&1"""


# strep_gatk_mark_duplicates
#
# Mark duplicate reads.
rule strep_gatk_mark_duplicates:
    input:
        bam='local/strep/temp/{accession}/gatk/align_init.bam'
    output:
        bam=temp('local/strep/temp/{accession}/gatk/align_marked.bam'),
        metrics='local/strep/results/{accession}/gatk/mark_dup.metrics',
        time='local/strep/results/{accession}/gatk/bm/mark_dup.time',
        trace='local/strep/results/{accession}/gatk/bm/mark_dup.trace'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_mark_duplicates.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.picard} """
            """MarkDuplicatesWithMateCigar """
            """MINIMUM_DISTANCE=500 """
            """INPUT={input.bam} OUTPUT={output.bam} """
            """METRICS_FILE={output.metrics} """
            """>{log} 2>&1"""

# strep_gatk_sort_align
#
# Create a sorted BAM from the aligned reads.
rule strep_gatk_sort_align:
    input:
        sam='local/strep/temp/{accession}/gatk/align_init_nosort.sam'
    output:
        bam=temp('local/strep/temp/{accession}/gatk/align_init.bam'),
        time='local/strep/results/{accession}/gatk/bm/sort_sam.time',
        trace='local/strep/results/{accession}/gatk/bm/sort_sam.trace'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_sort_align.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.picard} SortSam """
            """I={input.sam} O={output.bam} """
            """SO=coordinate CREATE_INDEX=true """
            """>{log} 2>&1"""

# strep_gatk_align
#
# Align sequences against the reference.
rule strep_gatk_align:
    input:
        fq_1='local/strep/samples/{accession}/{accession}_1.fastq.gz',
        fq_2='local/strep/samples/{accession}/{accession}_2.fastq.gz',
        ref=STREP_REF,
        index=STREP_BWA_INDEX_LIST
    output:
        sam=temp('local/strep/temp/{accession}/gatk/align_init_nosort.sam'),
        time='local/strep/results/{accession}/gatk/bm/bwa.time',
        trace='local/strep/results/{accession}/gatk/bm/bwa.trace'
    log:
        'local/strep/results/{accession}/gatk/log/strep_gatk_align.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/bwa mem -R "@RG\\tID:Strep{wildcards.accession}\\tSM:{wildcards.accession}" {input.ref} """
            """{input.fq_1} {input.fq_2} """
            """>{output.sam} 2>{log}"""
