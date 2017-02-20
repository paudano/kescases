"""
Rules for calling variants from the E. coli assemblies.
"""

from kescaseslib import variant
from kescaseslib import pileup
from kescaseslib import interval


#############
### Rules ###
#############

# ecoli_asm_pileup_to_vcf
#
# Call variants from the pileup, convert to VCF, and generate a BED file of loci where
# the assembly does not have a consensus (depth != 1).
rule ecoli_asm_pileup_to_vcf:
    input:
        pileup='local/ecoli/results/{accession}/assemble/pileup.tab',
        ref=ECOLI_REF,
        fai=ECOLI_REF_FAI,
        ref_bed=ECOLI_REF_BED
    output:
        vcf=temp('local/ecoli/temp/{accession}/assemble/variants.vcf'),
        bed='local/ecoli/results/{accession}/assemble/no_consensus.bed'
    run:
        interval_container = interval.IntervalContainer()
        interval_container.add_bed(input.ref_bed)

        variant_list = pileup.read_pileup_variants(input.pileup, wildcards.accession, interval_container, output.bed)
        variant.write_vcf(variant_list, wildcards.accession, output.vcf, input.ref)

# ecoli_asm_pileup
#
# Generate a pileup from the aligned contig BAM file.
rule ecoli_asm_pileup:
    input:
        bam='local/ecoli/results/{accession}/assemble/contig.bam',
        bai='local/ecoli/results/{accession}/assemble/contig.bam.bai',
        ref=ECOLI_REF
    output:
        pileup='local/ecoli/results/{accession}/assemble/pileup.tab',
        time='local/ecoli/results/{accession}/assemble/bm/pileup.time',
        trace='local/ecoli/results/{accession}/assemble/bm/pileup.trace'
    log:
        'local/ecoli/results/{accession}/assemble/log/mpileup.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/samtools mpileup """
            """-f {input.ref} """
            """{input.bam} """
            """-o {output.pileup} """
            """>{log} 2>&1"""

# ecoli_asm_bam_index
#
# Index aligned assembly BAM.
rule ecoli_asm_bam_index:
    input:
        bam='local/ecoli/results/{accession}/assemble/contig.bam'
    output:
        bai='local/ecoli/results/{accession}/assemble/contig.bam.bai',
        time='local/ecoli/results/{accession}/assemble/bm/indexbam.time',
        trace='local/ecoli/results/{accession}/assemble/bm/indexbam.trace'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/samtools index {input.bam}"""

# ecoli_asm_sort_sam
#
# Sort aligned contigs and index.
rule ecoli_asm_sort_sam:
    input:
        sam='local/ecoli/temp/{accession}/assemble/contig.sam'
    output:
        bam='local/ecoli/results/{accession}/assemble/contig.bam',
        time='local/ecoli/results/{accession}/assemble/bm/sortsam.time',
        trace='local/ecoli/results/{accession}/assemble/bm/sortsam.trace'
    log:
        'local/ecoli/results/{accession}/assemble/log/sort_contig_alignment.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.picard} """
            """SortSam """
            """I={input.sam} """
            """O={output.bam} """
            """SO=coordinate CREATE_INDEX=true """
            """2>{log}"""

# ecoli_asm_align
#
# Align assembled contigs to a reference.
rule ecoli_asm_align:
    input:
        scaffolds='local/ecoli/samples/{accession}.scaffolds_min500bp.fa.gz',
        ref=ECOLI_REF,
        index=ECOLI_BWA_INDEX_LIST
    output:
        sam=temp('local/ecoli/temp/{accession}/assemble/contig.sam'),
        time='local/ecoli/results/{accession}/assemble/bm/bwa.time',
        trace='local/ecoli/results/{accession}/assemble/bm/bwa.trace'
    log:
        'local/ecoli/results/{accession}/assemble/log/asm_align.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/bwa mem """
            """-B 1 """
            """-L 1000,1000 """
            """-R "@RG\\tID:Ecoli-{wildcards.accession}\\tSM:{wildcards.accession}" """
            """{input.ref} """
            """{input.scaffolds} """
            """> {output.sam} """
            """2>{log}"""
