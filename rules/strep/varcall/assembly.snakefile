"""
Rules for assembling contigs and calling variants from the assemblies.
"""

from kescaseslib import variant
from kescaseslib import pileup


#############
### Rules ###
#############

### Plots ###

# strep_asm_plot_tp_assembly
#
# Generate plots showing memory usage trends over time.
rule strep_asm_plot_tp_assembly:
    input:
        tp='local/strep/results/{accession}/assemble/bm/spades.trace'
    output:
        eps='local/strep/figures/sample/{accession}/assemble/assemble_rss.eps',
        pdf='local/strep/figures/sample/{accession}/assemble/assemble_rss.pdf'
    shell:
        """bin/Rscript scripts/strep/plots/assembly_trace_plot.R {input.tp} {output.eps}; """
        """bin/Rscript scripts/strep/plots/assembly_trace_plot.R {input.tp} {output.pdf}"""


### Run Pipeline ###

# strep_asm_pileup_to_vcf
#
# Call variants from the pileup, convert to VCF, and generate a BED file of loci where
# the assembly does not have a consensus (depth != 1).
rule strep_asm_pileup_to_vcf:
    input:
        pileup='local/strep/results/{accession}/assemble/pileup.tab',
        ref=STREP_REF,
        fai=STREP_REF_FAI
    output:
        vcf=temp('local/strep/temp/{accession}/assemble/variants.vcf'),
        bed='local/strep/results/{accession}/assemble/no_consensus.bed'
    run:
        variant_list = pileup.read_pileup_variants(input.pileup, wildcards.accession, STREP_INTERVAL_CONTAINER, output.bed)
        variant.write_vcf(variant_list, wildcards.accession, output.vcf, input.ref)

# strep_asm_pileup
#
# Generate a pileup from the aligned contig BAM file.
rule strep_asm_pileup:
    input:
        bam='local/strep/results/{accession}/assemble/contig.bam',
        bai='local/strep/results/{accession}/assemble/contig.bam.bai',
        interval=config['strep']['pbp_bed'],
        ref=STREP_REF
    output:
        pileup='local/strep/results/{accession}/assemble/pileup.tab',
        time='local/strep/results/{accession}/assemble/bm/pileup.time',
        trace='local/strep/results/{accession}/assemble/bm/pileup.trace'
    log:
        'local/strep/results/{accession}/assemble/log/mpileup.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/samtools mpileup """
            """-f {input.ref} """
            """-l {input.interval} """
            """{input.bam} """
            """-o {output.pileup} """
            """>{log} 2>&1"""

# strep_asm_bam_index
#
# Index aligned assembly BAM.
rule strep_asm_bam_index:
    input:
        bam='local/strep/results/{accession}/assemble/contig.bam'
    output:
        bai='local/strep/results/{accession}/assemble/contig.bam.bai',
        time='local/strep/results/{accession}/assemble/bm/indexbam.time',
        trace='local/strep/results/{accession}/assemble/bm/indexbam.trace'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/samtools index {input.bam}"""

# strep_asm_sort_sam
#
# Sort aligned contigs and index.
rule strep_asm_sort_sam:
    input:
        sam='local/strep/temp/{accession}/assemble/contig.sam'
    output:
        bam='local/strep/results/{accession}/assemble/contig.bam',
        time='local/strep/results/{accession}/assemble/bm/sortsam.time',
        trace='local/strep/results/{accession}/assemble/bm/sortsam.trace'
    log:
        'local/strep/results/{accession}/assemble/log/sort_contig_alignment.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -jar {tools.picard} """
            """SortSam """
            """I={input.sam} """
            """O={output.bam} """
            """SO=coordinate CREATE_INDEX=true """
            """2>{log}"""

# strep_asm_align
#
# Align assembled contigs to a reference.
rule strep_asm_align:
    input:
        scaffolds='local/strep/results/{accession}/assemble/spades/scaffolds.fasta',
        ref=STREP_REF,
        index=STREP_BWA_INDEX_LIST
    output:
        sam=temp('local/strep/temp/{accession}/assemble/contig.sam'),
        time='local/strep/results/{accession}/assemble/bm/bwa.time',
        trace='local/strep/results/{accession}/assemble/bm/bwa.trace'
    log:
        'local/strep/results/{accession}/assemble/log/strep_asm_align.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/bwa mem """
            """-B 1 """
            """-L 1000,1000 """
            """-R "@RG\\tID:Strep{wildcards.accession}\\tSM:{wildcards.accession}" """
            """{STREP_REF} """
            """{input.scaffolds} """
            """> {output.sam} """
            """2>{log}"""

# strep_asm_assemble
#
# Assemble contigs for a sample.
rule strep_asm_assemble:
    input:
        fq_1='local/strep/samples/{accession}/{accession}_1.fastq.gz',
        fq_2='local/strep/samples/{accession}/{accession}_2.fastq.gz'
    output:
        scaffolds='local/strep/results/{accession}/assemble/spades/scaffolds.fasta',
        time='local/strep/results/{accession}/assemble/bm/spades.time',
        trace='local/strep/results/{accession}/assemble/bm/spades.trace'
    log:
        'local/strep/results/{accession}/assemble/log/spades.log'
    shell:
        """rm -rf $(dirname {output.scaffolds}); """
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """bin/spades.py """
            """-1 {input.fq_1} -2 {input.fq_2} """
            """-o $(dirname {output.scaffolds}) """
            """> {log}"""
