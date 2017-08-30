"""
Rules for running Kestrel on E. coli samples.
"""


###################
### Definitions ###
###################

def _ecoli_get_art_seed(wildcards):
    """
    Get a random seed to a sample.
    """
    return {accession: seed for seed, accession in zip(SEED_LIST, ECOLI_ACCESSIONS)}[wildcards.accession]

#############
### Rules ###
#############

# ecoli_kestrel_call_all
#
# Get calls for all samples.
rule ecoli_kestrel_call_all:
    input:
        tab=expand('local/ecoli/results/{accession}/kestrel/variants.tab', accession=ECOLI_ACCESSIONS)
    run:
        pass

# ecoli_kestrel_call
#
# Call variants from the IKC file.
rule ecoli_kestrel_call:
    input:
        ikc='local/ecoli/results/{accession}/kestrel/kmertable.ikc',
        ref=ECOLI_REF,
        consensus_bed='local/ecoli/results/{accession}/kestrel/consensus_regions.bed'
    output:
        vcf=temp('local/ecoli/temp/{accession}/kestrel/variants.vcf'),
        sam=temp('local/ecoli/temp/{accession}/kestrel/haplotypes.sam'),
        bam='local/ecoli/results/{accession}/kestrel/haplotypes.bam',
        bai='local/ecoli/results/{accession}/kestrel/haplotypes.bam.bai',
        time='local/ecoli/results/{accession}/kestrel/bm/variants.time',
        trace='local/ecoli/results/{accession}/kestrel/bm/variants.trace'
    log:
        'local/ecoli/results/{accession}/kestrel/log/kestrel.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -Xmx24G -jar {tools.kestrel} """
            """-r {input.ref} """
            """-m vcf """
            """--noanchorboth """
            """--varfilter="coverage:0.5,5" """
            """--scanlimitfactor=7 """
            """--noambivar """
            """-i {input.consensus_bed} """
            """--loglevel=all --logfile={log} """
            """-o {output.vcf} """
            """-p {output.sam} """
            """-s{wildcards.accession} """
            """{input.ikc}; """
        """samtools view -b {output.sam} > {output.bam}; """
        """samtools index {output.bam}"""

# ecoli_kestrel_make_ikc
#
# Read FASTQ files and k-merize them to an IKC (Indexed K-mer Count) file. Kestrel
# will read k-mer data from this file.
rule ecoli_kestrel_make_ikc:
    input:
        fq_1='local/ecoli/results/{accession}/kestrel/reads/simulated_reads1.fq.gz',
        fq_2='local/ecoli/results/{accession}/kestrel/reads/simulated_reads2.fq.gz'
    output:
        ikc='local/ecoli/results/{accession}/kestrel/kmertable.ikc',
        time='local/ecoli/results/{accession}/kestrel/bm/kmertable.time',
        trace='local/ecoli/results/{accession}/kestrel/bm/kmertable.trace',
        seg_size='local/ecoli/results/{accession}/kestrel/bm/seg_size',
        seg_count='local/ecoli/results/{accession}/kestrel/bm/seg_count'
    log:
        'local/ecoli/results/{accession}/kestrel/log/kmertable.log'
    shell:
        """mkdir -p local/ecoli/temp/{wildcards.accession}/kanalyze; """
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -Xmx3G -jar {tools.kanalyze} """
            """count """
            """-k 31 """
            """--countfilter=kmercount:5 """
            """--quality=10 """
            """-m ikc """
            """--minsize 15 """
            """--temploc local/ecoli/temp/{wildcards.accession}/kanalyze """
            """--noautodelete """
            """-o {output.ikc} """
            """{input.fq_1} {input.fq_2} """
            """> {log}; """
        """ls -l --block-size=1024 local/ecoli/temp/{wildcards.accession}/kanalyze | """
            """head -n 1 | """
            """awk '{{print $2}}' > {output.seg_size}; """
        """/bin/ls --block-size=1024 local/ecoli/temp/{wildcards.accession}/kanalyze | """
            """wc -l > {output.seg_count}; """
        """rm -rf local/ecoli/temp/{wildcards.accession}/kanalyze"""

# ecoli_make_consensus_bed
#
# Get a BED file of consensus regions.
rule ecoli_make_consensus_bed:
    input:
        bam='local/ecoli/results/{accession}/assemble/contig.bam',
        bed='local/ecoli/results/{accession}/assemble/no_consensus.bed'
    output:
        tmp_con=temp('local/ecoli/temp/{accession}/kestrel/no_consensus_bed6.bed'),
        tmp_aln=temp('local/ecoli/temp/{accession}/kestrel/contigs_bed6.bed'),
        bed='local/ecoli/results/{accession}/kestrel/consensus_regions.bed'
    shell:
        """awk -vOFS="\\t" '{{print $1, $2, $3, $4, "0", "+"}}' {input.bed} > {output.tmp_con}; """
        """bamToBed -i {input.bam} > {output.tmp_aln}; """
        """bedtools subtract -a {output.tmp_aln} -b {output.tmp_con} > {output.bed}"""

# ecoli_kestrel_make_pileup
#
# Make pileup of the simulated reads.
rule ecoli_kestrel_make_pileup:
    input:
        sam='local/ecoli/results/{accession}/kestrel/reads/simulated_reads.bam',
        ref=ECOLI_REF
    output:
        pileup='local/ecoli/results/{accession}/kestrel/reads/pileup.tab'
    shell:
        """bin/samtools mpileup """
            """-f {input.ref} """
            """{input.sam} """
            """-o {output.pileup} """
            """>{log} 2>&1"""

# ecoli_compress_simulated_reads
#
# Compress output from simulated reads
rule ecoli_compress_simulated_reads:
    input:
        fq_1='local/ecoli/temp/{accession}/kestrel/reads/simulated_reads1.fq',
        fq_2='local/ecoli/temp/{accession}/kestrel/reads/simulated_reads2.fq',
        aln_1='local/ecoli/temp/{accession}/kestrel/reads/simulated_reads1.aln',
        aln_2='local/ecoli/temp/{accession}/kestrel/reads/simulated_reads2.aln',
        sam='local/ecoli/temp/{accession}/kestrel/reads/simulated_reads.sam'
    output:
        fq_1='local/ecoli/results/{accession}/kestrel/reads/simulated_reads1.fq.gz',
        fq_2='local/ecoli/results/{accession}/kestrel/reads/simulated_reads2.fq.gz',
        aln_1='local/ecoli/results/{accession}/kestrel/reads/simulated_reads1.aln.gz',
        aln_2='local/ecoli/results/{accession}/kestrel/reads/simulated_reads2.aln.gz',
        bam='local/ecoli/results/{accession}/kestrel/reads/simulated_reads.bam',
        bai='local/ecoli/results/{accession}/kestrel/reads/simulated_reads.bam.bai'
    shell:
        """gzip -c {input.fq_1} > {output.fq_1}; """
        """gzip -c {input.fq_2} > {output.fq_2}; """
        """gzip -c {input.aln_1} > {output.aln_1}; """
        """gzip -c {input.aln_2} > {output.aln_2}; """
        """samtools view -b {input.sam} | samtools sort > {output.bam}; """
        """samtools index {output.bam}"""

# ecoli_kestrel_simulate_reads
#
# Simulate reads from assemblies for variant-calling by Kestrel.
rule ecoli_kestrel_simulate_reads:
    input:
        contigs='local/ecoli/temp/{accession}/kestrel/reads/contigs.fa'
    output:
        fq_1=temp('local/ecoli/temp/{accession}/kestrel/reads/simulated_reads1.fq'),
        fq_2=temp('local/ecoli/temp/{accession}/kestrel/reads/simulated_reads2.fq'),
        aln_1=temp('local/ecoli/temp/{accession}/kestrel/reads/simulated_reads1.aln'),
        aln_2=temp('local/ecoli/temp/{accession}/kestrel/reads/simulated_reads2.aln'),
        sam=temp('local/ecoli/temp/{accession}/kestrel/reads/simulated_reads.sam')
    params:
        seed=_ecoli_get_art_seed
    log:
        'local/ecoli/results/{accession}/kestrel/log/simulated_reads.log'
    shell:
        """bin/art_illumina -sam -i {input.contigs} -rs {params.seed} -p -l 150 -ss HS25 -f 30 -m 200 -s 10 -o local/ecoli/temp/{wildcards.accession}/kestrel/reads/simulated_reads >{log} 2>&1"""

# ecoli_kestrel_uncompress_contigs
#
# Uncompress contigs for art_illumina
rule ecoli_kestrel_uncompress_contigs:
    input:
        contigs='local/ecoli/samples/{accession}.fa.gz'
    output:
        contigs=temp('local/ecoli/temp/{accession}/kestrel/reads/contigs.fa')
    shell:
        """gunzip -c {input.contigs} > {output.contigs}"""
