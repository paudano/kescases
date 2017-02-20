"""
Rules for running Kestrel on E. coli samples.
"""


#############
### Rules ###
#############

# ecoli_kestrel_call
#
# Call variants from the IKC file.
rule ecoli_kestrel_call:
    input:
        ikc='local/ecoli/results/{accession}/kestrel/kmertable.ikc',
        ref=ECOLI_REF
    output:
        vcf=temp('local/ecoli/temp/{accession}/kestrel/variants.vcf'),
        time='local/ecoli/results/{accession}/kestrel/bm/variants.time',
        trace='local/ecoli/results/{accession}/kestrel/bm/variants.trace'
    log:
        'local/ecoli/results/{accession}/kestrel/log/kestrel.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -Xmx2G -jar {tools.kestrel} """
            """-r {input.ref} """
            """-m vcf """
            """--noanchorboth """
            """--varfilter="coverage:0.5,5" """
            """--scanlimitfactor=max """
            """--loglevel=all --logfile={log} """
            """-o {output.vcf} """
            """-s{wildcards.accession} """
            """{input.ikc}"""

# ecoli_kestrel_make_ikc
#
# Read FASTQ files and k-merize them to an IKC (Indexed K-mer Count) file. Kestrel
# will read k-mer data from this file.
rule ecoli_kestrel_make_ikc:
    input:
        fq_1='local/ecoli/results/{accession}/kestrel/reads/simulated_reads1.fq',
        fq_2='local/ecoli/results/{accession}/kestrel/reads/simulated_reads2.fq'
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

# ecoli_kestrel_simulate_reads
#
# Simulate reads from assemblies for variant-calling by Kestrel.
rule ecoli_kestrel_simulate_reads:
    input:
        scaffolds='local/ecoli/temp/{accession}/kestrel/reads/scaffolds.fa'
    output:
        fq_1='local/ecoli/results/{accession}/kestrel/reads/simulated_reads1.fq',
        fq_2='local/ecoli/results/{accession}/kestrel/reads/simulated_reads2.fq'
    log:
        'local/ecoli/results/{accession}/kestrel/log/simulated_reads.log'
    shell:
        """bin/art_illumina -sam -i {input.scaffolds} -p -l 150 -ss HS25 -f 30 -m 200 -s 10 -o local/ecoli/results/{wildcards.accession}/kestrel/reads/simulated_reads >{log} 2>&1"""

# ecoli_kestrel_uncompress_scaffolds
#
# Uncompress scaffolds for art_illumina
rule ecoli_kestrel_uncompress_scaffolds:
    input:
        scaffolds='local/ecoli/samples/{accession}.scaffolds_min500bp.fa.gz'
    output:
        scaffolds=temp('local/ecoli/temp/{accession}/kestrel/reads/scaffolds.fa')
    shell:
        """gunzip -c {input.scaffolds} > {output.scaffolds}"""