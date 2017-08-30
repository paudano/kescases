"""
Rules for calling variants with Kestrel.
"""

#############
### Rules ###
#############


# na12878_kestrel_call
#
# Call variants from the IKC file.
rule na12878_kestrel_call:
    input:
        ikc='local/na12878/results/{accession}/kestrel/kmertable.ikc',
        interval=config['na12878']['hetero_bed'],
        ref=NA12878_REF,
    output:
        vcf=temp('local/na12878/temp/{accession,na12878}/kestrel/variants.vcf')
    log:
        time='local/na12878/results/{accession}/kestrel/bm/variants.time',
        trace='local/na12878/results/{accession}/kestrel/bm/variants.trace',
        log='local/na12878/results/{accession}/kestrel/log/kestrel.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {log.trace} """
        """java -Xmx2G -jar {log.kestrel} """
            """-r {input.ref} """
            """-m vcf """
            """--noanchorboth """
            """--varfilter="coverage:0.5,5" """
            """--scanlimitfactor=max """
            """-i {input.interval} """
            """--loglevel=all --logfile={log.log} """
            """-o {output.vcf} """
            """-s{wildcards.accession} """
            """{input.ikc}"""

# na12878_kestrel_make_ikc
#
# Read FASTQ files and k-merize them to an IKC (Indexed K-mer Count) file. Kestrel
# will read k-mer data from this file.
rule na12878_kestrel_make_ikc:
    input:
        fq_1='local/na12878/samples/{accession}/{accession}_1.fastq.gz',
        fq_2='local/na12878/samples/{accession}/{accession}_2.fastq.gz'
    output:
        ikc='local/na12878/results/{accession,na12878}/kestrel/kmertable.ikc',
        seg_size='local/na12878/results/{accession,na12878}/kestrel/bm/seg_size',
        seg_count='local/na12878/results/{accession,na12878}/kestrel/bm/seg_count'
    log:
        time='local/na12878/results/{accession,na12878}/kestrel/bm/kmertable.time',
        trace='local/na12878/results/{accession,na12878}/kestrel/bm/kmertable.trace',
        log='local/na12878/results/{accession}/kestrel/log/kmertable.log'
    shell:
        """mkdir -p local/na12878/temp/{wildcards.accession}/kanalyze; """
        """bin/time -p -o {log.time} """
        """bin/traceproc -o {log.trace} """
        """java -Xmx8G -jar {tools.kanalyze} """
            """count """
            """-k 31 """
            """--countfilter=kmercount:5 """
            """--quality=10 """
            """-m ikc """
            """--minsize 15 """
            """--temploc local/na12878/temp/{wildcards.accession}/kanalyze """
            """--noautodelete """
            """-o {output.ikc} """
            """{input.fq_1} {input.fq_2} """
            """> {log.log}; """
        """ls -l --block-size=1024 local/na12878/temp/{wildcards.accession}/kanalyze | """
            """head -n 1 | """
            """awk '{{print $2}}' > {output.seg_size}; """
        """/bin/ls --block-size=1024 local/na12878/temp/{wildcards.accession}/kanalyze | """
            """wc -l > {output.seg_count}; """
        """rm -rf local/na12878/temp/{wildcards.accession}/kanalyze"""
