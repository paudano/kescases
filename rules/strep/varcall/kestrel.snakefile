"""
Rules for calling variants with Kestrel.
"""

#############
### Rules ###
#############


# strep_kestrel_call
#
# Call variants from the IKC file.
rule strep_kestrel_call:
    input:
        ikc='local/strep/results/{accession}/kestrel/kmertable.ikc',
        interval=config['strep']['pbp_bed'],
        ref=STREP_REF,
    output:
        vcf=temp('local/strep/temp/{accession}/kestrel/variants.vcf'),
        time='local/strep/results/{accession}/kestrel/bm/variants.time',
        trace='local/strep/results/{accession}/kestrel/bm/variants.trace'
    log:
        'local/strep/results/{accession}/kestrel/log/kestrel.log'
    shell:
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -Xmx2G -jar {tools.kestrel} """
            """-r {input.ref} """
            """-m vcf """
            """--noanchorboth """
            """--varfilter="coverage:0.5,5" """
            """--scanlimitfactor=max """
            """-i {input.interval} """
            """--loglevel=all --logfile={log} """
            """-o {output.vcf} """
            """-s{wildcards.accession} """
            """{input.ikc}"""

# strep_kestrel_make_ikc
#
# Read FASTQ files and k-merize them to an IKC (Indexed K-mer Count) file. Kestrel
# will read k-mer data from this file.
rule strep_kestrel_make_ikc:
    input:
        fq_1='local/strep/samples/{accession}/{accession}_1.fastq.gz',
        fq_2='local/strep/samples/{accession}/{accession}_2.fastq.gz'
    output:
        ikc='local/strep/results/{accession}/kestrel/kmertable.ikc',
        time='local/strep/results/{accession}/kestrel/bm/kmertable.time',
        trace='local/strep/results/{accession}/kestrel/bm/kmertable.trace',
        seg_size='local/strep/results/{accession}/kestrel/bm/seg_size',
        seg_count='local/strep/results/{accession}/kestrel/bm/seg_count'
    log:
        'local/strep/results/{accession}/kestrel/log/kmertable.log'
    shell:
        """mkdir -p local/strep/temp/{wildcards.accession}/kanalyze; """
        """bin/time -p -o {output.time} """
        """bin/traceproc -o {output.trace} """
        """java -Xmx8G -jar {tools.kanalyze} """
            """count """
            """-k 31 """
            """--countfilter=kmercount:5 """
            """--quality=10 """
            """-m ikc """
            """--minsize 15 """
            """--temploc local/strep/temp/{wildcards.accession}/kanalyze """
            """--noautodelete """
            """-o {output.ikc} """
            """{input.fq_1} {input.fq_2} """
            """> {log}; """
        """ls -l --block-size=1024 local/strep/temp/{wildcards.accession}/kanalyze | """
            """head -n 1 | """
            """awk '{{print $2}}' > {output.seg_size}; """
        """/bin/ls --block-size=1024 local/strep/temp/{wildcards.accession}/kanalyze | """
            """wc -l > {output.seg_count}; """
        """rm -rf local/strep/temp/{wildcards.accession}/kanalyze"""
