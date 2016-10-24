"""
Rules for calling variants with Kestrel.
"""

#############
### Rules ###
#############

rule strep_kestrel_call:
    input:
        ikc='local/strep/results/{accession}/kestrel/kmertable.ikc',
        kestrel=config['tools']['kestrel'],
        interval=config['strep']['pbp_bed']
    output:
        vcf='local/strep/results/{accession}/kestrel/variants.vcf'
    log:
        'local/strep/results/sample/{accession}/kestrel/kestrel.log'
    shell:
        """java -jar {kestrel} """
        """-r {STREP_REF} """
        """-m vcf """
        """--noanchorboth """
        """--varfilter="coverage:0.5,5" """
        """--scanlimitfactor=max """
        """-i {interval} """
        """--loglevel=all --logfile={log} """
        """-o {out.vcf} """
        """{input.ikc}"""

rule strep_kestrel_make_ikc:
    input:
        fq_1='local/strep/samples/{accession}/{accession}_1.fastq.gz',
        fq_2='local/strep/samples/{accession}/{accession}_1.fastq.gz',
        kanalyze=config['tools']['kanalyze'],
        interval_bed=config['strep']['pbp_bed']
    output:
        ikc='local/strep/results/{accession}/kestrel/kmertable.ikc'
    log:
        'local/strep/results/sample/{accession}/kestrel/kanalyze.log'
    shell:
        """java -jar {input.kanalyze} """
        """count """
        """-k 31 """
        """--countfilter=kmercount:5 """
        """--quality=10 """
        """-m ikc """
        """--minsize 15 """
        """-o {output.ikc} """
        """{input.fq_1} {input.fq_2} """
        """> {log}"""

