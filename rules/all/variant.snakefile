"""
Rules for variant calling that multiple pipelines may use.
"""

#############
### Rules ###
#############

# vcf_to_tabix
#
# bgzip and index a variants VCF file.
rule vcf_to_tabix:
    input:
        vcf='local/{experiment}/temp/{accession}/{pipeline}/variants.vcf'
    output:
        vcf_gz='local/{experiment}/results/{accession}/{pipeline}/variants.vcf.gz',
        vcf_tabix='local/{experiment}/results/{accession}/{pipeline}/variants.vcf.gz.tbi'
    log:
        'local/{experiment}/results/{accession}/{pipeline}/log/tabix.log'
    shell:
        """bin/bgzip -c {input.vcf} > {output.vcf_gz} 2>{log}; """
        """bin/tabix -f {output.vcf_gz} 2>>{log}"""
