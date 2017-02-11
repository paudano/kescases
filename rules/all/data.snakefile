"""
Download and manage data.
"""

#############
### Rules ###
#############

# strep_get_sample
#
# Dump paired-end FASTQ files from SRA file.
rule data_get_sample_fastq:
    input:
        sra='local/{exp}/temp/samples/{accession}.sra'
    output:
        fastq_1='local/{exp,strep|mlst}/samples/{accession}/{accession}_1.fastq.gz',
        fastq_2='local/{exp,strep|mlst}/samples/{accession}/{accession}_2.fastq.gz'
    log:
        'local/{exp}/samples/log/{accession}/strep_get_sample_fastq.log'
    run:
        # Download
        out_dir = os.path.dirname(output.fastq_1)
        shell("""bin/fastq-dump --split-files --gzip {input.sra} --outdir {out_dir} >{log} 2>&1""")

# strep_get_sample_sra
#
# Download SRA file for sample.
rule data_get_sample_sra:
    output:
        sra=temp('local/{exp,strep|mlst}/temp/samples/{accession,(SRR|ERR)\d+}.sra')
    log:
        'local/{exp}/samples/log/{accession}/strep_get_sample_sra.log'
    run:
        sra_url = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{}/{}/{}/{}.sra'.format(
            wildcards.accession[0:3], wildcards.accession[0:6], wildcards.accession, wildcards.accession
        )

        dest_dir = os.path.dirname(output.sra)

        shell('wget -qc {sra_url} -P {dest_dir} -o {log}')
