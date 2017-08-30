"""
Run tests on NA12878.
"""


##################
### Defnitions ###
##################

### Accessions ###
NA12878_ACCESSIONS = list()

with open(config['na12878']['accessions'], 'r') as in_file:
    for line in in_file:
        line = line.strip()

        if not line:
            continue

        NA12878_ACCESSIONS.append(line)

### Sequence File List ###

def _na12878_get_input_fq(wildcards):
    """
    Get a list of all input sequence files.

    :param wildcards: Ignored.
    """
    seq_file_list = list()

    for acc in NA12878_ACCESSIONS:
        seq_file_list.append('local/na12878/samples/{}/{}_1.fastq.gz'.format(acc, acc))
        seq_file_list.append('local/na12878/samples/{}/{}_2.fastq.gz'.format(acc, acc))

    return seq_file_list

###############
### Include ###
###############

#include: 'varcall/kestrel.snakefile'


#############
### Rules ###
#############

# na12878_fetch
#
# Fetch NA12878 data.
rule na12878_fetch:
    input:
        fq=_na12878_get_input_fq
