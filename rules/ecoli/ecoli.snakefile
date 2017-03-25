"""
Master snakefile for the E. coli pipeline.
"""


###################
### Definitions ###
###################

# Read accessions
ECOLI_ACCESSIONS = list()

with open(config['ecoli']['accessions'], 'r') as in_file:
    for line in in_file:
        line = line.strip()

        if not line:
            continue

        ECOLI_ACCESSIONS.append(line)


###############
### Include ###
###############

include: 'data.snakefile'
include: 'varcall/assembly.snakefile'
include: 'varcall/kestrel.snakefile'
include: 'varcall/quality.snakefile'