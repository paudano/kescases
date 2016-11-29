"""
Rules and variables common to the strep experiment.
"""

from kescaseslib import interval

###################
### Definitions ###
###################

### Accessions ###
STREP_ACCESSIONS = list()

with open(config['strep']['accessions'], 'r') as in_file:
    for line in in_file:
        line = line.strip()

        if not line:
            continue

        STREP_ACCESSIONS.append(line)

STREP_INTERVAL_CONTAINER = interval.IntervalContainer()
STREP_INTERVAL_CONTAINER.add_bed(config['strep']['pbp_bed'])

#############
### Rules ###
#############

include: 'data.snakefile'
include: 'varcall/assembly.snakefile'
include: 'varcall/gatk.snakefile'
include: 'varcall/kestrel.snakefile'
include: 'varcall/quality.snakefile'
