"""
Run Kestrel test cases.
"""

import os

from kescaseslib import bedreader


##################
### Initialize ###
##################

### Read configuration ###
configfile: "config/kescases.json"

### Environment ###

# Save environment as dictionary
ENV = os.environ.copy()

# Home
if 'HOME' in ENV:
    HOME_DIR = ENV['HOME']
else:
    HOME_DIR = None


### Persistent Data: Strep ###

# Variant call regions
STREP_REGIONS = bedreader.bed_interval_to_dataframe(config['strep']['pbp_bed'])


#############
### Rules ###
#############

### Strep ###
include: 'rules/strep/data.snakefile'
include: 'rules/strep/varcall/assembly.snakefile'
include: 'rules/strep/varcall/gatk.snakefile'
include: 'rules/strep/varcall/kestrel.snakefile'
