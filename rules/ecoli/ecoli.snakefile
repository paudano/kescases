"""
Master snakefile for the E. coli pipeline.
"""


###################
### Definitions ###
###################



###############
### Include ###
###############

include: 'data.snakefile'
include: 'varcall/assembly.snakefile'
include: 'varcall/kestrel.snakefile'
