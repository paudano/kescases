"""
Run ANI (Average Nucleotide Identity) to determine relative phylogeny among the samples and reference sequence.
No outgroup is used (not a full phylogeny analysis).
"""

import os

# strep_ani_run_ani
#
# Run ANI on all samples and the reference.
rule strep_ani_run_ani:
    input:
        reference='local/strep/reference/NC_003028.fasta',
        scaffolds=expand('local/strep/results/{accession}/assemble/spades/scaffolds.fasta', accession=STREP_ACCESSIONS)
    output:
        tab_pct='local/strep/ani/ani_tables/ANIm_percentage_identity.tab'
    run:

        # Make output directory
        os.makedirs('local/strep/ani/ani_input', exist_ok=True)

        # Link reference
        shell("""ln -sf $(readlink -f {input.reference}) local/strep/ani/ani_input/$(basename {input.reference})""")

        # Link samples
        for accession in STREP_ACCESSIONS:
            shell("""ln -sf $(readlink -f local/strep/results/{accession}/assemble/spades/scaffolds.fasta) local/strep/ani/ani_input/{accession}.fasta""")

        # Run ANI
        shell("""bin/python scripts/pyani/average_nucleotide_identity.py -f -m ANIm -i local/strep/ani/ani_input/ -o local/strep/ani/ani_tables/""")
