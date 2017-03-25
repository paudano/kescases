"""
Rules to summarize data.
"""

rule ecoli_summary_get_call_stats:
    input:
        consexpand()
    output:
        tab='local/ecoli/summary/table/call_stats.tab'

rule ecoli_summary_get_consensus_size:
    input:
        bed=expand('local/ecoli/results/{accession}/kestrel/consensus_regions.bed', accession=ECOLI_ACCESSIONS)
    output:
        tab='local/ecoli/temp/summary/consensus_len.tab'
    shell:
        """
        """