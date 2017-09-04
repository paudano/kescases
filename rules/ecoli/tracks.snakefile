"""
Make UCSC tracks.
"""

from Bio import SeqIO
import pandas as pd
import numpy as np


#
# Track Hub
#

# ecoli_tracks_make_hub
#
# Make track hub.
rule ecoli_tracks_make_hub:
    input:
        bb_gene='local/ecoli/tracks/hub/NC_000913.3/gene.bb',
        bb_repeat='local/ecoli/tracks/hub/NC_000913.3/repeat.bb',
        bb_mobile='local/ecoli/tracks/hub/NC_000913.3/mobile.bb',
        bb_rna='local/ecoli/tracks/hub/NC_000913.3/rna.bb',
        bb_misc='local/ecoli/tracks/hub/NC_000913.3/misc.bb',
        bb_kes='local/ecoli/tracks/hub/NC_000913.3/kes_only.bb',
        bw_kes_var='local/ecoli/tracks/hub/NC_000913.3/var_depth_kestrel.bw',
        bw_gatk_var='local/ecoli/tracks/hub/NC_000913.3/var_depth_gatk.bw',
        bit='local/ecoli/tracks/hub/NC_000913.3/NC_000913.3.2bit',
        hub='data/ecoli/tracks/hub.txt',
        genomes='data/ecoli/tracks/genomes.txt',
        trackdb='data/ecoli/tracks/trackDb.txt',
    output:
        hub='local/ecoli/tracks/hub/hub.txt',
        genomes='local/ecoli/tracks/hub/genomes.txt',
        trackdb='local/ecoli/tracks/hub/NC_000913.3/trackDb.txt'
    shell:
        """cp -f {input.hub} {output.hub}; """
        """cp -f {input.genomes} {output.genomes}; """
        """cp -f {input.trackdb} {output.trackdb}"""

# ecoli_tracks_fa_to_2bit
#
# FASTA to 2bit.
rule ecoli_tracks_fa_to_2bit:
    input:
        fa=ECOLI_REF
    output:
        bit='local/ecoli/tracks/hub/NC_000913.3/NC_000913.3.2bit'
    shell:
        """faToTwoBit {input.fa} {output.bit}"""


#
# Variant call depth
#

# ecoli_tracks_var_depth_bw
#
# Variant depth to BigWig.
rule ecoli_tracks_var_depth_bw:
    input:
        bed='local/ecoli/summary/bed/variant_depth_{pipeline}.bed',
        fai=ECOLI_REF_FAI
    output:
        bw='local/ecoli/tracks/hub/NC_000913.3/var_depth_{pipeline,kestrel|gatk}.bw'
    shell:
        """bedGraphToBigWig -blockSize=20 {input.bed} {input.fai} {output.bw}"""

#
# Kestrel-only loci
#

# ecoli_tracks_kes_only_bb
#
# Kestrel-only BigBed.
rule ecoli_tracks_kes_only_bb:
    input:
        bed='local/ecoli/temp/tracks/hub/kes_only.bed',
        asq='data/ecoli/tracks/kes_only.as',
        fai=ECOLI_REF_FAI
    output:
        bb='local/ecoli/tracks/hub/NC_000913.3/kes_only.bb'
    shell:
        """bedToBigBed -type=bed9+ -as={input.asq} {input.bed} {input.fai} {output.bb}"""

# ecoli_tracks_kes_only_bed
#
# Make track BED.
rule ecoli_tracks_kes_only_bed:
    input:
        bed='local/ecoli/summary/bed/kes_only_merged.bed'
    output:
        bed=temp('local/ecoli/temp/tracks/hub/kes_only.bed')
    run:

        df = pd.read_table(input.bed, header=0)

        df['POS1'] = df['POS'] + 1
        df['ID'] = df.apply(lambda row: 'KESONLY-{POS1}-{LENGTH}-{DEPTH}'.format(**row), axis=1)
        del(df['POS1'])

        df['SCORE'] = 0
        df['STRAND'] = '.'

        df['THICK_START'] = df['POS']
        df['THICK_END'] = df['END']

        df['ITEM_RGB'] = '0,0,0'

        df = df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SCORE', 'STRAND', 'THICK_START', 'THICK_END', 'ITEM_RGB', 'LENGTH', 'DEPTH')]

        df.to_csv(output.bed, sep='\t', index=False)

#
# Reference Annotation BED
#

# ecoli_tracks_make_bb
#
# Make BigBed.
rule ecoli_tracks_make_bb:
    input:
        bed='local/ecoli/tracks/bed/{track}.bed',
        asq='local/ecoli/tracks/bed/{track}.as',
        fai=ECOLI_REF_FAI
    output:
        bb='local/ecoli/tracks/hub/NC_000913.3/{track}.bb'
    run:

        df = pd.read_table(input.bed, header=0)

        extra = df.shape[1] - 9

        shell("""bedToBigBed -type=bed9+{extra} -as={input.asq} {input.bed} {input.fai} {output.bb}""")

# ecoli_tracks_make_as
#
# Make AS file for BigBed.
rule ecoli_tracks_make_as:
    input:
        bed='local/ecoli/tracks/bed/{track}.bed'
    output:
        asq='local/ecoli/tracks/bed/{track}.as'
    params:
        chrom='NC_000913.3'
    run:

        df = pd.read_table(input.bed, header=0)
        cols = df.columns[9:]

        track_name_cap = wildcards.track
        track_name_cap = track_name_cap[0].upper() + track_name_cap[1:]

        with open(output.asq, 'w') as out_file:
            out_file.write('table {}\n'.format(wildcards.track))
            out_file.write('"{} annotations on {}"\n'.format(track_name_cap, params.chrom))
            out_file.write('(\n')
            out_file.write('  string  chrom;        "ASSEMBLY"\n')
            out_file.write('  uint    chromStart;   "POS"\n')
            out_file.write('  uint    chromEnd;     "END"\n')
            out_file.write('  string  name;         "NAME"\n')
            out_file.write('  uint    score;        "SCORE"\n')
            out_file.write('  char[1] strand;       "STRAND"\n')
            out_file.write('  uint    thickStart;   "CDS_POS"\n')
            out_file.write('  uint    thickEnd;     "CDS_END"\n')
            out_file.write('  uint    reserved;	    "COLOR"\n')

            for col in cols:
                out_file.write('  string  {};  "{}"\n'.format(col, col))

            out_file.write(')\n')

# ecoli_tracks_make_bed
#
# Make BED files.
rule ecoli_tracks_make_bed:
    input:
        gb=ECOLI_REF_GB
    output:
        gene='local/ecoli/tracks/bed/gene.bed',
        repeat='local/ecoli/tracks/bed/repeat.bed',
        mobile='local/ecoli/tracks/bed/mobile.bed',
        rna='local/ecoli/tracks/bed/rna.bed',
        misc='local/ecoli/tracks/bed/misc.bed'
    params:
        chrom='NC_000913.3'
    run:

        # Definitions
        def make_table(anno_list, rgb='0,0,0'):
            df = pd.concat(anno_list, axis=1).T
            df['#CHROM'] = params.chrom

            if 'ID' not in df.columns:
                df['ID'] = '.'

            if 'STRAND' not in df.columns:
                df['STRAND'] = '.'

            df['SCORE'] = 0
            df['THICK_START'] = df['POS']
            df['THICK_END'] = df['END']
            df['ITEM_RGB'] = rgb

            # Rearrange columns
            cols = ['#CHROM', 'POS', 'END', 'ID', 'SCORE', 'STRAND', 'THICK_START', 'THICK_END', 'ITEM_RGB']

            cols.extend(sorted([col for col in df.columns if col not in cols]))

            df = df.loc[:, cols]

            df.fillna('.', inplace=True)

            return df

        # Get lists of records
        gene_list = list()
        repeat_list = list()
        mobile_list = list()
        rna_list = list()
        misc_list = list()

        with open(input.gb, 'r') as in_file:
            for record in SeqIO.parse(in_file, 'genbank'):
                for feature in record.features:

                    # All
                    item = {attr.upper(): ','.join(val).replace(' ', '') for attr, val in feature.qualifiers.items()}
                    item = {attr: val.strip() for attr, val in item.items()}
                    item = {attr: val if val else '.' for attr, val in item.items()}
                    item = {attr: val if len(val) < 255 else val[:252] + '...' for attr, val in item.items()}


                    item['POS'] = feature.location.start.position - 1
                    item['END'] = feature.location.end.position

                    item['STRAND'] = '-' if feature.strand < 0 else '+'

                    # Get annotations for type
                    if feature.type == 'gene':
                        item['ID'] = item['GENE'] if 'GENE' in item else item['LOCUS_TAG']

                        gene_list.append(pd.Series(item))

                    elif feature.type == 'repeat_region':
                        repeat_list.append(pd.Series(item))

                    elif feature.type == 'mobile_element':
                        mobile_list.append(pd.Series(item))

                    elif feature.type[-3:] == 'RNA':
                        item['ID'] = item['GENE'] if 'GENE' in item else item['LOCUS_TAG']

                        rna_list.append(pd.Series(item))

                    elif feature.type in {'STS', 'misc_feature', 'rep_origin'}:
                        item['RECORD_TYPE'] = feature.type.upper()
                        misc_list.append(pd.Series(item))

        # Make tables
        df_gene = make_table(gene_list, '0,0,255')
        df_repeat = make_table(repeat_list, '255,0,255')
        df_mobile = make_table(mobile_list, '0,160,160')
        df_rna = make_table(rna_list, '255,0,0')
        df_misc = make_table(misc_list)

        # Write
        df_gene.to_csv(output.gene, sep='\t', index=False)
        df_repeat.to_csv(output.repeat, sep='\t', index=False)
        df_mobile.to_csv(output.mobile, sep='\t', index=False)
        df_rna.to_csv(output.rna, sep='\t', index=False)
        df_misc.to_csv(output.misc, sep='\t', index=False)
