import pandas as pd


def bed_interval_to_dataframe(bed_file_name):
    """
    Read a BED interval file, which is expected to be in BED6 format.

    Data frame fields:
    1. chrom: Chromosome name
    2. start: Start position (1-based, inclusive)
    3. end: End position in both 1-based inclusive (start) and 0-based exclusive (bed_start) coordinate systems.
    4. name: Name of this record (must be unique, same as row index)
    5. strand: + or -
    6. bed_start: Start position (0-based, inclusive)

    :param bed_file_name: Name of the BED6 file to read.

    :return: A Pandas dataframe with the fields described and row name set to the 'name' field.
    """

    # Read BED
    df = pd.read_table(bed_file_name, header=None)
    df.columns = ['chrom', 'bed_start', 'end', 'name', 'score', 'strand']

    # Set row names
    df = df.set_index(df['name'])

    # Translate to 1-based and inclusive (keep 0-based as bed_start, end is the same in eicher coordinate system)
    df['start'] = df['bed_start'] + 1

    # Rearrange columns (and drop the score)
    df = df.ix[:, ('chrom', 'start', 'end', 'name', 'strand', 'bed_start')]

    # Return data frame
    return df
