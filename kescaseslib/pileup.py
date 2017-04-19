"""
Tools for processing SAMtools pileup output.
"""

import re
from kescaseslib import variant
from kescaseslib import interval


def read_pileup_variants(pileup_file_name, sample_name, interval_container, no_call_bed, allow_n=True):
    """
    Read variants from the a pileup file.

    :param pileup_file_name: Assembled-aligned pileup file.
    :param sample_name: Name of the sample variants belong to.
    :param interval_container: Interval container.
    :param no_call_bed: BED file of regions where variants could not be called (BED4) or `None` if no
        file should be written.
    :param allow_n: Allow a variant to contain "N" in the reference or alternate sequence if `True`.

    :return: A list of variants.
    """

    # Create a list of variants and no-call regions
    var_list = list()

    no_call = interval.IntervalContainer()

    # Initialize
    chr_name = ''
    loc = 0

    read_base = ''
    read_base_org = ''
    ref = ''
    alt = ''

    line_count = 0

    # Initialize reference container and related variables
    current_interval = None
    ref_sequence = ''
    alt_sequence = ''
    depth_list = []

    nc_start = 0  # In a no-call interval when greater than 0 (assembly coverage is not 1)

    # Open file
    with open(pileup_file_name, 'r') as pileup_file:

        # Read all variants
        while True:

            # Get line if read_base is empty
            if not read_base:
                line = pileup_file.readline()

                if not line:  # EOF
                    break

                line = line.strip()

                line_count += 1

                if not line:
                    continue  # Blank line

                # Tokenize
                tok = line.strip().split('\t')

                if len(tok) < 5:
                    raise RuntimeError('Pileup line does not have at least 5 tab-delimited columns ({} at line {})'
                                       .format(pileup_file_name, line_count))

                # Get fields
                chr_name = tok[0]
                loc = int(tok[1])
                ref = tok[2].upper()
                depth = int(tok[3])
                read_base = tok[4].upper()

                read_base_org = read_base

                # Check for a new interval
                if current_interval is None or loc > current_interval.end:

                    # Add current sequence to the reference container
                    if current_interval is not None:

                        # Add no-call regions
                        if nc_start > 0:
                            no_call.add_interval(current_interval.chrom, None,
                                                 nc_start, current_interval.end)

                        # Reset
                        ref_sequence = ''
                        alt_sequence = ''
                        depth_list = []
                        nc_start = 0

                    # Get next interval
                    current_interval = interval_container.get_interval(chr_name, loc)

                    if current_interval is None:
                        raise RuntimeError('Variant in pileup file is not in a configured interval '
                                           '({} at line {}): {}:{} {}->{}'
                                           .format(pileup_file_name, line_count, chr_name, loc, ref, alt))

                # Append reference sequence
                ref_sequence += ref
                depth_list.append(depth)

                # Check coverage
                if nc_start > 0:
                    alt_sequence += ref

                    if depth == 1:
                        no_call.add_interval(current_interval.chrom, None, nc_start, loc - 1)
                        nc_start = 0

                    else:
                        read_base = ''
                        continue

                # Depth must be uniform 1
                if depth != 1:
                    alt_sequence += ref
                    nc_start = loc
                    read_base = ''
                    continue

                # Remove read start-end markers
                read_base = re.sub('[.,$]|(\^.)', '', read_base)

                # Check for deletion (handled in a previous line)
                if read_base == '*':
                    read_base = ''
                    continue

                # Stop if there is no variant at this locus
                if not read_base:
                    alt_sequence += ref
                    continue

            # Convert to variant
            if re.match('[ACGTN]', read_base):
                alt = read_base[0]

                if len(alt) > 1:
                    read_base = read_base[1:]
                else:
                    read_base = ''

            elif re.match('\+\d+[ACGTN]+', read_base):
                loc += 1
                ref = ''
                alt = re.sub('\+\d+([ACGTN]+)', '\\1', read_base)
                read_base = ''

            elif re.match('-\d+[ACGTN]+', read_base):
                loc += 1
                ref = re.sub('-\d+([ACGTN]+)', '\\1', read_base)
                alt = ''
                read_base = ''

            else:
                raise RuntimeError('Pileup line contains an unexpected read-base column ({} at line {}): {}'
                                   .format(pileup_file_name, line_count, read_base_org))

            # Skip N in reference if disabled.
            if not allow_n and ('N' in ref.upper() or 'N' in alt.upper()):
                continue

            # Append variant
            var_list.append(variant.Variant(
                chr_name, loc, ref, alt, current_interval, sample_name, var_depth=1, loc_depth=1
            ))

            # Add to alternate sequence
            alt_sequence += alt

        # Write last no-call
        if current_interval is not None and nc_start > 0:
            no_call.add_interval(current_interval.chrom, None, nc_start, loc)

        # Write BED file
        if no_call_bed is not None:
            no_call.write_bed(no_call_bed)

        # Return variants
        return var_list
