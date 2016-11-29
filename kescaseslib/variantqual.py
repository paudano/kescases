"""
Utilities for processing the results of vcfeval.
"""

import gzip
import pandas as pd
import numpy as np

from kescaseslib import interval

def summarize_calls(tp_file_name, fp_file_name, fn_file_name, bs_file_name,
                    regions_file_name, nc_file_name=None,
                    accession=None
                    ):

    # Get regions
    regions = interval.IntervalContainer()
    regions.add_bed(regions_file_name)

    # Get no-call (removed) regions
    if nc_file_name is not None and accession is not None:
        nc_regions = interval.IntervalContainer()
        nc_regions.add_nc_tab(nc_file_name, accession)

    # Init dataframe
    summary_df = pd.DataFrame(
        index=regions.interval_list,
        columns=('BL', 'TP', 'FP', 'FN', 'TPR', 'PPV', 'F1')
    )

    summary_df['TPR'] = summary_df['TPR'].astype(np.float)
    summary_df['PPV'] = summary_df['PPV'].astype(np.float)
    summary_df['F1'] = summary_df['F1'].astype(np.float)

    # Get baseline TP calls (actual calls)
    with gzip.open('bs_file_name', 'rb') as tp_file:

        line_number = 0

        # Process each line
        for line in tp_file:
            line_number += 1

            line = line.strip

            if not line or line.startswith('#'):
                continue




        interval_name = regions.get_interval()