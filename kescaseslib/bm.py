"""
Summarize benchmark tables.
"""

import pandas as pd


class RunTime:
    """
    Stores runtime (in seconds) output from the time utility.
    """

    def __init__(self, time_file_name):
        with open(time_file_name, 'r') as in_file:

            # Read times
            self.real = float(next(in_file).split(' ')[1])
            self.user = float(next(in_file).split(' ')[1])
            self.sys = float(next(in_file).split(' ')[1])

    def to_series(self, series_name):
        """
        Get a series from this runtime object.

        :param series_name: Name of this series (the step that generated the run-times).
        """

        time_series = pd.Series({
            'real': self.real,
            'user': self.user,
            'sys': self.sys,
        })

        time_series.name = series_name

        return time_series

class MaxMemUsage:
    """
    Stores the maximum memory usage per trace.
    """

    def __init__(self, mem_file_name):

        df = pd.read_table(mem_file_name)

        self.max_rss = max(df['rss_kb'])
        self.max_vss = max(df['vss_kb'])

    def to_series(self, series_name):
        """
        Get a series from this memory usage object.

        :param series_name: Name of this series (the step that consumed the memory).
        """

        mem_series = pd.Series({
            'rss': self.max_rss,
            'vss': self.max_vss
        })

        mem_series.name = series_name

        return mem_series

def summarize_runtime_table(tab_file_name, accession):
    """
    Generate a series by adding all the runtimes for each type (real, sys, user).

    :param tab_file_name: Name of the table file with runtimes for each step.
    :param accession: Sample accession.
    """

    summary_series = pd.read_table(tab_file_name, header=0, index_col=0).apply(sum, axis=1)

    summary_series.name = accession

    return summary_series


def summarize_trace_table(tab_file_name, accession):
    """
    Get a series describing the maximum vss and rss from a table of rss and vss for each step.

    :param tab_file_name: Table of rss and vss for each step.
    :param accession: Sample accession.
    """

    summary_series = pd.read_table(tab_file_name, header=0, index_col=0).apply(max, axis=1)

    summary_series.name = accession

    return summary_series
