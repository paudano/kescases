import pandas as pd

class Interval:
    """
    Defines an interval of the reference genome.

    :param chrom: Chromosome.
    :param name: Region name.
    :param start: Start base (1-based, inclusive).
    :param end: End base (1-based, inclusive).
    :param tag: Optional tag for the interval. For example, if the interval is a filtered region,
        the tag may contain information about why it was filtered.
    """

    def __init__(self, chrom, name, start, end, tag=None):

        if name is None:
            name = '{}-{}-{}'.format(chrom, start, end)

        self.chrom = chrom
        self.name = name
        self.start = int(start)
        self.end = int(end)
        self.tag = tag

    def is_in_interval(self, chrom, locus):
        """
        Determine if a location is in this interval.

        :param chrom: Chromosome name.
        :param locus: Location.

        :return: True if the location is in this interval, and False if it is not.
        """
        return chrom == chrom and self.start <= locus <= self.end

    def range_in_interval(self, chrom, start, end):
        """
        Determine if a range of bases overlaps this interval by at least one base.

        :param chrom: Chomosome.
        :param start: Start (1-based, inclusive)
        :param end: End (1-based, inclusive)

        :return: `True` if the range is in this interval.
        """

        if chrom != self.chrom:
            return False

        return start <= self.end and end >= self.start

    def __repr__(self):
        """
        Get a string representation of this interval.

        :return: String representation of this interval.
        """
        return 'Interval[name={name}, chrom={chrom}, loc={start}-{end}]'.format(**self.__dict__)

    def __str__(self):
        return self.name

    def __gt__(self, interval):
        """
        Test greater-than.

        :param interval: Interval.

        :return: Boolean result of test.
        """
        if self.chrom != interval.chrom:
            return self.chrom > interval.chrom

        if self.start != interval.start:
            return self.start > interval.start

        return self.end > interval.end

    def __lt__(self, interval):
        """
        Test less-than.

        :param interval: Interval.

        :return: Boolean result of test.
        """
        if self.chrom != interval.chrom:
            return self.chrom < interval.chrom

        if self.start != interval.start:
            return self.start < interval.start

        return self.end < interval.end

    def __ge__(self, interval):
        """
        Test greater-than or equal to.

        :param interval: Interval.

        :return: Boolean result of test.
        """
        if self.chrom != interval.chrom:
            return self.chrom > interval.chrom

        if self.start != interval.start:
            return self.start > interval.start

        return self.end >= interval.end

    def __le__(self, interval):
        """
        Test less-than or equal to.

        :param interval: Interval.

        :return: Boolean result of test.
        """
        if self.chrom != interval.chrom:
            return self.chrom < interval.chrom

        if self.start != interval.start:
            return self.start < interval.start

        return self.end <= interval.end

    def __eq__(self, interval):
        """
        Test equal-to.

        :param interval: Interval.

        :return: Boolean result of test.
        """
        return self.chrom == interval.chrom and self.start == interval.start and self.end == interval.end

    def __ne__(self, interval):
        """
        Test not equal-to.

        :param interval: Interval.

        :return: Boolean result of test.
        """
        return self.chrom != interval.chrom or self.start != interval.start or self.end != interval.end


class IntervalContainer:
    """
    A set of intervals.
    """

    def __init__(self):
        self.interval_map = dict()
        self.interval_list = []

    def add_interval(self, chrom, name, start, end, tag=None):
        """
        Add an interval to this container.

        :param chrom: Chromosome name.
        :param name: Interval name or ``None`` to build a default name.
        :param start: Start position (1-based, inclusive).
        :param end: end position (1-based, inclusive).
        :param tag: Optional tag describing this interval.
        """

        interval = Interval(chrom, name, start, end, tag)

        self.interval_map[interval.name] = interval
        self.interval_list.append(interval.name)

    def get_interval(self, chrom, pos, end=None):
        """
        Find an interval by the chromosome name and a position included within it.

        :param chrom: Chromosome name.
        :param pos: Position within the interval.
        :param end: End position to search, or `None` search for a single base at `pos`.

        :return: An Interval object or `None` if no matching interval was found.
        """

        if end is None:
            end = pos

        for interval_name in self.interval_list:

            interval = self.interval_map[interval_name]

            if interval.chrom != chrom:
                continue

            if pos <= interval.end and end >= interval.start:
                return interval

        return None

    def add_bed(self, bed_file_name, tag=None):
        """
        Add intervals by BED file.

        :param bed_file_name: Name of the BED file to add.
        :param tag: Optional tag describing the intervals in this BED file.
        """

        line_count = 0

        with open(bed_file_name, 'r') as bed_file:
            for line in bed_file:
                line_count += 1

                # Check line
                line = line.strip()

                if line == '' or line.startswith('#') or line.startswith('browser') or line.startswith('track'):
                    continue

                tok = line.split('\t')

                if len(tok) < 3:
                    raise RuntimeError('Bed file contains fewer than 3 tab-delimited fields on line {}: {}'
                                       .format(line_count, line))

                # Get fields
                chrom = tok[0]
                start = int(tok[1]) + 1
                end = int(tok[2])

                if len(tok) >= 4:
                    name = tok[3]
                else:
                    name = '{}-{}-{}'.format(chrom, start, end)

                # Add entry
                self.add_interval(chrom, name, start, end, tag)

    def add_blacklist(self, blacklist_file_name, accession):
        """
        Read a blacklist tab file and add any intervals that match `accession`.

        The first line of the blacklist file is a header. Any lines beginning with `#`
        and blank lines are also ignored. The rest are tab-separated fields.

        Fields:
        1) chr: Chromosome name.
        2) accession: Accession or sample name.
        3) start: Start position (1-based, inclusive)
        4) end: End position (1-based, inclusive)
        5) type: A keyword (such as "BADASM" for a bad assembly, "CONTAMINATION" for contaminated data, or
            "INCOMPLETE" if the sequence data is known to not cover the genome sufficiently). This becomes
            the filter tag.

        Fields beyond these (such as the "reason" field in the blacklist file) are not used.

        :param blacklist_file_name: Blacklist file name.
        :param accession: Extract regions for this accession.
        """

        # Read blacklist entries for this accession
        df_bl = pd.read_table(blacklist_file_name, header=0, skip_blank_lines=True, comment='#')
        df_bl = df_bl.ix[df_bl['accession'] == accession, :]

        if df_bl.shape[0] > 0:
            for row in df_bl.iterrows():
                row = row[1]

                self.add_interval(row['chr'], None, row['start'], row['end'], row['type'])

    def get_interval_count(self):
        """
        Get the number of intervals in this container.

        :return: The number of intervals in this container.
        """
        return len(self.interval_list)

    def write_bed(self, bed_file_name):
        """
        Write intervals as a BED4 file.

        :param bed_file_name: Name of the BED file.
        """

        with open(bed_file_name, 'w') as bed_file:
            for interval in sorted(self.interval_map.values()):
                bed_file.write('{}\t{}\t{}\t{}\n'.format(interval.chrom, interval.start - 1, interval.end, interval.name))
