class Interval:
    """
    Defines an interval of the reference genome.

    :param chrom: Chromosome.
    :param name: Region name.
    :param start: Start base (1-based, inclusive).
    :param end: End base (1-based, inclusive).
    """

    def __init__(self, chrom, name, start, end):

        if name is None:
            name = '{}-{}-{}'.format(chrom, start, end)

        self.chrom = chrom
        self.name = name
        self.start = int(start)
        self.end = int(end)

    def is_in_interval(self, chrom, locus):
        """
        Determine if a location is in this interval.

        :param chrom: Chromosome name.
        :param locus: Location.

        :return: True if the location is in this interval, and False if it is not.
        """
        return chrom == chrom and self.start <= locus <= self.end

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

    def add_interval(self, chrom, name, start, end):
        """
        Add an interval to this container.

        :param chrom: Chromosome name.
        :param name: Interval name or ``None`` to build a default name.
        :param start: Start position (1-based, inclusive).
        :param end: end position (1-based, inclusive).
        """

        interval = Interval(chrom, name, start, end)

        self.interval_map[interval.name] = interval
        self.interval_list.append(interval.name)

    def get_interval(self, chrom, position):
        """
        Find an interval by the chromosome name and a position included within it.

        :param chrom: Chromosome name.
        :param position: Position within the interval.

        :return: An Interval object or None if no matching interval was found.
        """

        for interval in self.interval_map.values():
            if interval.chrom == chrom and interval.start <= position <= interval.end:
                return interval

        return None

    def add_bed(self, bed_file_name):
        """
        Add intervals by BED file.

        :param bed_file_name: Name of the BED file to add.
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
                self.add_interval(chrom, name, start, end)

    def add_nc_tab(self, nc_tab_file, accession):
        """
        Read a no-call tab file and add any intervals that match `accession`

        The first line of the no-call tab file is a header, and it is ignored. Any lines beginning with `#`
        and blank lines are also ignored. The rest are tab-separated fields.

        Fields:
        1) Chromosome name.
        2) Accession or sample name.
        3) Start position (1-based, inclusive)
        4) End position (1-based, inclusive)
        5) Type: A keyword (such as "BADASM" for a bad assembly, "CONTAMINATION" for contaminated data, or
            "INCOMPLETE" if the sequence data is known to not cover the genome sufficiently).

        :param nc_tab_file: No-call tab file name.
        :param accession: Extract regions for this accession or sample name.
        """

        with open(nc_tab_file, 'r') as nc_file:

            line_number = 0

            # Read each line
            for line in nc_file:
                line_number += 1

                # Skip header
                if line_number == 1:
                    continue

                line = line.strip()

                # Skip empty and commented lines
                if not line or line.startswith('#'):
                    continue

                # Tokenize
                tok = line.split('\t')

                # Add region
                if tok[1] == accession:
                    self.add_interval(tok[0], '{0}-{2}-{4}'.format(*tok), int(tok[2]), int(tok[3]))

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
