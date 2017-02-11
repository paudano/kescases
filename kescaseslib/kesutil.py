"""
Miscellaneous utility functions.
"""

import gzip


class BinFileDecoder:

    def __init__(self, in_file_name):
        """
        Create a new decoder.

        :param in_file_name: Input file name.
        """

        self.in_file_name = in_file_name
        self.in_file = None
        self.do_encode = False

    def __enter__(self):
        """
        Open input file for reading.

        :return: `self`.
        """

        if self.in_file_name.endswith('.gz'):
            self.in_file = gzip.open(self.in_file_name, 'rb')
            self.do_encode = True
        else:
            self.in_file = open(self.in_file_name, 'r')
            self.do_encode = False

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Close open files.

        :param exc_type: Ignored.
        :param exc_value: Ignored.
        :param traceback: Ignored.
        """
        if self.in_file is not None:
            self.in_file.close()
            self.in_file = None
            self.do_encode = False

    def __iter__(self):
        """
        Get a reference to this object.

        :return: `self`.
        """
        return self

    def __next__(self):
        """
        Get the next line or raise `StopIteration`.

        :return: Next line.
        """

        if self.do_encode:
            out_str = self.in_file.readline().decode('utf-8').strip()
        else:
            out_str = self.in_file.readline().strip()

        if not out_str:
            raise StopIteration()

        return out_str


def has_uncommented_lines(in_file_name):
    """
    Determine if a file has uncommnted lines.

    :param in_file_name: Input file name. If the file ends with ".gz", it is read as a compressed text file.

    :return: `True` if the file has uncommented non-blank lines (whitespace only), and `False` otherwise.
    """

    with BinFileDecoder(in_file_name) as in_file:
        for line in in_file:
            if line != '' and not line.startswith('#'):
                return True

    return False
