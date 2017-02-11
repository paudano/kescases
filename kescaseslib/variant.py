"""
Data structures for representing variants.
"""

import pysam


class Variant:
    """
    Represents a variant.

    :param chrom: Chromosome name.
    :param start: Start location on the chromosome (1-based, inclusive).
    :param ref: Reference allele.
    :param alt: Alternate allele.
    :param interval: Interval this variant belongs to or `None`.
    :param sample_name: Name of the sample this variant belongs to.
    :param format: Format of the call field if this variant was in a VCF file.
    :param call: Call field if this variant was in a VCF file (last column; variant call string complying to FORMAT).
    :param filter: Filter if this variant passed some quality filter(s).
    :param var_depth: Estimated depth of reads supporting this variant.
    :param loc_depth: Estimated read depth at this locus.
    :param quality: Variant quality score or `None`.
    :param include_pad: Insertion and deletion variants include pad bases. If `False`, then `ref` or
        `alt` is an empty string for an insertion or a deletion (respectively). If `True`, `ref` or
        `alt` has more than one base. `start` is the 1-based inclusive position of the first base
        in `ref` regardless of the variant type.
    """

    def __init__(self, chrom, start, ref, alt, interval,
                 sample_name, format=None, call=None, filter=None,
                 var_depth=None, loc_depth=None, quality=None, include_pad=False):

        self.chrom = chrom
        self.start = int(start)
        self.ref = ref.upper()
        self.alt = alt.upper()
        self.interval = interval
        self.sample_name = sample_name
        self.format = format
        self.call = call
        self.filter = filter
        self.var_depth = int(var_depth) if var_depth is not None else None
        self.loc_depth = int(loc_depth) if var_depth is not None else None
        self.quality = float(quality) if quality is not None else None

        # Update ref and alt if a pad base was included
        if include_pad and (len(self.ref) > 1 or len(self.alt) > 1):

            # Both ref and alt cannot be more than one base (only one can be)
            if len(self.ref) > 1 and len(self.alt) > 1:
                raise ValueError(
                    (
                        'Reference and alternate alleles are both longer than one base: '
                        'ref="{}", alt="{}", chrom="{}", start="{}"'
                    ).format(self.ref, self.alt, self.chrom, self.start)
                )

            # Remove pad
            if self.ref[0] == self.alt[0]:
                self.ref = self.ref[1:]
                self.alt = self.alt[1:]
                self.start += 1

            elif self.ref[-1] == self.alt[-1]:
                self.ref = self.ref[:-1]
                self.alt = self.alt[:-1]

            else:
                raise ValueError(
                    (
                        'Variant is an indel with padding, but ref and alt alleles '
                        'do not start or do not end with the same base: '
                        'ref="{}", alt="{}", chrom="{}", start="{}"'
                    ).format(self.ref, self.alt, self.chrom, self.start)
                )

    def is_snp(self):
        """
        Determine if this variant is a SNP.

        :return: ``True`` if this variant is a SNP.
        """
        return len(self.ref) == 1 and len(self.alt) == 1

    def is_ins(self):
        """
        Determine if this variant is an insertion.

        :return: ``True`` if this variant is an insertion.
        """
        return len(self.ref) == 0 and len(self.alt) > 0

    def is_del(self):
        """
        Determine if this variant is a deletion.

        :return: ``True`` if this variant is a deletion.
        """
        return len(self.ref) > 0 and len(self.alt) == 0

    def length(self):
        """
        Get the length of this variant.

        :return: Number of bases affected by this variant. 1 for SNPs, and length of insertion or deletion for indels.
        """

        if self.is_snp():
            return 1

        if self.is_ins():
            return len(self.alt)

        return len(self.ref)

    def get_end(self):
        """
        Get the last reference base affected by this variant.

        :return: Last reference base affected by this variant.
        """

        if self.is_snp() or self.is_ins():
            return self.start

        return self.start + self.length() - 1

    def compare(self, other):
        """
        Compare this variant to another.

        :param other: Other variant.

        :return: a negative number, zero, or a positive number if this variant is less than, equal to, or greater
            than the other variant (respectively).
        """

        if self.chrom != other.chrom:
            return -1 if self.chrom < other.chrom else 1

        if self.start != other.start:
            return -1 if self.start < other.start else 1

        if self.ref != other.ref:
            return -1 if self.ref < other.ref else 1

        if self.alt != other.alt:
            return -1 if self.alt < other.alt else 1

        return 0

    def __lt__(self, other):
        """
        Determine if this variant is less than another.

        :param other: Other variant.

        :return: True if this variant is less than the other.
        """
        return self.compare(other) < 0

    def __le__(self, other):
        """
        Determine if this variant is less than or equal to another.

        :param other: Other variant.

        :return: True if this variant is less than or equal to  the other.
        """
        return self.compare(other) <= 0

    def __gt__(self, other):
        """
        Determine if this variant is greater than another.

        :param other: Other variant.

        :return: True if this variant is greater than the other.
        """
        return self.compare(other) > 0

    def __ge__(self, other):
        """
        Determine if this variant is greater than or equal to another.

        :param other: Other variant.

        :return: True if this variant is greater than or equal to the other.
        """
        return self.compare(other) >= 0

    def __eq__(self, other):
        """
        Determine if this variant is equal to another.

        :param other: Other variant.

        :return: True if this variant is equal to the other.
        """

        if other is None:
            return False

        return self.compare(other) == 0

    def __ne__(self, other):
        """
        Determine if this variant is not equal to another.

        :param other: Other variant.

        :return: True if this variant is not equal to the other.
        """

        if other is None:
            return True

        return self.compare(other) != 0

    def __str__(self):
        """
        Get a string representation of this variont.

        :return: String representation of this variant.
        """

        if len(self.ref) == 1 and len(self.alt) == 1:
            return '{start}{ref}>{alt}'.format(**self.__dict__)

        if len(self.ref) == 1:
            return '{start}ins{alt}'.format(**self.__dict__)

        if len(self.alt) == 1:
            return '{}-{}del'.format(self.start, self.start + len(self.ref) - 1)

        return 'UNKNOWN_VAR[start={start}, ref={ref}, alt={alt}]'.format(**self.__dict__)


def write_vcf(variant_list, sample_name, vcf_file_name, ref_file_name):
    """
    Write a VCF file.

    :param variant_list: List of variants or a single variant.
    :param sample_name: Name of this sample.
    :param vcf_file_name: VCF file to write.
    :param ref_file_name: Reference FASTA file. There must be a file of the same name with ".fai" appended (samtools
        index).
    """

    variant_list = sorted(list(variant_list))

    if sample_name is None:
        sample_name = 'UNKNOWN'

    fasta_file = pysam.FastaFile(ref_file_name)

    with open(vcf_file_name, 'w') as out_file:

        # Start header
        out_file.write('##fileformat=VCF4.2\n')
        out_file.write('##source=kescases\n')

        # Write contigs
        with open(ref_file_name + '.fai', 'r') as fai_file:

            line_number = 0

            for line in fai_file:
                line_number += 1

                line = line.strip()

                if not line:
                    continue

                tok = line.split('\t')

                if len(tok) < 2:
                    raise ValueError(
                        'Index file {} contains a record on line {} that does not have at least two fields'
                        .format(ref_file_name + '.fai', line_number)
                    )

                out_file.write('##contig=<ID={},length={}>\n'.format(tok[0], tok[1]))

        # Write FORMAT headers
        out_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        out_file.write('##FORMAT=<ID=GDP,Number=A,Type=Integer,Description='
                       '"Estimated depth of all haplotypes supporting the alternate variant">\n'
                       )

        out_file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description='
                       '"Estimated depth of all haplotypes in the variant active region">\n'
                       )

        # Write column heading
        out_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(sample_name))

        # Write variants
        for variant in variant_list:

            # Get record elements
            if variant.is_snp():
                ref = variant.ref
                alt = variant.alt
                start = variant.start

            elif variant.is_ins():
                if variant.start > 1:
                    start = variant.start - 1
                    ref = fasta_file.fetch(variant.chrom, start - 1, start)
                    alt = ref + variant.alt

                else:
                    start = variant.start
                    ref = fasta_file.fetch(variant.chrom, start - 1, start)
                    alt = variant.alt + ref

            elif variant.is_del():
                if variant.start > 1:
                    start = variant.start - 1
                    alt = fasta_file.fetch(variant.chrom, start - 1, start)
                    ref = alt + variant.ref

                else:
                    start = variant.start
                    end_loc = start + len(variant.ref)
                    alt = fasta_file.fetch(variant.chrom, end_loc - 1, end_loc)
                    ref = variant.ref + alt

            else:
                raise ValueError('Variant is an unrecognized type: {}'.format(str(variant)))

            # Write
            out_file.write('{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT:GDP:DP\t1:{}:{}\n'.format(
                           variant.chrom, start, ref, alt,
                           variant.var_depth, variant.loc_depth
                           ))
