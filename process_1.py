"""
Filter a Phylip-like alignment.

Columns are preserved if and only they meet these conditions:
1) The column is annotated with an asterisk (high confidence column).
2) Each character in the column is an unambiguous nucleotide (that is,
   not a gap character or a nucleotide ambiguity character).

The input consists of interleaved phylip paragraphs
each with the following format:
1) A position number indicating the alignment position of the first
   column of the paragraph.
2) For each taxon, the name of the taxon followed by up to 10x10
   alignment characters
3) An annotation line, where each position is annotated
   with either an asterisk (sites intended to be correctly aligned)
   or a space (sites not necessarily correctly aligned).

The output is a sequential (not interleaved) phylip alignment
in the same format as the example alignment on the phyml web page.

                      1
Mycobacter.lepra      ---------- ---------- ---------- ---------- ---------- ---------- --AAGTAAGT G-CCTAAGGG CGCATGGTGG ATGCCTT--G
Streptomyc.grise      ---------- ---------- ---------- ---------- ---------- ---------- --GGCCAAGT T-TTTAAGGG CGCACGGTGG ATGCCTT--G
[...]

"""
from __future__ import print_function, division


from StringIO import StringIO
import sys
import argparse
from itertools import izip_longest

from numpy.testing import assert_equal, assert_

PHYLIP_SCHEMA = 'phylip'
MASE_SCHEMA = 'mase'


def get_offset_of_first_nonspace_following_space(line):
    n = len(line)
    for i in range(n-1):
        if line[i] == ' ' and line[i+1] != ' ':
            return i+1
    return None


class Para:

    def __init__(self):
        self.first_pos = None
        self.name_indices = None
        self.seq_indices = None
        self.name_seq_pairs = []
        self.annotation = None

    def add_pos(self, first_pos):
        self.first_pos = first_pos

    def add_data_line(self, line):
        # nnnnnnnn______________ssssssssss_ssssssssss_ssssssssss_...

        # Find the offset of the first nonspace character that follows a space.
        k = get_offset_of_first_nonspace_following_space(line)
        if k is None:
            raise ValueError('failed to find name boundary in line "%s"' % line)

        # Find offsets of all name characters.
        # This includes offsets of potential trailing spaces.
        name_indices = range(k)

        # Find the offsets of sequence or annotation indices.
        seq_indices = [i for i, c in enumerate(line) if i >= k and c != ' ']

        # Check that the name indices and seq indices
        # are consistent within a paragraph.
        if self.name_seq_pairs:
            assert_equal(name_indices, self.name_indices)
            assert_equal(seq_indices, self.seq_indices)
        else:
            self.name_indices = name_indices
            self.seq_indices = seq_indices

        name = ''.join(line[i] for i in name_indices).strip()
        seq = ''.join(line[i] for i in seq_indices)
        name_seq_pair = (name, seq)
        self.name_seq_pairs.append(name_seq_pair)

    def add_annotation_line(self, line):
        n = len(line)
        annotation = []
        for i in self.seq_indices:
            if i < n:
                annotation.append(line[i])
            else:
                annotation.append(' ')
        self.annotation = ''.join(annotation)

    def get_annotated_gapfree_count(self):
        return len(self.get_annotated_gapfree_set())

    def get_annotated_gapfree_set(self):
        n = len(self.annotation)
        annotated_set = set(i for i, c in enumerate(self.annotation) if c == '*')
        annotated_gapfree_set = annotated_set.copy()
        for name, seq in self.name_seq_pairs:
            gap_indices = set(i for i, c in enumerate(seq) if c not in 'ACGT')
            annotated_gapfree_set -= gap_indices
        return annotated_gapfree_set

    def __str__(self):
        f = StringIO()
        print('name indices:', self.name_indices, file=f)
        print('seq indices:', self.seq_indices, file=f)
        print('annotation:', self.annotation, file=f)
        print('name sequence pairs:', file=f)
        print(self.name_seq_pairs, file=f)
        return f.getvalue()


def gen_paragraphs(fin):
    p = None
    for line in fin:
        stripped = line.strip()

        # Look for a position number to begin a new paragraph.
        if p is None:

            # If the stripped line is empty
            # then continue looking for a position number.
            if not stripped:
                continue

            # Check if the stripped line consists of an integer.
            # If so, this is the position number of a new paragraph.
            pos = None
            try:
                pos = int(stripped)
            except TypeError as e:
                pass
            if pos is not None:
                #print('creating paragraph for initial position', pos)
                p = Para()
                p.add_pos(pos)
                continue

        # Determine whether the line defines a sequence
        # or whether the line defines the paragraph annotation.
        elements = set(stripped)
        #print('elements:', elements)
        if elements <= {' ', '*'}:
            p.add_annotation_line(line.rstrip())
            yield p
            p = None
        else:
            p.add_data_line(line.rstrip())

    # Assert that all paragraphs have been yielded,
    # and that we are not in the process of constructing a paragraph.
    assert_(p is None)


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def meld(paragraphs, fout, schema):
    """
    The output will look like the following, not interleaved.
    The two numbers in the header line are the number of taxa
    and the number of columns respectively.

    phylip schema:
     60 500
    T25
    ACTATTGAAAGAAGGGGGTTCCTAGATATCTGCGAGTATAATCGTGCTTGGTCTCCTATC...
    T16
    ATTAATCAAAGTAGGCGGGGCGGCCGTAGATGCTAAGAAAATCGAGTTCGGTCACCTCCC...

    mase schema:
    ;; this alignment is in the 'mase' format
    ;comment line
    T25
    ACTATTGAAAGAAGGGGGTTCCTAGATATCTGCGAGTATAATCGTGCTTGGTCTCCTATC...
    ;comment line
    T16
    ATTAATCAAAGTAGGCGGGGCGGCCGTAGATGCTAAGAAAATCGAGTTCGGTCACCTCCC...


    """
    assert_(schema in (PHYLIP_SCHEMA, MASE_SCHEMA))
    name_count = None
    out_pairs = None
    for p in paragraphs:
        if name_count is None:
            name_count = len(p.name_seq_pairs)
        if out_pairs is None:
            out_pairs = [(name, []) for name, seq in p.name_seq_pairs]
        indices = sorted(p.get_annotated_gapfree_set())
        for pair_idx, (name, seq) in enumerate(p.name_seq_pairs):
            for i in indices:
                out_pairs[pair_idx][1].append(seq[i])

    if schema == PHYLIP_SCHEMA:
        s = ' {} {}'.format(name_count, len(out_pairs[0][1]))
        print(s, file=fout)
        for name, nucs in out_pairs:
            print(name, file=fout)
            print(''.join(nucs), file=fout)

    if schema == MASE_SCHEMA:
        print(";; this alignment is in the 'mase' format", file=fout)
        for name, nucs in out_pairs:
            print('; no description', file=fout)
            print(name, file=fout)
            for chunk in grouper(nucs, 60, fillvalue=''):
                print(''.join(chunk), file=fout)


def main(args, fin):
    """
    """
    paragraphs = []
    for p in gen_paragraphs(fin):
        paragraphs.append(p)

    if args.verbose:
        print('number of paragraphs:')
        print('first paragraph:')
        print(paragraphs[0])
        print('last paragraph:')
        print(paragraphs[-1])

        total = sum(len(p.seq_indices) for p in paragraphs)
        print('cumulative sequence length:', total)

        total = sum(p.annotation.count('*') for p in paragraphs)
        print('cumulative number of asterisks:', total)

        total = sum(p.get_annotated_gapfree_count() for p in paragraphs)
        print('cumulative count of annotated gapfree sites:', total)

        sets = [set(s) for p in paragraphs for n, s in p.name_seq_pairs]
        total = set.union(*sets)
        print('set of all sequence elements:', total)

    # write the phylip alignment if requested
    if args.phylip_out:
        with open(args.phylip_out, 'w') as fout:
            meld(paragraphs, fout, schema='phylip')

    # write the mase alignment if requested
    if args.mase_out:
        with open(args.mase_out, 'w') as fout:
            meld(paragraphs, fout, schema='mase')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--phylip-out', help='write this phylip alignment file')
    parser.add_argument('--mase-out', help='write this mase alignment file')
    args = parser.parse_args()
    main(args, sys.stdin)
