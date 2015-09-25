"""
This filter removes bootstrap values from a newick tree.

"""
from __future__ import print_function, division

import sys

import dendropy


def main():
    tree = dendropy.Tree.get(file=sys.stdin, schema='newick')
    tree.write(
            schema='newick',
            file=sys.stdout,
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True,
            suppress_annotations=True,
            suppress_item_comments=True,
            )


if __name__ == '__main__':
    main()
