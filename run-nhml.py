"""
Run nhml possibly with multiple option files.

"""
from __future__ import print_function, division

import argparse
import subprocess
import os


def main(args):
    # output filenames replace the final .foo with .nhml

    # create the subprocess commands
    output_filenames = []
    subprocess_arg_lists = []
    for opt in args.optfiles:
        root, ext = os.path.splitext(opt)
        output_filenames.append(root + '.nhml')
        lst = [args.executable, args.alignment, args.tree, opt]
        subprocess_arg_lists.append(lst)

    # create the procs
    procs = []
    for lst in subprocess_arg_lists:
        proc = subprocess.Popen(
                lst, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        procs.append(proc)

    # get the output of each proc and write it to the output file
    stdouts = []
    for proc, output_filename in zip(procs, output_filenames):
        stdoutdata, stderrdata = proc.communicate()
        with open(output_filename, 'w') as fout:
            print(stdoutdata, file=fout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', required=True,
            help='the nhml executable')
    parser.add_argument('--alignment', required=True,
            help='the alignment file in mase format')
    parser.add_argument('--tree', required=True,
            help='the tree file in newick format')
    parser.add_argument('optfiles', nargs='+',
            help='one or more nhml option files')
    args = parser.parse_args()
    main(args)
