#!/usr/bin/env python

import sys
import os
import screed


def main():
    '''
    pick seqs that are in list file
    '''
    if len(sys.argv) != 2:
        mes = ('*** Usage: python {} <seqfile.fa>\n')
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)

    seqfile = sys.argv[1]
    if seqfile == '-': 
        seqfile = '/dev/stdin'

    with screed.open(seqfile) as sp:
        for rec in sp:
            name = rec.name.split(None, 1)[0]
            seq = rec.sequence
            sys.stdout.write('{}\t{}\n'.format(name, len(seq)))

if __name__ == '__main__':
    main()
