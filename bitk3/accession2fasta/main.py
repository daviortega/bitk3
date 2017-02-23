#!/usr/bin/python3.5
import argparse
from accession2fasta import Accession2fasta

if __name__ == "__main__":
    if __package__ is None:
        from os import sys, path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    parser = argparse.ArgumentParser(
        prog='accession2fasta',
        usage='%(prog)s accessionList.txt',
        description='Make FASTA formated file from accession using bitk3tags',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'accessionList',
        metavar='accessionList.txt',
        type=argparse.FileType('r'),
        help='File with single refseq accession number per line'
    )
    parser.add_argument(
        '-o',
        '--outfile',
        metavar='filename.fa',
        type=str,
        help='Optional name for output.',
        default='accession2fasta.fa'
    )
    args = parser.parse_args()

    accessionList = args.accessionList
    outputfile = args.outfile

    ac = Accession2fasta(accessionList)

    with open(outputfile, 'w') as f:
        f.write(ac.fasta)