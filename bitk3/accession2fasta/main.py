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
        description='Make FASTA formated file from accession using bitk3tags'
    )
    parser.add_argument(
        'accessionList',
        metavar='accessionList.txt',
        type=argparse.FileType('r'),
        help='Valid fasta formated file'
    )
    args = parser.parse_args()

    accessionList = args.accessionList

    ac = Accession2fasta(accessionList)

    with open('accession2fasta.fa', 'w') as f:
        f.write(ac.fasta)