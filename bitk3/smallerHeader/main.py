#!/usr/bin/python3.5
import argparse
from smallerHeader import SmallerHeader
import bitk3
import json

if __name__ == "__main__":
    if __package__ is None:
        from os import sys, path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    parser = argparse.ArgumentParser(
        prog='smallHeader',
        usage='%(prog)s fastaFile.fa',
        description='Switches all headers of FASTA file  \
        to XXNXX where N is the sequence number. \
        It also outputs a file to restore the original name with restoreHeader',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'fastafile',
        metavar='fastafile.fa',
        type=argparse.FileType('r'),
        help='FASTA formated file'
    )
    args = parser.parse_args()
    fastaFile = str(args.fastafile.name)

    fasta = SmallerHeader(fastaFile)

    outputFile = bitk3.insertMethod(fastaFile, 'smallerHeader')
    with open(outputFile, 'w') as f:
        newFasta = str(fasta)
        f.write(newFasta)

    dictionaryFile = bitk3.changeExtension(outputFile, 'json')
    with open(dictionaryFile, 'w') as f:
        json.dump(fasta.transdic, f, indent=2)
