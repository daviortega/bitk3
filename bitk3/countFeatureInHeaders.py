#!/usr/bin/python3.5
from bitk3 import bitk3
import sys
import argparse

def main(sampleFile, pos):
    """ script to use counterFeaturesInHeaders """
    seqInfo = bitk3.fastaReader(sampleFile)
    countsDict = bitk3.countFeaturesInHeaders(seqInfo[0], pos)

    keys = list(countsDict.keys())
    keys.sort()

    for key in keys:
        print('{} : {}'.format(key, countsDict[key]))

    return countsDict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='countFeatureInTag', usage='%(prog)s fasta_file.fa position', description='Count features in headers of fasta formated files')
    parser.add_argument('sampleFile', metavar='fasta_file.fa', type=str, help='Valid fasta formated file')
    parser.add_argument('pos', metavar='position', type=int, help='Integer of the position of the information to be counted. Positions are divided by default using "|" as separator')
    args = parser.parse_args()

    sampleFile = args.sampleFile
    pos = args.pos

    main(sampleFile, pos)