#!/usr/bin/python3.5
import argparse


def main(sampleFile, pos, sep=''):
    """ script to use counterFeaturesInHeaders """
    seqInfo = bitk3.fastaReader(sampleFile)
    if sep == '':
        countsDict = bitk3.countFeaturesInHeaders(seqInfo[0], pos)
    else:
        countsDict = bitk3.countFeaturesInHeaders(seqInfo[0], pos, sep)

    keys = list(countsDict.keys())
    keys.sort()

    for key in keys:
        print('{} : {}'.format(key, countsDict[key]))

    return countsDict

if __name__ == "__main__":
    import bitk3
    parser = argparse.ArgumentParser(prog='countFeatureInTag', usage='%(prog)s fasta_file.fa position', description='Count features in headers of fasta formated files')
    parser.add_argument('sampleFile', metavar='fasta_file.fa', type=str, help='Valid fasta formated file')
    parser.add_argument('pos', metavar='position', type=int, help='Integer of the position of the information to be counted. Positions are divided by default using "|" as separator')
    parser.add_argument('--sep', metavar='sep', type=str, help='Change the character that separates the fields in header', default='')

    args = parser.parse_args()

    sampleFile = args.sampleFile
    pos = args.pos
    sep = args.sep

    main(sampleFile, pos, sep)
else:
    from bitk3 import bitk3