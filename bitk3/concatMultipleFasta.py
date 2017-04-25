#!/usr/bin/python3.5
import argparse


def _withSameNumberOfSequence(msaList=[]):
    first = True
    for seqDic, seqList in msaList:
        if first:
            msaLen = len(seqList)
            first = False
        else:
            if msaLen != len(seqList):
                raise Exception('The number of sequences \
in one of the files is different than the others')
    return True


def _sameOrganisms(msaList=[]):
    msaIndex = range(1, len(msaList))
    seqDic, seqList = msaList[0]
    for i, baseTag in enumerate(seqList):
        baseOrgID = baseTag.split(bitk3.BITKTAGSEP)[0]
        for j in msaIndex:
            tag = msaList[j][1][i]
            orgID = tag.split(bitk3.BITKTAGSEP)[0]
            if orgID != baseOrgID:
                print('The sequence {} in alignment number {} is \
different from alignment {}'.format(
                    tag,
                    j,
                    baseTag
                ))
                raise Exception('The alignments seems to be out of order')
    return True


def _buildPartFileRaxml(msaList=[]):
    """It will read the list of CCDs and build a partitioning file for RAxML
    Keyword arguments:
    msaList -- list of [dictionary, listOfheaders] msa.
    Returns string to be writen in file.
    """

    partFileRaxml = ''
    startCoord = 1
    partNum = 1

    for seqDic, seqList in msaList:
        endCoord = len(seqDic[seqList[0]]) - 1 + startCoord
        partFileRaxml += 'AUTO, part{} = {}-{}\n'.format(
            partNum,
            startCoord,
            endCoord
        )
        startCoord = endCoord + 1
        partNum += 1

    return partFileRaxml


def _parseData(dataSets=[]):

    output = {
        'partFileRaxml': '',
        'fastaConcat': '',
        'assocTable': [],
    }

    fastaString = ''
    assocTable = []

    for i, header in enumerate(dataSets[0][1]):
        fastaString += '>{}\n'.format(header)
        sequence = ''
        data = []
        for j in range(len(dataSets)):
            data.append(dataSets[j][1][i])
            sequence += dataSets[j][0][dataSets[j][1][i]]
        fastaString += '{}\n'.format(sequence)
        assocTable.append(data)

    partFileRaxml = _buildPartFileRaxml(dataSets)

    output['partFileRaxml'] = partFileRaxml
    output['fastaString'] = fastaString
    output['assocTable'] = assocTable

    return output


def main(fastaFilesHandles=[], noFiles=False):
    """
    It will build a concatenated multiple fasta

    Argument Keywords:
    fastaFilesHandles -- List of filenames

    Output:
    File with concatenated alignment of ABR
    JSON file with information

    Return:
    JSON with information.

    """

    dataSets = []
    for i, fastaFile in enumerate(fastaFilesHandles):
        seqDic, listOrder = bitk3.fastaReaderByHandle(fastaFile)
        if bitk3.isMSA(seqDic):
            dataSets.append([seqDic, listOrder])
        else:
            raise Exception('The {} file in the list is not an MSA.'.format(i))

    _sameOrganisms(dataSets)
    _withSameNumberOfSequence(dataSets)
    output = _parseData(dataSets)

    if not noFiles:
        with open('concat.fa', 'w') as f:
            f.write(output['fastaString'])

        with open('partitions.txt', 'w') as f:
            f.write(output['partFileRaxml'])

    return output


if __name__ == "__main__":
    import bitk3
    parser = argparse.ArgumentParser(
        prog='concatMultipleFasta',
        usage='%(prog)s fastaFile1.fa fastaFile2.fa fasta3File.fa... ',
        description='Makes concatenated alignment of multiple \
        fasta files. Trick here is that they must be in the same order'
    )
    parser.add_argument(
        'fastaFileHandles',
        metavar='file1.fa ...',
        type=argparse.FileType('r'),
        nargs='+',
        help='1 or more fasta files'
    )
    parser.add_argument(
        '-nf',
        '--noFiles',
        action='store_true',
        default=False,
        help='Do not output file'
    )

    args = parser.parse_args()

    fastaFileHandles = args.fastaFileHandles
    print("Let's concatenate the following files:\n{}".format(
        '\n'.join([i.name for i in fastaFileHandles])))

    main(fastaFileHandles, noFiles=args.noFiles)
    print('All done.')
else:
    from bitk3 import bitk3
