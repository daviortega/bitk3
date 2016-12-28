# -*- coding: utf-8 -*-
"""This is the BITK for python 3"""

#List of contants
BITKTAGSEP = '|' #TAG separator
BITKGENSEP = '_' #GENome separator
COGOUTPUTSEP = '\t' #Separator for cog output
MAXSIZEASEQ2SEQ = 100000

class InputError(Exception):
    """ Exception raised for errors of the input
    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        super(InputError, self).__init__(message)
        self.message = message

class Fasta:
    """ Ogun Adebali's approach to stream fasta files"""
    def __init__(self, input):
        self.file = input

    def stream(self, bufsize=4096):
        def chunk2seqDict(chunk):
            lines = chunk.split('\n')
            header = lines[0]
            del lines[0]
            sequence = ''.join(lines)
            seqObject = {'h': header, 's': sequence}
            return seqObject

        filein = open(self.file, 'r')
        delimiter = '\n>'
        buf = ''
        justStarted = True
        while True:
            newbuf = filein.read(bufsize)
            if not newbuf:
                yield chunk2seqDict(buf)
                return
            buf += newbuf
            sequenceChunks = buf.split(delimiter)
            for chunk in sequenceChunks[0:-1]:
                if justStarted and chunk.startswith('>'):
                    chunk = chunk[1:]
                    justStarted = False
                yield chunk2seqDict(chunk)
            buf = sequenceChunks[-1]

    def makeDictionary(self):
        """ Class method to make dictionary out of streaming """
        myDict = {}
        for seqObject in self.stream():
            myDict[seqObject['h']] = seqObject['s']
        return myDict


def fastaReader(datafile):
    """Reads fasta files into a dictionary where keys are the sequence headers
    and values are sequences and a list of the headers to preserved order"""
    listOrder = []
    seqDic = {}
    fastaBuffer = None

    fileHandle = open(datafile, 'w')
    line = fastaBuffer if fastaBuffer else fileHandle.readline()
    while line:
        if line[0] != '>':
            raise Exception('Invalid FASTA file. Header line must begin with a greater than symbol\nLine: ' + line + '\n\n')
        name = line[1:-1]
        listOrder.append(name)
        seqDic[name] = ''
        line = fileHandle.readline()
        while line:
            if line[0] != '>':
                seqDic[name] += line
                line = fileHandle.readline()
                continue
            fastaBuffer = line
            break
        seqDic[name] = seqDic[name].replace('\n', '').replace(' ', '')
    return seqDic, listOrder

def countFeaturesInHeaders(listOfHeaders, pos, sep = BITKTAGSEP):
    """ Return all unique info in a particular position of the tag and
        count in how many headers they appear"""
    result = {}
    for tag in listOfHeaders:
        info = tag.split(sep)[pos]
        if info not in result.keys():
            result[info] = 1
        else:
            result[info] += 1
    return result
