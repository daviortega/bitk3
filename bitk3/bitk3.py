# -*- coding: utf-8 -*-
"""This is the BITK for python 3"""
import SeqDepot
import pymongo
import sys
import re


# List of contants
BITKTAGSEP = '|'  # TAG separator
BITKGENSEP = '_'  # GENome separator
COGOUTPUTSEP = '\t'  # Separator for cog output
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


def isMSA(seqDic={}):
    """ test if it is a multiple sequence alignment \
    by checking the sequences for being the same length
    """
    ismsa = True
    first = True
    for header in seqDic.keys():
        if first:
            aligLen = len(seqDic[header])
        elif aligLen != len(seqDic[header]):
            ismsa = False
            print('Not an MSA - offending sequence:{}'.format(header))
            break
        first = False

    return ismsa


def fastaReader(datafile):
    """Reads fasta file into a dictionary where keys are the sequence headers
    and values are sequences and a list of the headers to preserved order"""

    fileHandle = open(datafile, 'r')
    seqDic, listOrder = fastaReaderByHandle(fileHandle)
    fileHandle.close()

    return seqDic, listOrder


def fastaReaderByHandle(fileHandle):
    """Takes a fasta file handle into a dictionary where keys are the
    sequence headers and values are sequences and a list of the headers
    to preserved order"""
    listOrder = []
    seqDic = {}
    fastaBuffer = None

    line = fastaBuffer if fastaBuffer else fileHandle.readline()
    while line:
        if line[0] != '>':
            raise Exception('Invalid FASTA file. Header line must begin with \
a greater than symbol\nLine: ' + line + '\n\n')
        name = line[1:-1]
        listOrder.append(name)
        seqDic[name] = ''
        line = fileHandle.readline()
        while line:
            line = line.replace(' ','')
            if line[0] != '>':
                seqDic[name] += line
                line = fileHandle.readline()
                continue
            fastaBuffer = line
            break
        seqDic[name] = seqDic[name].replace('\n', '').replace(' ', '')
    return seqDic, listOrder


def countFeaturesInHeaders(listOfHeaders, pos, sep=BITKTAGSEP):
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


def countFeaturesInHeaders_streaming(filename, pos, sep=BITKGENSEP):
    """ Return all unique info in a particular position of the tag and
        count in how many headers they appear using streaming"""

    result = {}
    for seqObject in Fasta(filename).stream():
        info = seqObject['h'].split(sep)[pos]
        if info not in result.keys():
            result[info] = 1
        else:
            result[info] += 1
    return result


def bitk3tagToAccession(bitk3tag=''):
    """ Extract accession number from bitk3 tag """
    return bitk3tag.split(BITKTAGSEP)[2]


def get_mist22_client():
    """ Get mist22 client - soon to be deprecated"""
    try:
        client = pymongo.MongoClient('localhost', 27019)
        client.mist22.genes.find_one()
    except TypeError:
        print("You must open a tunnel with ares.bio.utk.edu: ssh -p 32790 \
        -f -N -L 27019:localhost:27017 unsername@ares.bio.utk.edu")
        sys.exit()
    except pymongo.errors.ConnectionFailure:
        print("You must open a tunnel with ares.bio.utk.edu: ssh -p 32790 \
        -f -N -L 27019:localhost:27017 unsername@ares.bio.utk.edu")
        sys.exit()

    return client


def isValidRefSeqAccession(accession=''):
    """ Test if accession is a valid accession """
    pattern = '((NC|AC|NG|NT|NW|NZ|NM|NR|XM|XR|NP|AP|XP|YP|ZP)_[0-9]+)|(REF.*:[A-Z]{1,}.*_[0-9]+)'
    if re.match(pattern, str(accession)):
        match = True
    else:
        match = False
    return match


def cleanListOfAccessions(accessionList=[]):
    newList = []
    badTypes = []
    for ac in accessionList:
        if isValidRefSeqAccession(ac):
            newList.append(ac)
        else:
            badTypes.append(ac)

    return newList, badTypes


def accessionToMist22GeneInfo(accessionList=[]):
    """ Read a list of accession and returns a generator to run \
    over all info from the genes"""

    client = get_mist22_client()
    mist22 = client.mist22

    newList, badTypes = cleanListOfAccessions(accessionList)

    genes = mist22.genes.find(
        {
            'p.ac': {
                '$in': newList
            }
        }
    )

    geneInfo = []
    for gene in genes:
        geneInfo.append(gene)

    client.close
    return geneInfo, badTypes


def getAseqFromMist22Gene(gene={}):
    """ returns the Aseq of a gene information from Mist22 """
    aseq = None
    if 'p' in gene.keys():
        aseq = gene['p']['aid']
    return aseq


def getMistIDFromMist22Gene(gene={}):
    """ returns the internal ID of a gene information from Mist22 """
    return gene['_id']


def getAccessionFromMist22Gene(gene={}):
    """ returns the accession number of a gene information from Mist22 """
    ac = None
    if 'p' in gene.keys():
        ac = gene['p']['ac']
    return ac


def getLocusFromMist22Gene(gene={}):
    """ returns the locus number of a gene information from Mist22 """
    lo = None
    if 'lo' in gene.keys():
        lo = gene['lo']
    return lo


def getGenomeIDFromMist22Gene(gene={}):
    """ returns the internal genome ID of a gene information from Mist22 """
    mGid = gene['gid']
    return mGid


def getSpeciesNameFromMist22Genome(genome={}):
    """ returns the species name of a genome information from Mist22 """
    return genome['sp']


def getGenusNameFromMist22Genome(genome={}):
    """ returns the genus name of a genome information from Mist22 """
    return genome['g']


def getGenomeIDFromMist22Genome(genome={}):
    """ returns the internal genome ID of a genome information from Mist22 """
    return genome['_id']


def getGenomeInfoFromMistIDs(gids=[]):
    """
    Read a list of internal mist22 genome IDs and returns
    a dictionary with the ID as key and the genome info as value

    Argument Keywords:
    gids -- List of internal mist22 genome IDs

    Returns:
    {gids : genomeInfo}

    """
    client = get_mist22_client()
    mist22 = client.mist22
    genomes = mist22.genomes.find(
        {
            '_id': {
                '$in': gids
            }
        }
    )

    genDic = {}
    for gen in genomes:
        genDic[gen['_id']] = gen
    return genDic


def addBitk3tagToMist22GeneInfo(genes=[]):
    """ Build bitk3 tag for fasta sequences from list of mist22 genes """
    gids = []

    oldGenes = []

    for gene in genes:
        print(gene)
        if isinstance(gene, dict):
            mistGenomeId = getGenomeIDFromMist22Gene(gene)
            if mistGenomeId not in gids:
                gids.append(mistGenomeId)
            oldGenes.append(gene)
        else:
            oldGenes.append(None)

    genDic = getGenomeInfoFromMistIDs(gids)

    newGenes = []
    for gene in oldGenes:
        if gene:
            mistGenomeId = getGenomeIDFromMist22Gene(gene)
            lo = getLocusFromMist22Gene(gene)
            accession = getAccessionFromMist22Gene(gene)
            genus = genDic[mistGenomeId]['g']
            species = genDic[mistGenomeId]['sp']

            bitk3genID = (str(genus[:2]) + BITKGENSEP + str(species[:3]) +
                        BITKGENSEP + str(mistGenomeId))

            bitk3tag = (bitk3genID + BITKTAGSEP + str(lo) +
                        BITKTAGSEP + str(accession))

            gene['bitk3tag'] = bitk3tag
        newGenes.append(gene)


    return newGenes


def getSeqFromAseq(aseqs=[]):
    aseq2seq = {}
    sd = SeqDepot.new()
    seqs = sd.find(aseqs, {'fields': 's'})
    if seqs:
        for seq in seqs:
            aseq2seq[seq['data']['id']] = seq['data']['s']
    return aseq2seq


def addSeqToMist22GeneInfo(genes=[]):
    aseqs = []
    for gene in genes:
        aseq = getAseqFromMist22Gene(gene)
        if aseq:
            aseqs.append(aseq)
    aseq2seq = getSeqFromAseq(aseqs)

    newGenes = []
    for gene in genes:
        aseq = getAseqFromMist22Gene(gene)
        gene['seq'] = aseq2seq[aseq]
        newGenes.append(gene)

    return newGenes
