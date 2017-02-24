if __package__ == '':
    from os import sys, path
    sys.path.append(path.dirname(path.dirname(__file__)))
    import bitk3
else:
    from bitk3 import bitk3
import warnings
import math


class SmallerHeader:
    def __init__(self, filename):
        self.fasta, self.tags = bitk3.fastaReader(filename)
        self.transdic, self.newFastaDic, self.newTags = self._translate()

    def _numZeros(self):
        n = len(self.tags)
        return math.floor(math.log10(n))

    def _translate(self):
        numZeros = self._numZeros()
        newTags = ['XX{:0{numZ}d}XX'.format(
            i,
            numZ=numZeros
        ) for i in range(len(self.tags))]
        transdic = {}
        newFastaDic = {}

        for i, tag in enumerate(self.tags):
            transdic[newTags[i]] = tag
            newFastaDic[newTags[i]] = self.fasta[tag]

        return transdic, newFastaDic, newTags

    def __str__(self):
        fasta = ''
        for tag in self.newTags:
            fasta += '>{}\n{}\n'.format(tag, self.newFastaDic[tag])
        return fasta
