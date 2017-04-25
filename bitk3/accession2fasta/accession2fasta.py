from bitk3 import bitk3
import warnings


class Accession2fasta:
    def __init__(self, accessionList):
        self.fasta = self._make_fasta(accessionList)

    def _throw_warning(self, accession):
        warnings.warn('Invalid accession - skipping - {}'.format(
            accession), Warning)
        return 1

    def _make_fasta(self, accessionList):
        accessions = []
        for line in accessionList:
            accession = line.replace('\n', '').replace(' ', '')
            if bitk3.isValidRefSeqAccession(accession):
                accessions.append(accession)
            else:
                self._throw_warning(accession)

        genes, badTypes = bitk3.accessionToMist22GeneInfo(accessions)
        for badType in badTypes:
            self._throw_warning(badType)
        genes = bitk3.addSeqToMist22GeneInfo(genes)
        genes = bitk3.addBitk3tagToMist22GeneInfo(genes)

        fasta = ''
        for i, gene in enumerate(genes):
            if gene:
                fasta += '>{}\n{}\n'.format(gene['bitk3tag'], gene['seq'])
            else:
                self._throw_warning(accessions[i])
        return fasta



