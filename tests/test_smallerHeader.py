# -*- coding: utf-8 -*-

"""
test_smallerHeader.py
----------------------------------

Tests for `smallerHeader` script.
"""
import bitk3.smallerHeader.smallerHeader as sH
import os


myPath = os.getcwd()
dataPath = myPath + '/sampledata/'


sampleFile = dataPath + 'fasta.with.bitk.tags.fa'
expectedFile = dataPath + 'smallerHeader.fa'
result = sH.SmallerHeader(sampleFile)


class TestSmallHeader:

    def test_transdic(self):
        expected = [
            {'s': 'XX0XX', 'l': 'Org1|locus1|Acce1|B|C|D'},
            {'s': 'XX1XX', 'l': 'Org2|locus2|Acce2|B|C|D'},
            {'s': 'XX2XX', 'l': 'Org3|locus3|Acce3|A|C|D'},
            {'s': 'XX3XX', 'l': 'Org4|locus4|Acce4|A|C|E'}
        ]
        assert result.transdic == expected
        return 1

    def test_makingNewFastaDictionary(self):
        expected = {
            'XX0XX': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKA',
            'XX1XX': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
            'XX2XX': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
            'XX3XX': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
        }

        assert result.newFastaDic == expected
        return 1

    def test_makingFastaFormatedString(self):
        expectedFasta = ''
        with open(expectedFile, 'r') as f:
            for line in f:
                expectedFasta += line

        newFasta = str(result)
        assert newFasta == expectedFasta
        return 1
