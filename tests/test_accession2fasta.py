# -*- coding: utf-8 -*-

"""
test_accession2fasta.py
----------------------------------

Tests for `accession2fasta` script.
"""
import bitk3.accession2fasta.accession2fasta as ac2fa
import os
import pytest
import warnings

myPath = os.getcwd()
dataPath = myPath + '/sampledata/'


@pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason="Skipping this test on Travis CI."
)
class TestUsingMiST22:
    def test_accession2fasta(self):
        sampleFile = dataPath + 'accessionList.txt'
        expectedFile = dataPath + 'accession2fasta.fa'

        with warnings.catch_warnings(record=True) as w:
            with open(sampleFile, 'r') as f:
                ac = ac2fa.Accession2fasta(f)
            assert len(w) == 2
            assert issubclass(w[-1].category, Warning)
            assert 'Invalid accession' in str(w[-1].message)

        expectedFasta = ''
        with open(expectedFile, 'r') as f:
            for line in f:
                expectedFasta += line
        assert ac.fasta == expectedFasta

        return 1
