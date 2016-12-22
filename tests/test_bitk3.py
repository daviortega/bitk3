# -*- coding: utf-8 -*-

"""
test_bitk3
----------------------------------

Tests for `bitk3` module.
"""

import pytest

from contextlib import contextmanager
from click.testing import CliRunner

from bitk3 import bitk3
from bitk3 import cli
import os

def test_command_line_interface():
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert 'bitk3.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help  Show this message and exit.' in help_result.output

myPath = os.getcwd()
dataPath = myPath + '/sampledata/'

def test_fastaReader_lists():
    """Simple test for fasta_reader to see if it spits the correct list"""
    sampleFile = dataPath + 'fasta.with.bitk.tags.fa'
    seqInfo = bitk3.fastaReader(sampleFile)
    assert seqInfo[1] == [
        'Org1|locus1|Acce1|B|C|D',
        'Org2|locus2|Acce2|B|C|D',
        'Org3|locus3|Acce3|A|C|D',
        'Org4|locus4|Acce4|A|C|E'
    ]
    assert seqInfo[1] != [
        'Org2|locus2|Acce2|B|C|D',
        'Org3|locus3|Acce3|A|C|D',
        'Org4|locus4|Acce4|A|C|E',
        'Org1|locus1|Acce1|B|C|D',
    ]
    return None

def test_fastaReader_dictionary():
    """Simple test for fasta_reader to see if spits the correct dictionary"""
    sampleFile = dataPath + 'fasta.with.bitk.tags.fa'
    seqInfo = bitk3.fastaReader(sampleFile)
    assert seqInfo[0] == {
        'Org1|locus1|Acce1|B|C|D':'AAAAAAAAAAAAAKKKKKKKKKKKKKKKA',
        'Org2|locus2|Acce2|B|C|D':'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org3|locus3|Acce3|A|C|D':'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org4|locus4|Acce4|A|C|E':'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
    }
    return None

def test_fastaReader_fileMissingHeaderMarker():
    """Test for invalid fasta file - missing '>' marker"""
    sampleFile = dataPath + 'invalid_fasta.fa'
    line = 'Org1|locus1|Acce1|B|C|D'
    with pytest.raises(Exception) as exinfo:
        SeqInfo = bitk3.fastaReader(sampleFile)
    assert 'Invalid FASTA file. Header line must begin with a greater than symbol\nLine: ' + line + '\n\n' in str(exinfo.value)

def test_countFeatureInHeaders_counting():
    """Tests if it can count"""
    sampleFile = dataPath + 'fasta.with.bitk.tags.fa'
    seqInfo = bitk3.fastaReader(sampleFile)
    countsDict3 = bitk3.countFeaturesInHeaders(seqInfo[0], 3)
    assert countsDict3 == {
        'B': 2,
        'A': 2
    }
    return None