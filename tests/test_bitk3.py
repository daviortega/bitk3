# -*- coding: utf-8 -*-

"""
test_bitk3
----------------------------------

Tests for `bitk3` module.
"""

import pytest
import json

from bitk3 import bitk3
import os

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
        'Org1|locus1|Acce1|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKA',
        'Org2|locus2|Acce2|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org3|locus3|Acce3|A|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org4|locus4|Acce4|A|C|E': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
    }
    return None


def test_fastaReader_fileMissingHeaderMarker():
    """Test for invalid fasta file - missing '>' marker"""
    sampleFile = dataPath + 'invalid_fasta.fa'
    line = 'Org1|locus1|Acce1|B|C|D'
    with pytest.raises(Exception) as exinfo:
        bitk3.fastaReader(sampleFile)
    assert 'Invalid FASTA file. Header line must begin with a greater than \
symbol\nLine: ' + line + '\n\n' in str(exinfo.value)


def test_fasta_class():
    sampleFile = dataPath + 'fasta.with.bitk.tags.fa'
    myDict = {}
    for seqObject in bitk3.Fasta(sampleFile).stream():
        myDict[seqObject['h']] = seqObject['s']
    assert myDict == {
        'Org1|locus1|Acce1|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKA',
        'Org2|locus2|Acce2|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org3|locus3|Acce3|A|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org4|locus4|Acce4|A|C|E': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
    }
    return None


def test_fasta_class_makeDictionary():
    sampleFile = dataPath + 'fasta.with.bitk.tags.fa'
    myDict = bitk3.Fasta(sampleFile).makeDictionary()
    assert myDict == {
        'Org1|locus1|Acce1|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKA',
        'Org2|locus2|Acce2|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org3|locus3|Acce3|A|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org4|locus4|Acce4|A|C|E': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
    }
    return None


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


def test_countFeatureInHeaders_alternativeTagSeparator():
    """Tests if it can count"""
    sampleFile = dataPath + 'fasta.with.bitk.alttags.fa'
    seqInfo = bitk3.fastaReader(sampleFile)
    countsDict3 = bitk3.countFeaturesInHeaders(seqInfo[0], 3, sep='-')
    assert countsDict3 == {
        'B': 2,
        'A': 2
    }
    return None


def test_bitk3tagToAccession():
    """Test if it can extract accession from tag"""
    sampleTags = [
        'Ac_xyl_7922|AX27061_2425|REF_ETEC:AX27061_2425',
        'Al_die_2140|B5T_03398|YP_006822008.1',
    ]
    accession = [
        'REF_ETEC:AX27061_2425',
        'YP_006822008.1'
    ]

    for i, tag in enumerate(sampleTags):
        Xted = bitk3.bitk3tagToAccession(tag)
        assert accession[i] == Xted

    return None


def test_getAseqFromMist22Gene():
    """ Test if it can get Aseq from gene info in MiST22 and pass None if there \
    is no protein info """
    sampleFile = dataPath + 'mistGenes.json'
    expected = ['ulF-SXsxtnFn7TYkYb3hnw', None]
    with open(sampleFile, 'r') as f:
        genes = json.load(f)
    results = []
    for i, gene in enumerate(genes):
        result = bitk3.getAseqFromMist22Gene(gene)
        results.append(result)
    assert expected == results
    return 0


def test_getMistIDFromMist22Gene():
    """ Test if it can get internal protein gene ID from gene \
    info in MiST22 """
    sampleFile = dataPath + 'mistGenes.json'
    expected = [333724, 333725]
    with open(sampleFile, 'r') as f:
        genes = json.load(f)
    results = []
    for i, gene in enumerate(genes):
        result = bitk3.getMistIDFromMist22Gene(gene)
        results.append(result)
    assert expected == results
    return 0


def test_getAccessionFromMist22Gene():
    """ Test if it can get Accession from gene info in MiST22 \
    and pass None if there is no protein info"""
    sampleFile = dataPath + 'mistGenes.json'
    expected = ['YP_003444191.1', None]
    with open(sampleFile, 'r') as f:
        genes = json.load(f)
    results = []
    for i, gene in enumerate(genes):
        result = bitk3.getAccessionFromMist22Gene(gene)
        results.append(result)
    assert expected == results
    return 0


def test_getLocusFromMist22Gene():
    """ Test if it can get Locus from gene info in MiST22 """
    sampleFile = dataPath + 'mistGenes.json'
    expected = ['Alvin_2240', 'Alvin_2241']
    with open(sampleFile, 'r') as f:
        genes = json.load(f)
    results = []
    for i, gene in enumerate(genes):
        result = bitk3.getLocusFromMist22Gene(gene)
        results.append(result)
    assert expected == results
    return 0


def test_getGenomeIDFromMist22Gene():
    """ Test if it can get GenomeID from gene info in MiST22 """
    sampleFile = dataPath + 'mistGenes.json'
    expected = [90, 90]
    with open(sampleFile, 'r') as f:
        genes = json.load(f)
    results = []
    for i, gene in enumerate(genes):
        result = bitk3.getGenomeIDFromMist22Gene(gene)
        results.append(result)
    assert expected == results
    return 0


def test_isValidRefSeqAccession():
    """ Test is input is a valid Accession """
    fixtures = [
        [True, 'REF_ETEC:AX27061_2429'],
        [True, 'YP_001452647.1'],
        [True, 'YP_574071.1'],
        [False, None]
    ]

    for fixture in fixtures:
        expected, inp = fixture
        assert expected == bitk3.isValidRefSeqAccession(inp)

    return 0


def test_getSeqFromAseq():
    aseqs = [
        'laERqmSRb11rhvmY2Did7Q',
        '33EXyIiXMk5ugL1mWicTqw'
    ]

    seqs = [
        'MSINMAEFHQVFFEESHEHLENMEQLLIAINLQSPDPEELNTIFRAAHSIKG\
GSGIFGFTALSSVTHVMENLLDKVRKGTFELSSGIIDLLLKTVDTLSHILSLY\
REEEPIDWQQVEFAKNQLVAALNGDPFLTTSPEAIPATPTKTEITQTAVVVSD\
ASHQDDIGFGFFEDDVELTLAIENQDFGFFDEAFTAKEISVEDINKTLLDSQL\
VTNTEPDDDLGFGFFESLTPESVDSEINALVHIQANKAIAAIEPLSKPVNSPS\
VPKIPRYTGNTAESTKVAPADPVDPSPVLNKSVASGTKSTPTAATKKGNTSTQ\
DATLRVETSKIDTLVNLAGELVITQSMLTLIGNEISGELGERLKTALVELERN\
TREMQEAVMSVRMLPVSFVFNRFHRLVRDLSDQLGKNVNLVIEGGNTEIDKGM\
IEKLVDPLTHLVRNSLDHGIEKPEVRRLLGKAEIAQLSLRASQRGGNIVIAVH\
DDGAGLHREKILQKARENNMSVTDNMPDKQIWQLIFAAGFSTAKEITDVSGRG\
VGMDVVRRNIEALGGRIDIDSVAGQGATFEIQLPLTLAIVDGMSVSVGKQIYI\
LPLVHIIESIQPHTEQLKYLAQERLIRVREEYLPLLNLHQLMEITPFAKCPEE\
GIVVLLESNNKRFGLCVDALVGQQQVVIKSLEKHYRRIPGVSGATIMGDGSVA\
LILDVESLAQQIKN',
        'MGILDIFGNTEERQHQDEERQRYFQLLDNSRNNFMIADSDRNIIYANKAVLNMLSEAEAD\
IRKQLPQFSVARVIGSNIDIFHVNPAHQRNMLERLTQSHTAQISIGKRIFKLILTPIISRD\
NKHLGTGVEWIDRTESIDAERATQRILEALNNTSTNVMIADANRTIIYMNRSVEAMLRRSE\
SEIRQVLPHFSVDKILGSSMDIFHRNPAHQASLLDKLDRKYESQIQVASCHFRLTASPIIL\
TSGERLGSVVEWLDRTEEVQVEQEIARIVAAAAAGDFSQRVDSHGKQGFFLMLANSLNSLI\
ETSDRGLQDVARVLMAMAEGDLTTRIYNEYEGTFNDLKNYSNQTAEKLSYMIRDIQKAADT\
INTASSEIAQGNADLSSRTEEQASSLEETSASMEELTGTVKLNADNASQANALASKAADVA\
EDGGELIQQVVQTMASINESARKIADIIGVIDGIAFQTNILALNAAVEAARAGEQGRGFAV\
VASEVRSLAQRSANAAKDIKALISDSVSKIDSGNNLVGKSGDTMKEIVIAIKRVNDIMAEI\
ASASNEQAIGIDEIGKAVVQMDEMTQQNAALVEEAAAAAESMQSQAQQLADSVANFKVDEE\
TRSARHTTELKKIPQKIPTLARVTPKPKAMTPKLNKADQDEWEEF'
    ]

    aseq2seq = bitk3.getSeqFromAseq(aseqs)

    for i, aseq in enumerate(aseqs):
        assert aseq2seq[aseq] == seqs[i]

    return 0


def test_isMSA_False():
    seqInfo = {
        'Org1|locus1|Acce1|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKA',
        'Org2|locus2|Acce2|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKKA',
        'Org3|locus3|Acce3|A|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org4|locus4|Acce4|A|C|E': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
    }

    assert not bitk3.isMSA(seqInfo)

    return 0


def test_isMSA_True():
    seqInfo = {
        'Org1|locus1|Acce1|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKA',
        'Org2|locus2|Acce2|B|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org3|locus3|Acce3|A|C|D': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK',
        'Org4|locus4|Acce4|A|C|E': 'AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
    }

    assert bitk3.isMSA(seqInfo)

    return 0

@pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason="Skipping this test on Travis CI."
)
class TestUsingMist22:
    def test_accession2GeneInfo(self):
        accessions = [
            'REF_ETEC:AX27061_2429',
            'YP_001452647.1',
            'YP_574071.1',
            None
        ]

        genes, badTypes = bitk3.accessionToGeneInfo(accessions)
        assert isinstance(genes, list)
        assert isinstance(badTypes, list)

        listOfRetrievedAC = [item['p']['ac'] for item in genes if item]
        listOfRetrievedAC += badTypes

        for ac in accessions:
            assert ac in listOfRetrievedAC

        return 0

    def test_addBitk3tagTomist22GeneInfo(self):
        """ Test if it can generate a bitk3 tag from gene info in MiST22 """
        sampleFile = dataPath + 'mistGenes.json'
        expected = [
            'Al_vin_90|Alvin_2240|YP_003444191.1',
            'Al_vin_90|Alvin_2241|None'
        ]
        with open(sampleFile, 'r') as f:
            genes = json.load(f)
        genes = bitk3.addBitk3tagTomist22GeneInfo(genes)

        listOfBitk3tag = [gene['bitk3tag'] for gene in genes]

        assert set(expected) == set(listOfBitk3tag)
        return 0
