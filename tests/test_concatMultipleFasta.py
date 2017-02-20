# -*- coding: utf-8 -*-

"""
test_buildCheABRmakeFasta.py
----------------------------------

Tests for `test_buildCheABRmakeFasta` script.
"""
import pytest
import json
import os
import io

from bitk3 import bitk3
from bitk3 import concatMultipleFasta

myPath = os.getcwd()
dataPath = myPath + '/sampledata/'


class TestNotUsingMiST22:
    def test__parseData(self):
        fasta1 = '>Org1|locus1|Acce1|B|C|D\n \
                  AAAAAAAAAAAAAKKKKKKKKKKKKKKKA\n \
                  >Org2|locus2|Acce2|B|C|D\n \
                  AAAAAAAAAAAAAKKKKKKKKKKKKKKKK\n \
                  >Org3|locus3|Acce3|A|C|D\n \
                  AAAAAAAAAAAAAKKKKKKKKKKKKKKKK\n \
                  >Org4|locus4|Acce4|A|C|E\n \
                  AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
        fasta2 = '>Org5|locus5|Acce5|B|C|D\n \
                  RRRAAAAAAAAAAAKKKKKKKKKKKKKKKA\n \
                  >Org6|locus6|Acce6|B|C|D\n \
                  RRAAAAAAAAAAAAKKKKKKKKKKKKKKKK\n \
                  >Org7|locus7|Acce7|A|C|D\n \
                  RAARAAAAAAAEAAKKKKKKKKKKKKKKKK\n \
                  >Org8|locus8|Acce8|A|C|E\n \
                  RAAARAAAAAAAAAKKKKKKKKKKKKKKKK'

        f1 = io.StringIO(fasta1)
        f2 = io.StringIO(fasta2)
        fs = [f1, f2]

        output = concatMultipleFasta.main(fs, noFiles=True)

        fastaConcat = '>Org1|locus1|Acce1|B|C|D\nAAAAAAAAAAAAAKKKKKKKKKKKKKKKARRRAAAAAAAAAAAKKKKKKKKKKKKKKKA\n>Org2|locus2|Acce2|B|C|D\nAAAAAAAAAAAAAKKKKKKKKKKKKKKKKRRAAAAAAAAAAAAKKKKKKKKKKKKKKKK\n>Org3|locus3|Acce3|A|C|D\nAAAAAAAAAAAAAKKKKKKKKKKKKKKKKRAARAAAAAAAEAAKKKKKKKKKKKKKKKK\n>Org4|locus4|Acce4|A|C|E\nAAAAAAAAAAAAAKKKKKKKKKKKKKKKKRAAARAAAAAAAAAKKKKKKKKKKKKKKKK\n'
        assert fastaConcat == output['fastaString']

        assocTable = [
            [
                "Org1|locus1|Acce1|B|C|D",
                "Org5|locus5|Acce5|B|C|D"
            ],
            [
                "Org2|locus2|Acce2|B|C|D",
                "Org6|locus6|Acce6|B|C|D"
            ],
            [
                "Org3|locus3|Acce3|A|C|D",
                "Org7|locus7|Acce7|A|C|D"
            ],
            [
                "Org4|locus4|Acce4|A|C|E",
                "Org8|locus8|Acce8|A|C|E"
            ]
        ]
        assert assocTable == output['assocTable']

        partFileRaxml = "AUTO, part1 = 1-29\nAUTO, part2 = 30-59\n" 
        assert partFileRaxml == output['partFileRaxml']
        return 0


    def test__parseData_withDiffNumOfSeqMustFail(self):
        fasta1 = '>Org1|locus1|Acce1|B|C|D\n \
                  AAAAAAAAAAAAAKKKKKKKKKKKKKKKA\n \
                  >Org2|locus2|Acce2|B|C|D\n \
                  AAAAAAAAAAAAAKKKKKKKKKKKKKKKK\n \
                  >Org3|locus3|Acce3|A|C|D\n \
                  AAAAAAAAAAAAAKKKKKKKKKKKKKKKK\n \
                  >Org4|locus4|Acce4|A|C|E\n \
                  AAAAAAAAAAAAAKKKKKKKKKKKKKKKK'
        fasta2 = '>Org5|locus5|Acce5|B|C|D\n \
                  RRRAAAAAAAAAAAKKKKKKKKKKKKKKKA\n \
                  >Org6|locus6|Acce6|B|C|D\n \
                  RRAAAAAAAAAAAAKKKKKKKKKKKKKKKK\n \
                  >Org7|locus7|Acce7|B|C|D\n \
                  RRAAAAAAAAAAAAKKKKKKKKKKKKKKKK\n \
                  >Org8|locus8|Acce8|A|C|D\n \
                  RAARAAAAAAAEAAKKKKKKKKKKKKKKKK\n \
                  >Org9|locus9|Acce9|A|C|E\n \
                  RAAARAAAAAAAAAKKKKKKKKKKKKKKKK'

        f1 = io.StringIO(fasta1)
        f2 = io.StringIO(fasta2)
        fs = [f1, f2]

        with pytest.raises(Exception) as exinfo:
            concatMultipleFasta.main(fs)
        assert 'The number of sequences in one of the files is different than the others' in str(exinfo.value)

        return 0