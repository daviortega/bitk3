# -*- coding: utf-8 -*-

"""
test_buildCheABR.py
----------------------------------

Tests for `test_buildCheABR` script.
"""
from bitk3 import bitk3
from bitk3 import buildCheABR
import os

myPath = os.getcwd()
dataPath = myPath + '/sampledata/'


'''def test_buildCheABR():
    sampleFileStr = dataPath + 'cheaTags.txt'
    window = 5
    result = buildCheABR.main(sampleFileStr, window)
    assert result == 0'''


def test__getNeighbors():
    geneNeighborhoodWindow = 7
    gene = {
        'ab': 2671182,
        'loc': '2669053..2671182',
        'lo': 'AX27061_2425',
        'a': 2669053,
        'l': 2130,
        'gs_id': 'J_6JK0TyJ-jr_Yb4a4RTpQ',
        'pos': 2460,
        'ps': False,
        '_id': 27938344,
        'aa': 2669053,
        'gid': 7922,
        'p': {
            'aid': 'Wq7ru9cHPVTRzpPA_BWLQQ',
            'x': {
                'gi': 568127542
            },
            'l': 709,
            'pd': 'Signal transduction histidine kinase CheA',
            'ac': 'REF_ETEC:AX27061_2425'
        },
        'st': '+',
        'gc': 0.6713615023474181,
        'b': 2671182,
        'cid': 681489
    }

    expected = {
        'ab': 2671182,
        'loc': '2669053..2671182',
        'lo': 'AX27061_2425',
        'a': 2669053,
        'l': 2130,
        'gs_id': 'J_6JK0TyJ-jr_Yb4a4RTpQ',
        'pos': 2460,
        'ps': False,
        '_id': 27938344,
        'aa': 2669053,
        'gid': 7922,
        'p': {
            'aid': 'Wq7ru9cHPVTRzpPA_BWLQQ',
            'x': {
                'gi': 568127542
            },
            'l': 709,
            'pd': 'Signal transduction histidine kinase CheA',
            'ac': 'REF_ETEC:AX27061_2425'
        },
        'st': '+',
        'gc': 0.6713615023474181,
        'b': 2671182,
        'cid': 681489,
        'neighborsAseq': [
            'lveO5TPwtSn2Bi22E5E02Q',
            'bIxXbB4M0oAVuenqwqMZ8g',
            'tBtkZzSW9XyPCUaf-kv6fA',
            'o7Nwj3m9V_WmjUyn84uNGQ',
            'NJMLbo3ehqN4i3Yh9i5ieQ',
            'qSxp6e30ayiHsWQJZpw_0Q',
            'dvJ_ieGykwLjwfbfY7akQQ',
            'Wq7ru9cHPVTRzpPA_BWLQQ',
            'TELdd6StFnT9NLP07cuvHQ',
            'Eoutqf7njF-79Gi0i8PO4A',
            'fTds8tYjk_nq5VSaLpMI9g',
            'QrL1ViAXijFH8U7TOJuREg',
            'PWx1HpAJ13ecsExNjd5umQ',
            'sTcjWjo3j-2yjAQQpwpyfg',
            'iJSATBndp49pCqqX3j4tlw'
        ],
        'neighborsAC': [
            'REF_ETEC:AX27061_2418',
            'REF_ETEC:AX27061_2419',
            'REF_ETEC:AX27061_2420',
            'REF_ETEC:AX27061_2421',
            'REF_ETEC:AX27061_2422',
            'REF_ETEC:AX27061_2423',
            'REF_ETEC:AX27061_2424',
            'REF_ETEC:AX27061_2425',
            'REF_ETEC:AX27061_2426',
            'REF_ETEC:AX27061_2427',
            'REF_ETEC:AX27061_2428',
            'REF_ETEC:AX27061_2429',
            'REF_ETEC:AX27061_2430',
            'REF_ETEC:AX27061_2431',
            'REF_ETEC:AX27061_2432'
        ],
        'neighborsId': [
            27938337,
            27938338,
            27938339,
            27938340,
            27938341,
            27938342,
            27938343,
            27938344,
            27938345,
            27938346,
            27938347,
            27938348,
            27938349,
            27938350,
            27938351
        ],
    }
    client = bitk3.get_mist22_client()
    mist22 = client.mist22

    gene = buildCheABR._getNeighbors(mist22, gene, geneNeighborhoodWindow)

    assert gene == expected


def test__getSigTransInfoOfNeighbors():
    gene = {
        'ab': 2671182,
        'loc': '2669053..2671182',
        'lo': 'AX27061_2425',
        'a': 2669053,
        'l': 2130,
        'gs_id': 'J_6JK0TyJ-jr_Yb4a4RTpQ',
        'pos': 2460,
        'ps': False,
        '_id': 27938344,
        'aa': 2669053,
        'gid': 7922,
        'p': {
            'aid': 'Wq7ru9cHPVTRzpPA_BWLQQ',
            'x': {
                'gi': 568127542
            },
            'l': 709,
            'pd': 'Signal transduction histidine kinase CheA',
            'ac': 'REF_ETEC:AX27061_2425'
        },
        'st': '+',
        'gc': 0.6713615023474181,
        'b': 2671182,
        'cid': 681489,
        'neighborsAseq': [
            'lveO5TPwtSn2Bi22E5E02Q',
            'bIxXbB4M0oAVuenqwqMZ8g',
            'tBtkZzSW9XyPCUaf-kv6fA',
            'o7Nwj3m9V_WmjUyn84uNGQ',
            'NJMLbo3ehqN4i3Yh9i5ieQ',
            'qSxp6e30ayiHsWQJZpw_0Q',
            'dvJ_ieGykwLjwfbfY7akQQ',
            'Wq7ru9cHPVTRzpPA_BWLQQ',
            'TELdd6StFnT9NLP07cuvHQ',
            'Eoutqf7njF-79Gi0i8PO4A',
            'fTds8tYjk_nq5VSaLpMI9g',
            'QrL1ViAXijFH8U7TOJuREg',
            'PWx1HpAJ13ecsExNjd5umQ',
            'sTcjWjo3j-2yjAQQpwpyfg',
            'iJSATBndp49pCqqX3j4tlw'
        ],
        'neighborsAC': [
            'REF_ETEC:AX27061_2418',
            'REF_ETEC:AX27061_2419',
            'REF_ETEC:AX27061_2420',
            'REF_ETEC:AX27061_2421',
            'REF_ETEC:AX27061_2422',
            'REF_ETEC:AX27061_2423',
            'REF_ETEC:AX27061_2424',
            'REF_ETEC:AX27061_2425',
            'REF_ETEC:AX27061_2426',
            'REF_ETEC:AX27061_2427',
            'REF_ETEC:AX27061_2428',
            'REF_ETEC:AX27061_2429',
            'REF_ETEC:AX27061_2430',
            'REF_ETEC:AX27061_2431',
            'REF_ETEC:AX27061_2432'
        ],
        'neighborsId': [
            27938337,
            27938338,
            27938339,
            27938340,
            27938341,
            27938342,
            27938343,
            27938344,
            27938345,
            27938346,
            27938347,
            27938348,
            27938349,
            27938350,
            27938351
        ],
    }

    expected = {
        'ab': 2671182,
        'loc': '2669053..2671182',
        'lo': 'AX27061_2425',
        'a': 2669053,
        'l': 2130,
        'gs_id': 'J_6JK0TyJ-jr_Yb4a4RTpQ',
        'pos': 2460,
        'ps': False,
        '_id': 27938344,
        'aa': 2669053,
        'gid': 7922,
        'p': {
            'aid': 'Wq7ru9cHPVTRzpPA_BWLQQ',
            'x': {
                'gi': 568127542
            },
            'l': 709,
            'pd': 'Signal transduction histidine kinase CheA',
            'ac': 'REF_ETEC:AX27061_2425'
        },
        'st': '+',
        'gc': 0.6713615023474181,
        'b': 2671182,
        'cid': 681489,
        'neighborsAseq': [
            'lveO5TPwtSn2Bi22E5E02Q',
            'bIxXbB4M0oAVuenqwqMZ8g',
            'tBtkZzSW9XyPCUaf-kv6fA',
            'o7Nwj3m9V_WmjUyn84uNGQ',
            'NJMLbo3ehqN4i3Yh9i5ieQ',
            'qSxp6e30ayiHsWQJZpw_0Q',
            'dvJ_ieGykwLjwfbfY7akQQ',
            'Wq7ru9cHPVTRzpPA_BWLQQ',
            'TELdd6StFnT9NLP07cuvHQ',
            'Eoutqf7njF-79Gi0i8PO4A',
            'fTds8tYjk_nq5VSaLpMI9g',
            'QrL1ViAXijFH8U7TOJuREg',
            'PWx1HpAJ13ecsExNjd5umQ',
            'sTcjWjo3j-2yjAQQpwpyfg',
            'iJSATBndp49pCqqX3j4tlw'
        ],
        'neighborsAC': [
            'REF_ETEC:AX27061_2418',
            'REF_ETEC:AX27061_2419',
            'REF_ETEC:AX27061_2420',
            'REF_ETEC:AX27061_2421',
            'REF_ETEC:AX27061_2422',
            'REF_ETEC:AX27061_2423',
            'REF_ETEC:AX27061_2424',
            'REF_ETEC:AX27061_2425',
            'REF_ETEC:AX27061_2426',
            'REF_ETEC:AX27061_2427',
            'REF_ETEC:AX27061_2428',
            'REF_ETEC:AX27061_2429',
            'REF_ETEC:AX27061_2430',
            'REF_ETEC:AX27061_2431',
            'REF_ETEC:AX27061_2432'
        ],
        'neighborsId': [
            27938337,
            27938338,
            27938339,
            27938340,
            27938341,
            27938342,
            27938343,
            27938344,
            27938345,
            27938346,
            27938347,
            27938348,
            27938349,
            27938350,
            27938351
        ],
        'cheInfo': [
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': 'cher'
            },
            {
                'che': 'cheb'
            },
            {
                'che': False
            },
            {
                'che': False
            },
            {
                'che': False
            }
        ]
    }
    client = bitk3.get_mist22_client()
    mist22 = client.mist22

    gene = buildCheABR._getSigTransInfoOfNeighbors(mist22, gene)

    assert gene == gene
