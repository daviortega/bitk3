# -*- coding: utf-8 -*-

"""
test_buildCheABRmakeFasta.py
----------------------------------

Tests for `test_buildCheABRmakeFasta` script.
"""
import pytest
import json
import os

from bitk3 import bitk3
from bitk3 import buildCheABRmakeFasta

myPath = os.getcwd()
dataPath = myPath + '/sampledata/'


class TestNotUsingMiST22:
    def test_parseCheInfo(self):
        with open(dataPath + 'mistGenes.more.json', 'r') as f:
            genes = json.load(f)
        with open(dataPath + 'seqInfo.sample.new.json', 'r') as f:
            expected = json.load(f)

        seqInfo = buildCheABRmakeFasta._parseCheInfo(genes)

        assert seqInfo == expected


@pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason="Skipping this test on Travis CI."
)
class TestUsingMiST22:
    def test_buildCheABRmakeFasta_execution(self):
        sampleFileStr = dataPath + 'cheaTags.txt'
        window = 7
        buildCheABRmakeFasta.main(
            sampleFileStr,
            window,
            verbose=False,
            noFiles=True
        )

    def test__getNeighbors(self):
        geneNeighborhoodWindow = 7

        fixtures = [
            [
                {
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
                },
                {
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
            ],
            [
                {
                    "_id": 1113403,
                    "cid": 633,
                    "gid": 319,
                    "l": 2163,
                    "a": 1052440,
                    "p": {
                    "ac": "NP_233475.1",
                    "l": 720,
                    "aid": "AtBd-8mQw2yws-ZrpoEQWA",
                    "x": {
                        "gi": 15601844
                    },
                    "pd": "chemotaxis protein CheA",
                    "n": "similar to GB:M34669 SP:P07363 PID:145529 PID:145530 GB:U00096; identified by sequence similarity"
                    },
                    "ab": 1052440,
                    "pos": 1006,
                    "lo": "VCA1095",
                    "ps": False,
                    "aa": 1050278,
                    "st": "-",
                    "b": 1050278,
                    "loc": "complement(1050278..1052440)",
                    "gc": 0.490522422561257,
                    "gs_id": "ohUew9_0Y21HMztN2UkYBw"
                },
                {
                    "b": 1050278,
                    "a": 1052440,
                    "ps": False,
                    "l": 2163,
                    "gid": 319,
                    "gs_id": "ohUew9_0Y21HMztN2UkYBw",
                    "pos": 1006,
                    "aa": 1050278,
                    "_id": 1113403,
                    "lo": "VCA1095",
                    "neighborsAC": [
                        "NP_233469.1",
                        None,
                        "NP_233470.1",
                        "NP_233471.1",
                        "NP_233472.1",
                        "NP_233473.1",
                        "NP_233474.1",
                        "NP_233475.1",
                        "NP_233476.1",
                        "NP_233477.1",
                        "NP_233478.1",
                        "NP_233479.1",
                        "NP_233480.1",
                        "NP_233481.1",
                        "NP_233482.1"
                    ],
                    "st": "-",
                    "neighborsAseq": [
                        "QpctL60-ANyBk1aGp3W-hA",
                        None,
                        "jc8tmp9XtR6lTPMX3oSAlg",
                        "xj70OvOV0Ok4UHkJQriVvA",
                        "ConBlJawh0S9eNcwpUdStA",
                        "aXjcWvPeespvQZp3aap5LQ",
                        "qyHHt7NIZE1oStynXSwaLA",
                        "AtBd-8mQw2yws-ZrpoEQWA",
                        "MBG9AyB9JZ0g2yFduk70_w",
                        "6msMWTG01YKLCBo5dqqL4A",
                        "gA7NzCfTVtXu1kVp8hPZWg",
                        "x0OMwjtP2hGSk9pS-3btXg",
                        "d4_i8687rdGbPEaz-nj7_w",
                        "8DbYhpuU7Wa9FcYd6Ebnhg",
                        "fIzeM58TDyJN6daUHJ3o7g"
                    ],
                    "gc": 0.490522422561257,
                    "ab": 1052440,
                    "loc": "complement(1050278..1052440)",
                    "neighborsId": [
                        1113396,
                        1113397,
                        1113398,
                        1113399,
                        1113400,
                        1113401,
                        1113402,
                        1113403,
                        1113404,
                        1113405,
                        1113406,
                        1113407,
                        1113408,
                        1113409,
                        1113410
                    ],
                    "cid": 633,
                    "p": {
                        "pd": "chemotaxis protein CheA",
                        "ac": "NP_233475.1",
                        "l": 720,
                        "aid": "AtBd-8mQw2yws-ZrpoEQWA",
                        "n": "similar to GB:M34669 SP:P07363 PID:145529 PID:145530 GB:U00096; identified by sequence similarity",
                        "x": {
                            "gi": 15601844
                        }
                    }
                }
            ],
            [
                {
                    "pos": 1368,
                    "gid": 65,
                    "a": 1594111,
                    "p": {
                        "aid": "BhwWfyWEGmZ3vyzusAAn3g",
                        "pd": "chemotaxis protein chea",
                        "ac": "YP_003375864.1",
                        "x": {
                            "gi": 285018153
                        },
                        "ec": "2.7.13.3",
                        "l": 666,
                        "f": [
                            [
                                "function",
                                "cell processes go:0009987; motility (incl. chemotaxis, energytaxis, aerotaxis, redoxtaxis) go:0042330"
                            ],
                            [
                                "inference",
                                "ab initio prediction:AutomaticValidationOfFrameDPredictionUsedWithbl a stxHitsSPTR"
                            ]
                        ]
                    },
                    "st": "-",
                    "cid": 133,
                    "l": 2001,
                    "gc": 0.646676661669165,
                    "lo": "XALc_1369",
                    "_id": 241292,
                    "ns": [
                        "cheA1"
                    ],
                    "b": 1592111,
                    "aa": 1592111,
                    "gs_id": "17p1jLVtEiNmBkev7qMUpA",
                    "loc": "complement(1592111..1594111)",
                    "ps": False,
                    "ab": 1594111
                },
                {
                    "ab": 1594111,
                    "b": 1592111,
                    "a": 1594111,
                    "ps": False,
                    "l": 2001,
                    "gid": 65,
                    "gs_id": "17p1jLVtEiNmBkev7qMUpA",
                    "lo": "XALc_1369",
                    "_id": 241292,
                    "aa": 1592111,
                    "neighborsAC": [
                        "YP_003375857.1",
                        "YP_003375858.1",
                        "YP_003375859.1",
                        "YP_003375860.1",
                        "YP_003375861.1",
                        "YP_003375862.1",
                        "YP_003375863.1",
                        "YP_003375864.1",
                        "YP_003375865.1",
                        "YP_003375866.1",
                        "YP_003375867.1",
                        "YP_003375868.1",
                        "YP_003375869.1",
                        "YP_003375870.1",
                        "YP_003375871.1"
                    ],
                    "ns": [
                        "cheA1"
                    ],
                    "pos": 1368,
                    "neighborsAseq": [
                        "vUAr5zs1PNT3_HmchtLDhw",
                        "gtDTgyPqWtBhOe7OvuQSxg",
                        "X5kaYb0eM1g2oePd9Vvs5Q",
                        "1VUoZ2RskO3Vza5ugyKhFw",
                        "9UguluIOplPu5F8fD6aUKw",
                        "qPySXeWOkgpYI9wTg5ozww",
                        "iBWZ0ms8qlrWP1hEmSrwQg",
                        "BhwWfyWEGmZ3vyzusAAn3g",
                        "B2S6964apP1AmxbFDVm7bQ",
                        "dLAsA9_MLAaMdtn7fyfzqQ",
                        "yJ_ChjwChhs8-vn7O4JpkQ",
                        "YgmGCDWP5oztQgfZd8Oitg",
                        "Q6ekvhB4N2reTt1Xz5cwGQ",
                        "moA_Wlbmqv4MnW2NB0SlXQ",
                        "o3o3B3ed_YKzjINKSuDgeg"
                    ],
                    "gc": 0.646676661669165,
                    "st": "-",
                    "loc": "complement(1592111..1594111)",
                    "neighborsId": [
                        241285,
                        241286,
                        241287,
                        241288,
                        241289,
                        241290,
                        241291,
                        241292,
                        241293,
                        241294,
                        241295,
                        241296,
                        241297,
                        241298,
                        241299
                    ],
                    "cid": 133,
                    "p": {
                        "pd": "chemotaxis protein chea",
                        "ac": "YP_003375864.1",
                        "l": 666,
                        "aid": "BhwWfyWEGmZ3vyzusAAn3g",
                        "f": [
                        [
                            "function",
                            "cell processes go:0009987; motility (incl. chemotaxis, energytaxis, aerotaxis, redoxtaxis) go:0042330"
                        ],
                        [
                            "inference",
                            "ab initio prediction:AutomaticValidationOfFrameDPredictionUsedWithbl a stxHitsSPTR"
                        ]
                        ],
                        "ec": "2.7.13.3",
                        "x": {
                            "gi": 285018153
                        }
                    }
                }
            ]
        ]
        client = bitk3.get_mist22_client()
        mist22 = client.mist22

        for gene, expect in fixtures:
            gene = buildCheABRmakeFasta._getNeighbors(
                mist22,
                gene,
                geneNeighborhoodWindow
            )

            assert gene == expect

        client.close()

    def test__getSigTransInfoOfNeighbors(self):
        fixtures = [
            [
                {
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
                },
                {
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
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': 'thisChea'},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': 'cheb'},
                        {'che': False},
                        {'che': False},
                        {'che': False}
                    ]
                }
            ],
            [
                {
                    "b": 1050278,
                    "a": 1052440,
                    "ps": False,
                    "l": 2163,
                    "gid": 319,
                    "gs_id": "ohUew9_0Y21HMztN2UkYBw",
                    "pos": 1006,
                    "aa": 1050278,
                    "_id": 1113403,
                    "lo": "VCA1095",
                    "neighborsAC": [
                        "NP_233469.1",
                        None,
                        "NP_233470.1",
                        "NP_233471.1",
                        "NP_233472.1",
                        "NP_233473.1",
                        "NP_233474.1",
                        "NP_233475.1",
                        "NP_233476.1",
                        "NP_233477.1",
                        "NP_233478.1",
                        "NP_233479.1",
                        "NP_233480.1",
                        "NP_233481.1",
                        "NP_233482.1"
                    ],
                    "st": "-",
                    "neighborsAseq": [
                        "QpctL60-ANyBk1aGp3W-hA",
                        None,
                        "jc8tmp9XtR6lTPMX3oSAlg",
                        "xj70OvOV0Ok4UHkJQriVvA",
                        "ConBlJawh0S9eNcwpUdStA",
                        "aXjcWvPeespvQZp3aap5LQ",
                        "qyHHt7NIZE1oStynXSwaLA",
                        "AtBd-8mQw2yws-ZrpoEQWA",
                        "MBG9AyB9JZ0g2yFduk70_w",
                        "6msMWTG01YKLCBo5dqqL4A",
                        "gA7NzCfTVtXu1kVp8hPZWg",
                        "x0OMwjtP2hGSk9pS-3btXg",
                        "d4_i8687rdGbPEaz-nj7_w",
                        "8DbYhpuU7Wa9FcYd6Ebnhg",
                        "fIzeM58TDyJN6daUHJ3o7g"
                    ],
                    "gc": 0.490522422561257,
                    "ab": 1052440,
                    "loc": "complement(1050278..1052440)",
                    "neighborsId": [
                        1113396,
                        1113397,
                        1113398,
                        1113399,
                        1113400,
                        1113401,
                        1113402,
                        1113403,
                        1113404,
                        1113405,
                        1113406,
                        1113407,
                        1113408,
                        1113409,
                        1113410
                    ],
                    "cid": 633,
                    "p": {
                        "pd": "chemotaxis protein CheA",
                        "ac": "NP_233475.1",
                        "l": 720,
                        "aid": "AtBd-8mQw2yws-ZrpoEQWA",
                        "n": "similar to GB:M34669 SP:P07363 PID:145529 PID:145530 GB:U00096; identified by sequence similarity",
                        "x": {
                            "gi": 15601844
                        }
                    }
                },
                {
                    "b": 1050278,
                    "a": 1052440,
                    "ps": False,
                    "l": 2163,
                    "gid": 319,
                    "gs_id": "ohUew9_0Y21HMztN2UkYBw",
                    "pos": 1006,
                    "aa": 1050278,
                    "_id": 1113403,
                    "lo": "VCA1095",
                    "neighborsAC": [
                        "NP_233469.1",
                        None,
                        "NP_233470.1",
                        "NP_233471.1",
                        "NP_233472.1",
                        "NP_233473.1",
                        "NP_233474.1",
                        "NP_233475.1",
                        "NP_233476.1",
                        "NP_233477.1",
                        "NP_233478.1",
                        "NP_233479.1",
                        "NP_233480.1",
                        "NP_233481.1",
                        "NP_233482.1"
                    ],
                    "st": "-",
                    "neighborsAseq": [
                        "QpctL60-ANyBk1aGp3W-hA",
                        None,
                        "jc8tmp9XtR6lTPMX3oSAlg",
                        "xj70OvOV0Ok4UHkJQriVvA",
                        "ConBlJawh0S9eNcwpUdStA",
                        "aXjcWvPeespvQZp3aap5LQ",
                        "qyHHt7NIZE1oStynXSwaLA",
                        "AtBd-8mQw2yws-ZrpoEQWA",
                        "MBG9AyB9JZ0g2yFduk70_w",
                        "6msMWTG01YKLCBo5dqqL4A",
                        "gA7NzCfTVtXu1kVp8hPZWg",
                        "x0OMwjtP2hGSk9pS-3btXg",
                        "d4_i8687rdGbPEaz-nj7_w",
                        "8DbYhpuU7Wa9FcYd6Ebnhg",
                        "fIzeM58TDyJN6daUHJ3o7g"
                    ],
                    "gc": 0.490522422561257,
                    "ab": 1052440,
                    "loc": "complement(1050278..1052440)",
                    "neighborsId": [
                        1113396,
                        1113397,
                        1113398,
                        1113399,
                        1113400,
                        1113401,
                        1113402,
                        1113403,
                        1113404,
                        1113405,
                        1113406,
                        1113407,
                        1113408,
                        1113409,
                        1113410
                    ],
                    "cid": 633,
                    "p": {
                        "pd": "chemotaxis protein CheA",
                        "ac": "NP_233475.1",
                        "l": 720,
                        "aid": "AtBd-8mQw2yws-ZrpoEQWA",
                        "n": "similar to GB:M34669 SP:P07363 PID:145529 PID:145530 GB:U00096; identified by sequence similarity",
                        "x": {
                            "gi": 15601844
                        }
                    },
                    'cheInfo': [
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': 'cher'},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': 'thisChea'},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False}
                    ]
                }
            ],
            [
                {
                    "ab": 1594111,
                    "b": 1592111,
                    "a": 1594111,
                    "ps": False,
                    "l": 2001,
                    "gid": 65,
                    "gs_id": "17p1jLVtEiNmBkev7qMUpA",
                    "lo": "XALc_1369",
                    "_id": 241292,
                    "aa": 1592111,
                    "neighborsAC": [
                        "YP_003375857.1",
                        "YP_003375858.1",
                        "YP_003375859.1",
                        "YP_003375860.1",
                        "YP_003375861.1",
                        "YP_003375862.1",
                        "YP_003375863.1",
                        "YP_003375864.1",
                        "YP_003375865.1",
                        "YP_003375866.1",
                        "YP_003375867.1",
                        "YP_003375868.1",
                        "YP_003375869.1",
                        "YP_003375870.1",
                        "YP_003375871.1"
                    ],
                    "ns": [
                        "cheA1"
                    ],
                    "pos": 1368,
                    "neighborsAseq": [
                        "vUAr5zs1PNT3_HmchtLDhw",
                        "gtDTgyPqWtBhOe7OvuQSxg",
                        "X5kaYb0eM1g2oePd9Vvs5Q",
                        "1VUoZ2RskO3Vza5ugyKhFw",
                        "9UguluIOplPu5F8fD6aUKw",
                        "qPySXeWOkgpYI9wTg5ozww",
                        "iBWZ0ms8qlrWP1hEmSrwQg",
                        "BhwWfyWEGmZ3vyzusAAn3g",
                        "B2S6964apP1AmxbFDVm7bQ",
                        "dLAsA9_MLAaMdtn7fyfzqQ",
                        "yJ_ChjwChhs8-vn7O4JpkQ",
                        "YgmGCDWP5oztQgfZd8Oitg",
                        "Q6ekvhB4N2reTt1Xz5cwGQ",
                        "moA_Wlbmqv4MnW2NB0SlXQ",
                        "o3o3B3ed_YKzjINKSuDgeg"
                    ],
                    "gc": 0.646676661669165,
                    "st": "-",
                    "loc": "complement(1592111..1594111)",
                    "neighborsId": [
                        241285,
                        241286,
                        241287,
                        241288,
                        241289,
                        241290,
                        241291,
                        241292,
                        241293,
                        241294,
                        241295,
                        241296,
                        241297,
                        241298,
                        241299
                    ],
                    "cid": 133,
                    "p": {
                        "pd": "chemotaxis protein chea",
                        "ac": "YP_003375864.1",
                        "l": 666,
                        "aid": "BhwWfyWEGmZ3vyzusAAn3g",
                        "f": [
                        [
                            "function",
                            "cell processes go:0009987; motility (incl. chemotaxis, energytaxis, aerotaxis, redoxtaxis) go:0042330"
                        ],
                        [
                            "inference",
                            "ab initio prediction:AutomaticValidationOfFrameDPredictionUsedWithbl a stxHitsSPTR"
                        ]
                        ],
                        "ec": "2.7.13.3",
                        "x": {
                            "gi": 285018153
                        }
                    }
                },
                {
                    "ab": 1594111,
                    "b": 1592111,
                    "a": 1594111,
                    "ps": False,
                    "l": 2001,
                    "gid": 65,
                    "gs_id": "17p1jLVtEiNmBkev7qMUpA",
                    "lo": "XALc_1369",
                    "_id": 241292,
                    "aa": 1592111,
                    "neighborsAC": [
                        "YP_003375857.1",
                        "YP_003375858.1",
                        "YP_003375859.1",
                        "YP_003375860.1",
                        "YP_003375861.1",
                        "YP_003375862.1",
                        "YP_003375863.1",
                        "YP_003375864.1",
                        "YP_003375865.1",
                        "YP_003375866.1",
                        "YP_003375867.1",
                        "YP_003375868.1",
                        "YP_003375869.1",
                        "YP_003375870.1",
                        "YP_003375871.1"
                    ],
                    "ns": [
                        "cheA1"
                    ],
                    "pos": 1368,
                    "neighborsAseq": [
                        "vUAr5zs1PNT3_HmchtLDhw",
                        "gtDTgyPqWtBhOe7OvuQSxg",
                        "X5kaYb0eM1g2oePd9Vvs5Q",
                        "1VUoZ2RskO3Vza5ugyKhFw",
                        "9UguluIOplPu5F8fD6aUKw",
                        "qPySXeWOkgpYI9wTg5ozww",
                        "iBWZ0ms8qlrWP1hEmSrwQg",
                        "BhwWfyWEGmZ3vyzusAAn3g",
                        "B2S6964apP1AmxbFDVm7bQ",
                        "dLAsA9_MLAaMdtn7fyfzqQ",
                        "yJ_ChjwChhs8-vn7O4JpkQ",
                        "YgmGCDWP5oztQgfZd8Oitg",
                        "Q6ekvhB4N2reTt1Xz5cwGQ",
                        "moA_Wlbmqv4MnW2NB0SlXQ",
                        "o3o3B3ed_YKzjINKSuDgeg"
                    ],
                    "gc": 0.646676661669165,
                    "st": "-",
                    "loc": "complement(1592111..1594111)",
                    "neighborsId": [
                        241285,
                        241286,
                        241287,
                        241288,
                        241289,
                        241290,
                        241291,
                        241292,
                        241293,
                        241294,
                        241295,
                        241296,
                        241297,
                        241298,
                        241299
                    ],
                    "cid": 133,
                    "p": {
                        "pd": "chemotaxis protein chea",
                        "ac": "YP_003375864.1",
                        "l": 666,
                        "aid": "BhwWfyWEGmZ3vyzusAAn3g",
                        "f": [
                            [
                                "function",
                                "cell processes go:0009987; motility (incl. chemotaxis, energytaxis, aerotaxis, redoxtaxis) go:0042330"
                            ],
                            [
                                "inference",
                                "ab initio prediction:AutomaticValidationOfFrameDPredictionUsedWithbl a stxHitsSPTR"
                            ]
                        ],
                        "ec": "2.7.13.3",
                        "x": {
                            "gi": 285018153
                        }
                    },
                    'cheInfo': [
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': 'thisChea'},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': False},
                        {'che': 'chea'}
                    ]
                }
            ]
        ]

        client = bitk3.get_mist22_client()
        mist22 = client.mist22

        for gene, expected in fixtures:
            gene = buildCheABRmakeFasta._getSigTransInfoOfNeighbors(mist22, gene)

            # assert gene['cheInfo'] == expected['cheInfo']
            assert gene == expected

        client.close()
    
