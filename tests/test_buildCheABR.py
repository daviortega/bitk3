# -*- coding: utf-8 -*-

"""
test_buildCheABR.py
----------------------------------

Tests for `test_buildCheABR` script.
"""

from bitk3 import buildCheABR
import os

myPath = os.getcwd()
dataPath = myPath + '/sampledata/'


def test_buildCheABR():
    sampleFileStr = dataPath + 'cheaTags.txt'
    window = 5
    result = buildCheABR.main(sampleFileStr, window)
    assert result == 0
