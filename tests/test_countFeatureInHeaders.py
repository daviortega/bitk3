# -*- coding: utf-8 -*-

"""
test_countFeatureInHeaders
----------------------------------

Tests for `countFeatureInHeaders` script.
"""

import pytest

from bitk3 import countFeatureInHeaders
import os

myPath = os.getcwd()
dataPath = myPath + '/sampledata/'

	
#@pytest.fixture
def test_countFeatureInHeaders():
	sampleFileStr = dataPath + 'fasta.with.bitk.tags.fa'
	position = 3
	result = countFeatureInHeaders.main(sampleFileStr, position)
	assert result == {
		'A': 2,
		'B': 2
	}

def test_countFeatureInHeaders_passingAltSeparator():
	sampleFileStr = dataPath + 'fasta.with.bitk.alttags.fa'
	position = 3
	sep = '-'
	result = countFeatureInHeaders.main(sampleFileStr, position, sep)
	assert result == {
		'A': 2,
		'B': 2
	}