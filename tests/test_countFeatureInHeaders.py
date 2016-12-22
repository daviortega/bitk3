# -*- coding: utf-8 -*-

"""
test_countFeatureInHeaders
----------------------------------

Tests for `countFeatureInHeaders` script.
"""

import pytest

from contextlib import contextmanager
from click.testing import CliRunner

from bitk3 import countFeatureInHeaders
from bitk3 import cli
import os

myPath = os.getcwd()
dataPath = myPath + '/sampledata/'
sampleFileStr = dataPath + 'fasta.with.bitk.tags.fa'
	
#@pytest.fixture
def test_countFeatureInHeaders():
	position = 3
	result = countFeatureInHeaders.main(sampleFileStr, position)
	assert result == {
		'A': 2,
		'B': 2
	}