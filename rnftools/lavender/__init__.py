# -*- coding: utf-8 -*-
"""
LAVEnder
~~~~~~~~

The LAVEnder Next-Generation Sequencing mappers evaluation tool.

:copyright: Copyright 2015 by Karel Brinda.
:license: MIT, see LICENSE for details.

.. module:: lavender
	:platform: Unix
	:synopsis: The LAVEnder Next-Generation Sequencing mappers evaluation tool.

.. moduleauthor:: Karel Brinda <karel.brinda@gmail.com>

"""

import lavender
import os

__docformat__ = 'reStructuredText'

# todo:
#  - left-right: cigar correction
#  - threshold as a parameter

include = os.path.join( os.path.dirname(__file__), "lavender.snake")
input = []


lavender._PANNELS_ = []
lavender._BAMS_    = []
lavender._REPORTS_ = []

from lavender.report import *
from lavender.pannel import *
from lavender.bam import *
