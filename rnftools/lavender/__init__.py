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

import os

__docformat__ = 'reStructuredText'
__all__=["Bam","Report","Pannel"]

# todo:
#  - left-right: cigar correction
#  - threshold as a parameter

include = os.path.join( os.path.dirname(__file__), "lavender.snake")
input = []


_PANNELS_ = []
_BAMS_    = []
_REPORTS_ = []

from .report import *
from .pannel import *
from .bam import *


