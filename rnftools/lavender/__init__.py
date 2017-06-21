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

__docformat__ = 'reStructuredText'
# __all__=["Bam","Report","Panel"]

DEFAULT_ALLOWED_DELTA = 5
MAXIMAL_MAPPING_QUALITY = 255


# todo:
#  - left-right: cigar correction
#  - threshold as a parameter

def include():
	return os.path.join(os.path.dirname(__file__), "lavender.snake")


__INPUT__ = []


def input():
	return __INPUT__


def add_input(input):
	__INPUT__.append(input)


__PANELS__ = []


def panels():
	return __PANELS__


def add_panel(panel):
	__PANELS__.append(panel)


__BAMS__ = []


def bams():
	return __BAMS__


def add_bam(bam):
	__BAMS__.append(bam)


__REPORTS__ = []


def reports():
	return __REPORTS__


def add_report(report):
	__REPORTS__.append(report)


def _default_gp_style_func(i, number):
	colors = ["red", "green", "blue", "goldenrod", "black"]
	color = colors[i % len(colors)]
	return 'set style line {i} lt 1 pt {i} lc rgb "{color}";'.format(color=color, i=i + 1)


from .Bam import *
from .Panel import *
from .Report import *


##########################

def _format_xxx(xxx):
	return xxx.format(
		M="$2",
		m="$3",
		w="$4",
		P="$5",
		U="$6",
		u="$7",
		T="$8",
		t="$9",
		x="$10",
		all="$11",
	)


def _default_gp_style(i):
	colors = ["#ff0000", "#00ff00", "#888888", "#0000ff", "#daa520", "#000000"]
	color = colors[i % len(colors)]
	return 'set style line {i} lt 1 pt {i} lc rgb "{color}";'.format(color=color, i=i + 1)
