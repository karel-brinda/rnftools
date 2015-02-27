
from .sample import *

from .source import *

from .artIllumina import *
from .cuReSim import *
from .dwgSim import *
from .wgSim import *

from ..errors import *

import os

__all__=[
	"ArtIllumina",
	"CuReSim",
	"DwgSim",
	"WgSim",
]


__INCLUDE__=os.path.join( os.path.dirname(__file__), "mishmash.snake")

def include():
	return __INCLUDE__

__INPUT__=[]

def add_input(input):
	__INPUT__.append(input)

def input():
	return __INPUT__


__SAMPLES__ = []

def add_sample(sample):
	__SAMPLES__.append(sample)

def samples():
	return __SAMPLES__

def current_sample():
	return __SAMPLES__[-1]


__SOURCES__ = []

def add_source(sample):
	if len(samples()) == 0:
		die("No sample defined")
	__SOURCES__.append(sample)

def sources():
	return __SOURCES__


"""
	Create a new sample
"""
def sample(name,ends):
	if name in [sample.get_name() for sample in __SAMPLES__]:
		die("More samples have the same name. Each sample must have a unique name.")

	Sample(name=name,ends=ends)

	# FIX!!!!!!!!!!!!!!!
	#__INPUT__=[sample.fq_fns() for sample in samples()]
	add_input(current_sample().fq_fns())

