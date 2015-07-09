import smbl

from .Sample import *

from .Source import *

from .ArtIllumina import *
from .MasonIllumina import *
from .CuReSim import *
from .DwgSim import *
from .WgSim import *

import os

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
		smbl.messages.error("No sample defined",program="RNFtools",subprogram="MIShmash",exception=ValueError)
	__SOURCES__.append(sample)

def sources():
	return __SOURCES__


"""
	Create a new sample
"""
def sample(name,reads_in_tuple):
	if name in [sample.get_name() for sample in __SAMPLES__]:
		smbl.messages.error("More samples have the same name. Each sample must have a unique name.",program="RNFtools",subprogram="MIShmash",exception=ValueError)

	Sample(name=name,reads_in_tuple=reads_in_tuple)
	add_input(current_sample().fq_fns())

