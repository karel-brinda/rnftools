
from .sample import *

from .source import *

from .artIllumina import *
from .cuReSim import *
from .dwgSim import *
from .wgSim import *

import os

__all__=[
	"ArtIllumina",
	"CuReSim",
	"DwgSim",
	"WgSim",
]

include = os.path.join( os.path.dirname(__file__), "mishmash.snake")

input=[]

_SAMPLES_ = []
_SOURCES_ = []

"""
	Create a new sample
"""
def sample(name,ends):
	if len(mishmash._SAMPLES_) > 0 and len(mishmash._SAMPLES_[-1].get_sources()) == 0:
		del mishmash._SAMPLES_[-1]

	if name in [sample.get_name() for sample in mishmash._SAMPLES_]:
		raise ValueError ("More samples have the same name. Each sample must have a unique name.")

	Sample(name=name,ends=ends)
	input=[sample.fq_fns() for sample in mishmash._SAMPLES_]

sample("sample",ends=1)

