import mishmash

from mishmash.sample import Sample

from mishmash.source import Source

from mishmash.artIllumina import ArtIllumina
from mishmash.cuReSim import CuReSim
from mishmash.dwgSim import DwgSim
from mishmash.wgSim import WgSim

import os

include = os.path.join( os.path.dirname(__file__), "mishmash.snake")

mishmash.input=[]

mishmash._SAMPLES_ = []
mishmash._SOURCES_ = []

"""
	Create a new sample
"""
def sample(name,ends):
	if len(mishmash._SAMPLES_) > 0 and len(mishmash._SAMPLES_[-1].get_sources()) == 0:
		del mishmash._SAMPLES_[-1]

	if name in [sample.get_name() for sample in mishmash._SAMPLES_]:
		raise ValueError ("More samples have the same name. Each sample must have a unique name.")

	mishmash.Sample(name=name,ends=ends)
	mishmash.input=[sample.fq_fns() for sample in mishmash._SAMPLES_]

mishmash.sample("sample",ends=1)

