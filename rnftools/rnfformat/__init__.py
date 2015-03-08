#! /usr/bin/env python3

from .formatter import *
from .read import *
from .segment import *

__all__=[
	"Segment",
	"Read",
	"formatter",
]


if __name__ == "__main__":
	print ("Test of the module")

	print ()
	print ("1) SEGMENT")
	print ()

	for segment_repr in ["(0,0,N,0,0)","(0,1,F,56,59)","(0,1,R,01,59)"]:
		print("Original segment string", segment_repr)
		segment=Segment()
		segment.destringize(segment_repr)
		print("Segment after destringization", segment.stringize())

	print ()
	print ("2) READ")
	print ()


	read=Read()

	for read_name_test in [
				"__000324a__(2,3,R,34,3643),(4,3,F,1,56)__",
				"__00000000__(01,1,F,0390501,0000000)__"
			]:
		read.destringize(read_name_test)
		print()
		print("Original read string", read_name_test)
		print("Read after destringization", read.stringize ())
