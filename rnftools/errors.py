import sys,os
from termcolor import colored, cprint

width=80

def die(message):
	whole_message = [
			width*"=",
			(int(width/9))*"  ERROR  ",
			"",
			message,
			"",
			"If it is a bug, please report it on http://github.com/karel-brinda/rnftools/",
			"or to karel.brinda@gmail.com. Thank you.",
			"",
			(int(width/9))*"  ERROR  ",
			width*"=",
		]
	whole_message=map(lambda x: x.ljust(width," "),whole_message)

	cprint(
		os.linesep.join(whole_message),
		"white",
		"on_blue",
		attrs=['bold'],
	)
	sys.exit(1)
