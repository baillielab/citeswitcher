'''
Read yaml from every file in a directory
'''
import os
import json
import oyaml as yaml
#-------------------
import citefunctions
#-------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dirpath',  default=None,   help='whole directory')
parser.add_argument('-t', '--taglist',    action='append', default=["tags"], help='use this to append as many values as you want')
args = parser.parse_args()
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#-----------------------------

yamlfound = {}

for filename in os.listdir(args.dirpath):
	thispath = os.path.join(args.dirpath,filename)
	if filename[0] in ["_", "."] or os.path.isdir(thispath):
		continue
	with open(thispath) as f:
		text = f.read()
		h,r = citefunctions.readheader(text)
	y = citefunctions.getyaml(h)
	for t in args.taglist:
		try:
			y[t]
		except:
			continue
		taglist = y[t]
		if isinstance(taglist, str):
			taglist = [taglist] # make it a list to keep things simple
		if taglist:
			for tag in taglist:
				try:
					yamlfound[tag].append(filename)
				except:
					yamlfound[tag]=[]
					yamlfound[tag].append(filename)

for tag in sorted(yamlfound.keys()):
	print ("\t{}:\n\t\t{}\n\n".format(tag, "\n\t\t".join(yamlfound[tag])))


