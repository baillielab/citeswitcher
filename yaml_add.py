'''
Add or append new yaml to every file in a directory
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
parser.add_argument('-y', '--newyaml',  default='{"projects": ["isaric4c", "odap"]}',   help='specify new yaml in json format')
args = parser.parse_args()
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#-----------------------------
for filename in os.listdir(args.dirpath):
	thispath = os.path.join(args.dirpath,filename)
	if filename[0] in ["_", "."] or os.path.isdir(thispath):
		continue
	with open(thispath) as f:
		text = f.read()
		h,r = citefunctions.readheader(text)
	workingyaml = citefunctions.mergeyaml(citefunctions.getyaml(h), json.loads(args.newyaml))
	with open(thispath,"w") as o:
		o.write('---\n{}\n---'.format(yaml.dump(workingyaml)).replace("\n\n","\n"))
		o.write(r)
