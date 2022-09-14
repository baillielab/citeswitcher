#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

'''
search through all subdirs
find all "cs.bib" files
find any that don't overlap with main bib file
put them in an "update" file ready for import into zotero
'''

#-------------------
import os
import sys
#-------------------
import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bibfile',    help='bibfile', default=config['default_bibfile'])
parser.add_argument('-o', '--outbib',    help='update (output) file', default=config["update_bibfile"])
parser.add_argument('-d', '--directory',    help='directory to search', default=config["writing_directory"])
args = parser.parse_args()
bibfile = os.path.expanduser(args.bibfile)
outbib = os.path.expanduser(args.outbib)
wdir = os.path.expanduser(args.directory)
#-------------------
import io
import json
import subprocess
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath, 'python-bibtexparser-master/'))
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
#-------------------
citefiles = []
for root, dirs, files in os.walk(wdir):
    for file in files:
        if file == "cs.bib":
            newfile = os.path.join(root, file)
            citefiles.append(newfile)
            print(newfile)
updatedat = citefunctions.read_bib_files(citefiles)
#-----------------
bibdat = citefunctions.read_bib_files([bibfile])
for entry in updatedat.entries:
    

#-----------------
# save update bibliography
with open(outbib, 'w') as bf:
    bibtexparser.dump(updatedat, bf)
#-----------------
























