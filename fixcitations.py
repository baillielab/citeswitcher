#!/usr/bin/env python
# encoding: utf-8

#-------------------
import os
import io
import sys
import json
import copy
import subprocess
#-------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath, 'python-bibtexparser-master/'))
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
sys.path.append(os.path.join(scriptpath, 'dependencies/'))
import oyaml as yaml
#-------------------
import bib
bib.init()
import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
#---- essential arguments
#parser.add_argument('-filepath', default="../genomicc-grant/caseforsupport.qmd", help='filepath') # - *** essential ***
#---- additional files to specify
parser.add_argument('-b', '--bibfile', default=config['default_bibfile'], help='bibfile')
parser.add_argument('-y', '--yaml', default='auto', help='use this yaml file as a source; use "normal" or "fancy" to use templates')
#---- other options
parser.add_argument('-w', '--wholereference', action="store_true", default=False, help='try to match whole references.')
parser.add_argument('-o', '--outputstyle', type=str, choices=['md','markdown','tex','latex','pubmed','pmid','inline'], default='null', help='output references format')
parser.add_argument('-ow', '--overwrite', action="store_true", default=False, help='overwrite input file with new markdown version')
parser.add_argument('-l', '--localbibonly', action="store_true", default=False, help='use only local bib file')
parser.add_argument('-flc', '--force_lowercase_citations', action="store_true", default=False, help='force all citation references into lowercase')
parser.add_argument('-d', '--outputsubdir', default=config['outputsubdirname'], help='outputdir - always a subdir of the working directory')
parser.add_argument('-v', '--verbose', action="store_true", default=False)
args = parser.parse_args()
#-------------------
'''
if args.filepath:
    args.filepath = os.path.abspath(os.path.expanduser(args.filepath))
else:
    print ("\nNo input file specified. Try a test file using this command:")
    print ("python fixcitations.py -f {} -o md".format(os.path.abspath(config['testfile'])))
    print ("and navigate to this directory to see the output: {}\n".format(os.path.split(config['testfile'])[0]))
    sys.exit()
'''
#-------------------
# manual settings for pythonista
args.filepath = "../genomicc-manuscript/manuscript.tex"
args.outputstyle = "tex"

#-------------------
sourcepath, filename = os.path.split(args.filepath)
outpath = os.path.join(sourcepath, args.outputsubdir)
if outpath != '':
    citefunctions.check_dir(outpath)
filestem = '.'.join(filename.split('.')[:-1])
#-------------------
# read input file
with io.open(args.filepath, "r", encoding="utf-8") as f:
	  text = f.read()
# Find bib, csl, yaml and template files
#-------------------
# ===> READ YAML according the YAML hierarchy: infile, local.yaml, other
if not args.yaml.endswith(".yaml"):
    args.yaml = args.yaml+".yaml"
yamlfile = os.path.join(sourcepath, filestem+".yaml")
infileyaml = citefunctions.getyaml(text) # READ FROM INPUT FILE FIRST
workingyaml = copy.copy(infileyaml)
if os.path.exists(yamlfile):
    with open(yamlfile) as f:
        workingyaml = citefunctions.mergeyaml(workingyaml, citefunctions.getyaml(f.read()))
if args.yaml in os.listdir(config['yamldir']):
    with open(os.path.join(config['yamldir'], args.yaml)) as f:
        workingyaml = citefunctions.mergeyaml(workingyaml, citefunctions.getyaml(f.read()))
# ===> CSL - hierarchy - yaml-specified, sup-files
if 'csl' in workingyaml:
    if not workingyaml['csl'].endswith('.csl'):
        workingyaml['csl'] = workingyaml['csl'] + '.csl'
    if workingyaml['csl'].startswith("http"):
        pass
    elif not os.path.exists(workingyaml['csl']):
        # try to find it in the sup-files folder
        newcsl = os.path.relpath(os.path.join(config['csldir'], os.path.split(workingyaml['csl'])[-1]))
        if os.path.exists(newcsl):
            workingyaml['csl'] = newcsl
    if 'csl' not in workingyaml:
        workingyaml['csl'] = os.path.relpath(os.path.join(config["csldir"], config["csldefault"])) # HARD OVERWRITE instead of merge
# ===> BIB - just read them all and copy into one local version
args.bibfile = os.path.abspath(os.path.expanduser(args.bibfile))
if 'bibliography' in workingyaml.keys():
    print ('using yaml-specified bib: {}'.format(workingyaml['bibliography']))
else:
    workingyaml['bibliography'] = config['default_localbib']  # HARD OVERWRITE instead of merge
print ("Using {} as bibout".format(workingyaml['bibliography']))
#-------------------
if args.verbose:
    print("Filepath:", args.filepath)
input_file_extension = args.filepath.split('.')[-1]
if args.filepath.endswith(".md") or args.filepath.endswith(".txt") or args.filepath.endswith(".qmd"):
    if args.outputstyle=='null':
        args.outputstyle = 'md'
elif args.filepath.endswith(".tex"):
    if args.outputstyle=='null':
        args.outputstyle = 'tex'
#-------------------
# name output file
if args.overwrite:
    print("overwriting original file")
    citelabel = "."
elif args.outputstyle == 'pubmed' or args.outputstyle == 'pmid':
    print("outputstyle = pmid")
    citelabel = ".citepmid."
elif args.outputstyle == 'markdown' or args.outputstyle == 'md':
    print("outputstyle = md")
    citelabel = ".citemd."
elif args.outputstyle == 'latex' or args.outputstyle == 'tex':
    print("outputstyle = tex")
    citelabel = ".citetex."
elif args.outputstyle == 'inline':
    print("outputstyle = inline")
outputfile = os.path.join(outpath, filestem+citelabel+input_file_extension)
#-------------------
text = citefunctions.make_unicode(text)
text = citefunctions.readheader(text)[1]
if args.verbose:
    print ("read input file: ", args.filepath)
#-------------------
localbibpath = os.path.join(sourcepath, workingyaml['bibliography'])
bib.db = citefunctions.read_bib_files([localbibpath])
if args.localbibonly:
    print ("\n*** reading local bib file only: {} ***\n".format(localbibpath))
    bib.full_bibdat = bib.db
else:
    if args.verbose: print ("reading bibfiles:", args.bibfile, localbibpath)
    bib.full_bibdat = citefunctions.read_bib_files([args.bibfile, localbibpath])
if args.force_lowercase_citations:
    print ("forcing lowercase citations")
    bib.id_to_lower()
bib.make_alt_dicts()
#-----------------
# replace the ids in the text with the outputstyle
text = citefunctions.replace_blocks(text, args.outputstyle, use_whole=args.wholereference, flc=args.force_lowercase_citations)
#-----------------
# save bibliography
print ('\nsaving bibliography for this file here:', localbibpath)
outbib = bibtexparser.dumps(bib.db)
outbib = citefunctions.make_unicode(outbib)
with open(localbibpath, "w", encoding="utf-8") as bf:
    bf.write(outbib)
#-----------------
# save new text file
with io.open(outputfile, 'w', encoding='utf-8') as file:
    if args.outputstyle =="md":
        file.write('---\n{}\n---'.format(yaml.dump(workingyaml)).replace("\n\n","\n"))
    file.write(text+"\n\n")
    print ("outputfile:", outputfile)

