﻿#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# encoding: utf-8

'''

python fix_citations.py -f test/eme1.md -o pmid

Intended use:
Feed in a .md or .tex document
references in square [] or curly {} or round brackets will be searched for PMID

If there is a mixture of reference types in one set of parentheses, ids with the format PMID:NNNNN will be handled correctly but spaces etc will be ignored.



The default or specified .bib file will then be searched for the relevant citations.
Citations will be replaced with a citaion in either .md or .tex or PMID format.

Any citations not present in the master .bib file will be downloaded from pubmed and added to the supplementary.bib file

Two new files will be written:
1. a new .md or .tex file with correct citations
2. a local .bib file containing all cited reference details.


export PATH=$PATH":$HOME/Dropbox/3_scripts_and_programs/citeswitcher/"


sudo ln -s "~/Dropbox/3_scripts_and_programs/citeswitcher/fixcitations.py" /usr/local/bin/fix


'''

#-------------------
import os
import io
import sys
import subprocess
#-------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath, 'python-bibtexparser-master/'))
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
#-------------------
import bib
bib.init()
import include
import wordcount
import citefunctions
scriptpath = os.path.dirname(os.path.realpath(__file__))
config = citefunctions.getconfig(os.path.join(scriptpath, 'config.json'))
#-------------------
import argparse
parser = argparse.ArgumentParser()
# - essential
parser.add_argument('-f', '--filepath',    help='filepath', default=config['testfile'])
parser.add_argument('-o', '--outputstyle',    type=str, choices=['md','markdown','tex','latex','pubmed','pmid','inline'], default='null', help='output references format')
# - optional
parser.add_argument('-p', '--pandoc_outputs',    action='append', default=[], help='append as many pandoc formats as you want: pdf docx html txt md')
parser.add_argument('-s', '--stripcomments', action="store_true", default=False, help='stripcomments')
parser.add_argument('-v', '--verbose', action="store_true", default=False, help='verbose')
parser.add_argument('-y', '--yaml',  default="thisfile", help='use this yaml file with pandoc; use "default" for config set one')
parser.add_argument('-i', '--include', action="store_false", default=True, help='do NOT include files')
parser.add_argument('-l', '--localbibonly', action="store_true", default=False, help='use only local bib file')
parser.add_argument('-r', '--customreplace', action="store_true", default=False, help='run custom find/replace commands specified in config file')
parser.add_argument('-m', '--messy', action="store_true", default=False, help='disable clean up of intermediate files')
parser.add_argument('-d', '--outputsubdir',    help='outputdir - always a subdir of the working directory', default=config['outputsubdirname'])
parser.add_argument('-b', '--bibfile',    help='bibfile', default=config['default_bibfile'])
parser.add_argument('-c', '--cslfile',    help='csl citation styles file', default=config['cslfile'])
parser.add_argument('-img', '--imagedir',    help='imagedirectoryname', default=config['imagedir'])
args = parser.parse_args()
#-------------------
args.filepath = os.path.abspath(os.path.expanduser(args.filepath))
with citefunctions.cd(os.path.split(args.filepath)[0]):
    args.bibfile = os.path.abspath(os.path.expanduser(args.bibfile))
    outpath, filename = os.path.split(args.filepath)
    pandocoutpath = os.path.join(outpath, args.outputsubdir)
    if pandocoutpath != '':
        citefunctions.check_dir(pandocoutpath)
    filestem = '.'.join(filename.split('.')[:-1])
    bibout = os.path.join(outpath, filestem+".bib")
    yamlinstruction = ""
    if args.yaml == "default":
        yamlinstruction = os.path.abspath(os.path.join(scriptpath, config['yamlfile']))
        yamlsource = yamlinstruction
    elif args.yaml == "thisfile":
        yamlsource = args.filepath
    else:
        yamlinstruction = os.path.abspath(os.path.join(outpath, args.yaml))
        yamlsource = yamlinstruction
    yamldata = citefunctions.getyaml(yamlsource, do_includes=args.include)
    if 'bibliography' in yamldata.keys():
        bibout = yamldata['bibliography']
        print ('using yaml-specified bib from {}: {}'.format(yamlsource, yamldata['bibliography']))
        yamlbib = True
    bibout = os.path.abspath(bibout)
#-------------------
print("Filepath:", args.filepath)
filetype = "unknown filetype"
if args.filepath.endswith(".md") or args.filepath.endswith(".txt"):
    filetype="md"
    if args.outputstyle=='null':
        args.outputstyle = 'md'
elif args.filepath.endswith(".tex"):
    filetype="tex"
    if args.outputstyle=='null':
        args.outputstyle = 'tex'
#-------------------
# name output file
if args.outputstyle == 'pubmed' or args.outputstyle == 'pmid':
    print("outputstyle = pmid")
    outputfile = os.path.join(outpath, filestem+".citepmid."+filetype)
elif args.outputstyle == 'markdown' or args.outputstyle == 'md':
    print("outputstyle = md")
    outputfile = os.path.join(outpath, filestem+".citemd."+filetype)
elif args.outputstyle == 'latex' or args.outputstyle == 'tex':
    print("outputstyle = tex")
    outputfile = os.path.join(outpath, filestem+".citetex."+filetype)
elif args.outputstyle == 'inline':
    print("outputstyle = inline")
    outputfile = os.path.join(outpath, filestem+".citeinline."+filetype)
#-------------------
# read input file
if args.include:
    text = include.parse_includes(args.filepath)
    #text = ''.join(lines)
else:
    with io.open(args.filepath, "r", encoding="utf-8") as f:
        text = f.read()
text = citefunctions.make_unicode(text)
if len(args.pandoc_outputs)>0 or args.stripcomments:
    text = include.stripcomments(text)
if args.verbose:
    print ("read input file: ", args.filepath)
#-------------------
bib.db = citefunctions.read_bib_files([bibout])
if args.localbibonly:
    if args.verbose: print ("reading localbibonly")
    bib.full_bibdat = bib.db
else:
    if args.verbose: print ("reading bibfiles:", args.bibfile, bibout)
    bib.full_bibdat = citefunctions.read_bib_files([args.bibfile, bibout])
bib.make_alt_dicts()
#-----------------
# replace the ids in the text with the outputstyle
text = citefunctions.replace_blocks(text, args.outputstyle)
#-----------------
# save remote bibliography
print ('saving remote bibliography for this file here:', bibout)
with open(bibout, 'w') as bf:
    bibtexparser.dump(bib.db, bf)
#-----------------
'''
if len(bib.newpmids)>0:
    print ("\nURL for batch import to reference manager [CMD+dbl_click]:")
    print ("https://www.ncbi.nlm.nih.gov/pubmed/?term={}\n".format\
        ('&term='.join(["{}[pmid]".format(x) for x in pmids_to_add_to_inputdb])))
'''
#-----------------
if args.customreplace:
    text = citefunctions.findreplace(text, config['custom_find_replace'])
#-----------------
print ("Word count:", wordcount.wordcount(text))
#-----------------
# save new text file
cslpath = os.path.abspath(os.path.join(os.path.dirname(__file__), args.cslfile))
text = citefunctions.addheader(text, os.path.abspath(bibout), cslpath)
with io.open(outputfile, 'w', encoding='utf-8') as file:
    file.write(text+"\n\n")
#-----------------
if len(args.pandoc_outputs)>0:
    # run pandoc
    workingdir, filename = os.path.split(os.path.abspath(outputfile))
    os.chdir(workingdir)
    for thisformat in args.pandoc_outputs:
        citefunctions.callpandoc(filename, '.{}'.format(thisformat), yaml=yamlinstruction, out_dir=pandocoutpath)
    if len(args.pandoc_outputs)>0 and not args.messy:
        #then clean up because the intermediate files are probably not wanted
        cmd = "rm {}".format(outputfile)
        print (cmd)
        subprocess.call(cmd, shell=True)





















