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
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
#-------------------
import bib
bib.init()
import include
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
parser.add_argument('-y', '--yaml', action="store_true", default=False, help='use default yaml file with pandoc')
parser.add_argument('-i', '--include', action="store_false", default=True, help='do NOT include files')
parser.add_argument('-l', '--localbibonly', action="store_true", default=False, help='use only local bib file')
parser.add_argument('-r', '--customreplace', action="store_true", default=False, help='run custom find/replace commands specified in config file')
parser.add_argument('-m', '--messy', action="store_true", default=False, help='disable clean up of intermediate files')
parser.add_argument('-d', '--outputsubdir',    help='outputdir - always a subdir of the working directory', default=config['outputsubdirname'])
parser.add_argument('-b', '--bibfile',    help='bibfile', default=config['default_bibfile'])
parser.add_argument('-u', '--updatebibfile',    help='bibfile', default=config['default_updatebibfile'])
parser.add_argument('-c', '--cslfile',    help='csl citation styles file', default=config['cslfile'])
parser.add_argument('-img', '--imagedir',    help='imagedirectoryname', default=config['imagedir'])
args = parser.parse_args()
#-------------------
args.bibfile = os.path.expanduser(args.bibfile)
args.updatebibfile = os.path.expanduser(args.updatebibfile)
args.filepath = os.path.expanduser(args.filepath)
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
outpath, filename = os.path.split(args.filepath)
outpath = os.path.join(outpath, args.outputsubdir)
if outpath != '':
    citefunctions.check_dir(outpath)
filestem = '.'.join(filename.split('.')[:-1])
bibout = os.path.join(outpath, filestem+".bib")
yamlinstruction = ""
if args.yaml:
    yamlinstruction = os.path.abspath(os.path.join(scriptpath, config['yamlfile']))
    yamlsource = yamlinstruction
else:
    yamlsource = args.filepath
yamldata = citefunctions.getyaml(yamlsource)
if 'bibliography' in yamldata.keys():
    bibout = yamldata['bibliography']
    print ('using yaml-specified bib from {}: {}'.format(yamlsource, yamldata['bibliography']))
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
    lines = include.parse_includes(args.filepath)
    text = ''.join(lines)
else:
    with io.open(args.filepath, "r", encoding="utf-8") as f:
        text = f.read()
text = citefunctions.make_unicode(text)
if len(args.pandoc_outputs)>0 or args.stripcomments:
    text = include.stripcomments(text)
if args.verbose:
    print ("read input file: ", args.filepath)
#-------------------
if args.localbibonly:
    bib.full_bibdat = citefunctions.read_bib_files([bibout])
else:
    bib.full_bibdat = citefunctions.read_bib_files([args.bibfile, bibout])
bib.make_alt_dicts()
if args.verbose:
    print ("read bib files: ", [args.bibfile, bibout])
#-----------------
# replace the ids in the text with the outputstyle
text = citefunctions.replace_blocks(text, args.outputstyle)
#-----------------
# save remote bibliography
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
# save new text file
cslpath = os.path.abspath(os.path.join(os.path.dirname(__file__), args.cslfile))
text = citefunctions.addheader(text, os.path.abspath(bibout), cslpath)
with io.open(outputfile, 'w', encoding='utf-8') as file:
    file.write(text+"\n\n")
#-----------------
if len(args.pandoc_outputs)>0:
    # make symlink to image dir (saves the hassle of converting all refs)
    home_img = os.path.abspath(os.path.join(os.path.split(args.filepath)[0], args.imagedir))
    subdir_img = os.path.abspath(os.path.join(outpath,args.imagedir))
    if os.path.exists(home_img) and not os.path.exists(subdir_img):
        cmd = "ln -s {} {}".format(home_img, subdir_img)
        print (cmd)
        subprocess.call(cmd, shell=True)
    # run pandoc
    workingdir, filename = os.path.split(os.path.abspath(outputfile))
    os.chdir(workingdir)
    for thisformat in args.pandoc_outputs:
        citefunctions.callpandoc(filename, '.{}'.format(thisformat), yaml=yamlinstruction)
    if len(args.pandoc_outputs)>0 and not args.messy:
        #then clean up because the intermediate files are probably not wanted
        cmd = "rm {}".format(outputfile)
        print (cmd)
        subprocess.call(cmd, shell=True)
    '''
    citefunctions.callpandoc(filename, '.pdf', yaml=yamlinstruction)
    #citefunctions.callpandoc(filename, '.tex', yaml=yamlinstruction)
    #citefunctions.callpandoc(filename, '.html', yaml=yamlinstruction)
    #citefunctions.callpandoc(filename, '.txt', yaml=yamlinstruction)
    '''




















