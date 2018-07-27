#!/usr/bin/python
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


'''

#-------------------
import string,os,sys
import io
import bibtexparser
import subprocess
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
#-------------------
cwd = os.getcwd()
sys.path.append(cwd)
import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
# - essential
parser.add_argument('-f', '--filepath',    help='filepath', default=config['testfile'])
parser.add_argument('-o', '--outputstyle',    type=str, choices=['md','markdown','tex','latex','pubmed','pmid'], default='null', help='output references format')
# - optional
parser.add_argument('-p', '--pandoc', action="store_true", default=False, help='also run pandoc')
parser.add_argument('-r', '--customreplace', action="store_true", default=False, help='run custom find/replace commands specified in config file')
parser.add_argument('-d', '--outputsubdir',    help='outputdir - always a subdir of the working directory', default=config['outputsubdirname'])
parser.add_argument('-b', '--bibfile',    help='bibfile', default=config['default_bibfile'])
parser.add_argument('-u', '--updatebibfile',    help='bibfile', default=config['default_updatebibfile'])
parser.add_argument('-c', '--cslfile',    help='csl citation styles file', default=config['cslfile'])
parser.add_argument('-i', '--imagedir',    help='imagedirectoryname', default=config['imagedir'])
args = parser.parse_args()
#-------------------
# correct ~ for the right path to home dir on linux or mac
home_dir_abspath = citefunctions.get_home_dir()
args.filepath = args.filepath.replace("~",home_dir_abspath)
args.bibfile = args.bibfile.replace("~",home_dir_abspath)
args.updatebibfile = args.updatebibfile.replace("~",home_dir_abspath)
#-------------------
print(args.filepath)
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
citefunctions.check_dir(outpath)
filestem = '.'.join(filename.split('.')[:-1])
bibout = os.path.join(outpath, filestem+".bib")
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
#-------------------
# read input file
with io.open(args.filepath, "r", encoding="utf-8") as my_file:
     text = my_file.read()
text = citefunctions.make_unicode(text)
print ("read input file")
#-------------------
# read bibtex file
with open(args.bibfile) as bibtex_file:
    bibdat = bibtexparser.bparser.BibTexParser(common_strings=True, homogenize_fields=True).parse_file(bibtex_file)
print ("read bib file")
#-------------------
# make dictionaries
pmids, ids = citefunctions.make_dictionaries(bibdat.entries)
#-------------------
# find all citations
idcitations = citefunctions.get_all_id_citations(text)
pmidcitations = citefunctions.get_all_pmid_citations(text)
print(idcitations)
print(pmidcitations)
#-------------------
# make new output database
db = BibDatabase()
for thisid in idcitations:
    print(thisid)
    try:
        ids[thisid]
    except:
        bestmatchingkey = citefunctions.find_similar_keys(thisid, ids)
        print(("biblatex id not found in biblatex file: {}. Best match in database is {}.".format(thisid, bestmatchingkey)))
        continue
    citefunctions.bibadd(db,ids[thisid])
#-------------------
# search through all the cited papers in this file and add them to the db 
update = BibDatabase()
for thispmid in pmidcitations:
    try:
        pmids[thispmid]
    except:
        b = citefunctions.p2b(thispmid)
        if b != 'null' and b != None:
            pmids[thispmid] = b
            print ("PMID:{} found online".format(thispmid))
            citefunctions.bibadd(update,pmids[thispmid])
        else:
            print ("PMID:{} NOT FOUND ON PUBMED".format(thispmid))
            continue
    citefunctions.bibadd(db,pmids[thispmid])
#-----------------
# replace the ids in the text with the outputstyle
text = citefunctions.replace_blocks(text, pmids, ids, args.outputstyle)
#-----------------
# save remote bibliography
with open(bibout, 'w') as biblatex_file:
    bibtexparser.dump(db, biblatex_file)
# save update bibliography
with open(args.updatebibfile, 'w') as biblatex_file:
    bibtexparser.dump(update, biblatex_file)
#-----------------
if args.customreplace:
    text = citefunctions.findreplace(text, config['custom_find_replace'])
#-----------------
# save new text file
cslpath = os.path.abspath(os.path.join(os.path.dirname(__file__), args.cslfile))
text = citefunctions.addheader(text, os.path.abspath(bibout), cslpath)
with io.open(outputfile, 'w', encoding='utf-8') as file:
    file.write(text)
#-----------------
if args.pandoc:
    # make symlink to image dir (saves the hassle of converting all refs)
    home_img = os.path.abspath(os.path.join(os.path.split(args.filepath)[0], args.imagedir))
    subdir_img = os.path.abspath(os.path.join(outpath,args.imagedir))
    if not os.path.exists(subdir_img):
        cmd = "ln -s {} {}".format(home_img, subdir_img)
        print (cmd)
        subprocess.call(cmd, shell=True)
    # run pandoc
    pandocscript = os.path.join(os.path.dirname(__file__), 'run_pandoc.py')
    cmd = "python {} -f {}".format(pandocscript, outputfile)
    subprocess.call(cmd, shell=True)
























