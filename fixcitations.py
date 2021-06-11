#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# encoding: utf-8

#-------------------
import os
import io
import sys
import yaml
import json
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
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
#---- essential arguments
parser.add_argument('-f', '--filepath', default=None,                               help='filepath') # - *** essential ***
#---- additional files to specify
parser.add_argument('-b', '--bibfile', default=config['default_bibfile'],           help='bibfile')
parser.add_argument('-y', '--yaml', default='auto'  ,                               help='use this yaml file as a source; use "normal" or "fancy" to use templates')
#---- other options
parser.add_argument('-o', '--outputstyle', type=str, choices=['md','markdown','tex','latex','pubmed','pmid','inline'], default='null', help='output references format')
parser.add_argument('-p', '--pandoc_outputs',    action='append', default=[],       help='append as many pandoc formats as you want: pdf docx html txt md tex')
parser.add_argument('-l', '--localbibonly', action="store_true", default=False,     help='use only local bib file')
parser.add_argument('-d', '--outputsubdir', default=config['outputsubdirname'],     help='outputdir - always a subdir of the working directory')
parser.add_argument('-img', '--imagedir', default=config['imagedir'],               help='imagedirectoryname')
parser.add_argument('-i', '--include', action="store_false", default=True,          help='do NOT include files')
parser.add_argument('-m', '--messy', action="store_true", default=False,            help='disable clean up of intermediate files')
parser.add_argument('-mf', '--move_figures', action="store_true", default=False,    help='move all figures to the end and create captions section for submission to journal')
parser.add_argument('-pm', '--pandoc_mermaid', action="store_true", default=False,    help='use pandoc-mermaid-filter')
parser.add_argument('-lt', '--latex_template', default='', help='a latex template to be passed to pandoc') # e.g. a letter format
parser.add_argument('-ptp', '--pathtopandoc', default='pandoc', help='specify a particular path to pandoc if desired')
parser.add_argument('-redact', '--redact', action="store_true", default=False,      help='redact between <!-- STARTREDACT --> <!-- ENDREDACT --> tags')
parser.add_argument('-s', '--stripcomments', action="store_true", default=False,    help='stripcomments in html format')
parser.add_argument('-v', '--verbose', action="store_true", default=False,          help='verbose')
parser.add_argument('-x', '--xelatex', action="store_true", default=False,          help='use xelatex in pandoc build')
parser.add_argument('-w', '--wholereference', action="store_true", default=False,   help='try to match whole references.')
parser.add_argument('-ch', '--chaptermode', action="store_true", default=False,     help='use pandoc --top-level-division=chapter')
parser.add_argument('-svg', '--convert_svg', action="store_true", default=False,          help='convert svg images to pdf - replaces any pdf files with the same name')
parser.add_argument('-ui', '--uncomment_images', action="store_true", default=False,          help='include images commented out using html syntax <!--![]()\{\} -->')
parser.add_argument('-flc', '--force_lowercase_citations', action="store_true", default=False,          help='force all citation references into lowercase')
args = parser.parse_args()
#-------------------
default_localbib = "cs.bib"
#-------------------
if args.filepath:
    args.filepath = os.path.abspath(os.path.expanduser(args.filepath))
else:
    print ("\nNo input file specified. Try a test file using this command:")
    print ("python fixcitations.py -f {} -o md".format(os.path.abspath(config['testfile'])))
    print ("and navigate to this directory to see the output: {}\n".format(os.path.split(config['testfile'])[0]))
    sys.exit()
sourcepath, filename = os.path.split(args.filepath)
pandocoutpath = os.path.join(sourcepath, args.outputsubdir)
if pandocoutpath != '':
    citefunctions.check_dir(pandocoutpath)
filestem = '.'.join(filename.split('.')[:-1])
# Find bib, csl and yaml files
# YAML according the YAML hierarchy: infile, local.yaml, other
yamlsources = ['.'.join(x.split(".")[:-1]) for x in os.listdir(config['yamldir'])]
yamlfile = os.path.join(sourcepath, filestem+".yaml")
workingyaml = citefunctions.getyaml(args.filepath) # READ FROM INPUT FILE FIRST
print (workingyaml)
if os.path.exists(yamlfile):
    workingyaml = citefunctions.mergeyaml(workingyaml, citefunctions.getyaml(yamlfile))
    print (workingyaml)
if args.yaml in yamlsources:
    workingyaml = citefunctions.mergeyaml(workingyaml, citefunctions.getyaml(os.path.join(config['yamldir'],args.yaml+".yaml")))
    print (workingyaml)
# CSL
if ('csl' not in workingyaml) or not os.path.exists(workingyaml['csl']):
    workingyaml['csl'] = config["cslfile"] # HARD OVERWRITE instead of merge
# BIB
args.bibfile = os.path.abspath(os.path.expanduser(args.bibfile))
if 'bibliography' in workingyaml.keys():
    print ('using yaml-specified bib: {}'.format(workingyaml['bibliography']))
else:
    workingyaml['bibliography'] = default_localbib  # HARD OVERWRITE instead of merge
print ("Using {} as bibout".format(workingyaml['bibliography']))
with open(yamlfile,"w") as o:
    o.write('---\n{}\n---'.format(yaml.dump(workingyaml)).replace("\n\n","\n"))

#-------------------
if len(args.pandoc_outputs) > 0:  # if -p is set
    if args.outputstyle not in ['null','md','markdown']: # and -o is not md
        print ("FAIL: both -o and -p options are set. If requesting pandoc output with -p, the citation format output (-o) must be markdown (it defaults to this).")
        sys.exit()
    args.outputstyle = 'md' # MANDATE THAT OUTPUTSTYLE IS MD
    if 'md' in args.pandoc_outputs:
        args.messy = True # mustn't delete the requested outputfile!
    print ("pandoc output types: {}".format(args.pandoc_outputs))
# IF NO PANDOC OUTPUT IS REQUESTED, guess at output style if it isn't already set.
if args.verbose:
    print("Filepath:", args.filepath)
input_file_extension = "unrecognisedinputfile"
if args.filepath.endswith(".md") or args.filepath.endswith(".txt"):
    input_file_extension="md"
    if args.outputstyle=='null':
        args.outputstyle = 'md'
elif args.filepath.endswith(".tex"):
    input_file_extension="tex"
    if args.outputstyle=='null':
        args.outputstyle = 'tex'
#-------------------
# name output file
if args.outputstyle == 'pubmed' or args.outputstyle == 'pmid':
    print("outputstyle = pmid")
    outputfile = os.path.join(sourcepath, filestem+".citepmid."+input_file_extension)
elif args.outputstyle == 'markdown' or args.outputstyle == 'md':
    print("outputstyle = md")
    outputfile = os.path.join(sourcepath, filestem+".citemd."+input_file_extension)
elif args.outputstyle == 'latex' or args.outputstyle == 'tex':
    print("outputstyle = tex")
    outputfile = os.path.join(sourcepath, filestem+".citetex."+input_file_extension)
elif args.outputstyle == 'inline':
    print("outputstyle = inline")
    outputfile = os.path.join(sourcepath, filestem+".citeinline."+input_file_extension)
#-------------------
# read input file
if args.include:
    text = include.parse_includes(args.filepath)
    #text = ''.join(lines)
else:
    with io.open(args.filepath, "r", encoding="utf-8") as f:
        text = f.read()
text = citefunctions.make_unicode(text)
if args.redact:
    text = include.redact(text)
if args.uncomment_images:
    text = citefunctions.uncomment_images(text)
if len(args.pandoc_outputs)>0 or args.stripcomments:
    text = include.stripcomments(text)
if args.move_figures:
    text = citefunctions.move_figures(text)
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
with open(localbibpath, 'w') as bf:
    bibtexparser.dump(bib.db, bf)
#-----------------
'''
if len(bib.newpmids)>0:
    print ("\nURL for batch import to reference manager [CMD+dbl_click]:")
    print ("https://www.ncbi.nlm.nih.gov/pubmed/?term={}\n".format\
        ('&term='.join(["{}[pmid]".format(x) for x in pmids_to_add_to_inputdb])))
'''
#-----------------
if "replace" in workingyaml.keys():
    text = citefunctions.findreplace(text, workingyaml["replace"])
replacefile = os.path.join(sourcepath, "replace.json")
if os.path.exists(replacefile):
    with open(replacefile) as rf:
        replacedic = json.load(rf)
    text = citefunctions.findreplace(text, replacedic)
#-----------------
print ("Word count: {}\n".format(wordcount.wordcount(text)))
#-----------------
svgs = citefunctions.find_svgs(text)
if len(svgs)>0:
    if args.convert_svg:
        text = citefunctions.replace_svgs(text, sourcepath)
    elif "pdf" in args.pandoc_outputs:
        print("svg files cannot be included in native form in pdf files:")
        for x in svgs:
            print ("\t{}".format(x))
        print("use the -svg argument to convert these to pdfs with the same name and include them\n")
#-----------------
# save new text file
with io.open(outputfile, 'w', encoding='utf-8') as file:
    file.write(text+"\n\n")
    print ("outputfile:", outputfile)
#-----------------
if len(args.pandoc_outputs)>0:
    # run pandoc
    workingdir, filename = os.path.split(os.path.abspath(outputfile))
    os.chdir(workingdir)
    for thisformat in args.pandoc_outputs:
        pargstring = ""
        if thisformat == "pdf" or thisformat =="tex":
            if len(args.latex_template)>0:
                pargstring+="--template {}".format(args.latex_template) # e.g. a letter format
                args.xelatex = True
            if args.pandoc_mermaid:
                pargstring+="--filter pandoc-mermaid"
                args.xelatex = True
        citefunctions.callpandoc(os.path.abspath(outputfile),
                    '.{}'.format(thisformat),
                    pargs=pargstring,
                    yaml=yamlfile,
                    out_dir=pandocoutpath,
                    x=args.xelatex,
                    ch=args.chaptermode,
                    pathtopandoc=args.pathtopandoc
                    )
    if len(args.pandoc_outputs)>0 and not args.messy:
        #then clean up because the intermediate files are probably not wanted
        cmd = "rm {}".format(outputfile)
        print (cmd)
        subprocess.call(cmd, shell=True)





















