#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# encoding: utf-8

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
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
#---- essential arguments
parser.add_argument('-f', '--filepath', default=None,                               help='filepath') # - *** essential ***
#---- additional files to specify
parser.add_argument('-b', '--bibfile', default=config['default_bibfile'],           help='bibfile')
parser.add_argument('-c', '--cslfile', default=config['cslfile'],                   help='csl citation styles file', )
parser.add_argument('-y', '--yaml', default='choose',                               help='use this yaml file with pandoc; use "normal" or "fancy"\
                    for config set one; the default, "choose", indicates that normal will be used only if there is no yaml in the current file')
#---- other options
parser.add_argument('-o', '--outputstyle', type=str, choices=['md','markdown','tex','latex','pubmed','pmid','inline'], default='null', help='output references format')
parser.add_argument('-p', '--pandoc_outputs',    action='append', default=[],          help='append as many pandoc formats as you want: pdf docx html txt md tex')
parser.add_argument('-l', '--localbibonly', action="store_true", default=False,     help='use only local bib file')
parser.add_argument('-d', '--outputsubdir', default=config['outputsubdirname'],     help='outputdir - always a subdir of the working directory', )
parser.add_argument('-img', '--imagedir', default=config['imagedir'],               help='imagedirectoryname')
parser.add_argument('-i', '--include', action="store_false", default=True,          help='do NOT include files')
parser.add_argument('-m', '--messy', action="store_true", default=False,            help='disable clean up of intermediate files')
parser.add_argument('-mf', '--move_figures', action="store_true", default=False,    help='move all figures to the end and create captions section for submission to journal')
parser.add_argument('-pm', '--pandoc_mermaid', action="store_true", default=False,    help='use pandoc-mermaid-filter')
parser.add_argument('-lt', '--latex_template', default='', help='a latex template to be passed to pandoc')
parser.add_argument('-ptp', '--pathtopandoc', default='pandoc', help='specify a particular path to pandoc if desired') # e.g. a letter format
parser.add_argument('-redact', '--redact', action="store_true", default=False,      help='redact between <!-- STARTREDACT --> <!-- ENDREDACT --> tags')
parser.add_argument('-s', '--stripcomments', action="store_true", default=False,    help='stripcomments in html format')
parser.add_argument('-v', '--verbose', action="store_true", default=False,          help='verbose')
parser.add_argument('-x', '--xelatex', action="store_true", default=False,          help='use xelatex in pandoc build')
parser.add_argument('-w', '--wholereference', action="store_true", default=False,          help='try to match whole references.')
parser.add_argument('-ch', '--chaptermode', action="store_true", default=False,          help='use pandoc --top-level-division=chapter')
args = parser.parse_args()
#-------------------
if args.filepath:
    args.filepath = os.path.abspath(os.path.expanduser(args.filepath))
else:
    print ("\nNo input file specified. Try a test file using this command:")
    print ("python fixcitations.py -f {} -o md".format(os.path.abspath(config['testfile'])))
    print ("and navigate to this directory to see the output: {}\n".format(os.path.split(config['testfile'])[0]))
    sys.exit()
with citefunctions.cd(os.path.split(args.filepath)[0]):
    args.bibfile = os.path.abspath(os.path.expanduser(args.bibfile))
    sourcepath, filename = os.path.split(args.filepath)
    pandocoutpath = os.path.join(sourcepath, args.outputsubdir)
    if pandocoutpath != '':
        citefunctions.check_dir(pandocoutpath)
    filestem = '.'.join(filename.split('.')[:-1])
    bibout = os.path.join(sourcepath, filestem+".bib")
    yamlsource = args.filepath
    if args.yaml=="choose":
        if len(citefunctions.getyaml(args.filepath, do_includes=args.include)) == 0:
            # then there is no yaml in this input file. use normal.
            args.yaml = "normal"
    if not args.yaml.endswith(".yaml"):
        # then this isn't a user-defined yaml file. Search config files.
        configyamlpath = os.path.join(config['yamldir'], args.yaml+".yaml")
        if os.path.exists(configyamlpath):
            yamlsource = configyamlpath
    else:
        specyamlpath = os.path.abspath(os.path.join(sourcepath, args.yaml))
        if os.path.exists(specyamlpath): # read user-specified yaml file
            yamlsource = specyamlpath
        else:
            print ("YAML file ({}) not found.\nProceeding with in-file YAML.".format(configyamlpath))
    yamldata = citefunctions.getyaml(yamlsource, do_includes=args.include)
    if 'bibliography' in yamldata.keys():
        bibout = yamldata['bibliography']
        print ('using yaml-specified bib from {}: {}'.format(yamlsource, yamldata['bibliography']))
        yamlbib = True
    bibout = os.path.abspath(bibout)
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
if len(args.pandoc_outputs)>0 or args.stripcomments:
    text = include.stripcomments(text)
if args.verbose:
    print ("read input file: ", args.filepath)
#-------------------
bib.db = citefunctions.read_bib_files([bibout])
if args.localbibonly:
    print ("\n*** reading local bib file only: {} ***\n".format(bibout))
    bib.full_bibdat = bib.db
else:
    if args.verbose: print ("reading bibfiles:", args.bibfile, bibout)
    bib.full_bibdat = citefunctions.read_bib_files([args.bibfile, bibout])
bib.make_alt_dicts()
#-----------------
# replace the ids in the text with the outputstyle
text = citefunctions.replace_blocks(text, args.outputstyle, use_whole=args.wholereference)
#-----------------
# save remote bibliography
print ('\nsaving remote bibliography for this file here:', bibout)
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
replacefile = os.path.join(sourcepath, "replace.json")
if os.path.exists(replacefile):
    text = citefunctions.findreplace(text, replacefile)
#-----------------
print ("Word count: {}\n".format(wordcount.wordcount(text)))
#-----------------
# save new text file
cslpath = os.path.abspath(os.path.join(os.path.dirname(__file__), args.cslfile))
if input_file_extension in ['md', 'markdown']:
    text = citefunctions.addheader(text, os.path.abspath(bibout), cslpath)
    # check for "# References" header or equivalent
    lines = text.replace("\r","\n").strip().split("\n")
    if not lines[-1].startswith("#") and not lines[-1].startswith("="):
        # then the last line isn't a header, so add one.
        if len(bib.db.entries) > 0:
            text += "<!--automatically added by fixcitations-->\n\n\n# References"
with io.open(outputfile, 'w', encoding='utf-8') as file:
    file.write(text+"\n\n")
#-----------------
if len(args.pandoc_outputs)>0:
    # run pandoc
    workingdir, filename = os.path.split(os.path.abspath(outputfile))
    os.chdir(workingdir)
    yamlinstruction = ""
    if yamlsource != args.filepath:
        yamlinstruction = yamlsource #NB this is concatenated with the input file by pandoc
    for thisformat in args.pandoc_outputs:
        pargstring = ""
        if thisformat == "pdf" or thisformat =="tex":
            if len(args.latex_template)>0:
                pargstring+="--template {}".format(args.latex_template) # e.g. a letter format
                args.xelatex = True
            if args.pandoc_mermaid:
                pargstring+="--filter pandoc-mermaid"
                args.xelatex = True
        citefunctions.callpandoc(filename,
                    '.{}'.format(thisformat),
                    pargs=pargstring,
                    yaml=yamlinstruction,
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





















