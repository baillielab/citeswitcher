#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# encoding: utf-8

#-------------------
import os
import io
import sys
import oyaml as yaml
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
#-------------------
import bib
bib.init()
import include
import wordcount
import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
#---- essential arguments
parser.add_argument('filepath', default=None, help='filepath') # - *** essential ***
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
#---- options for the pandoc wrapper (unnecessary in the quarto era)
parser.add_argument('-p', '--pandoc_outputs', action='append', default=[], help='append as many pandoc formats as you want: pdf docx html txt md tex')
parser.add_argument('-img', '--imagedir', default=config['imagedir'], help='imagedirectoryname')
parser.add_argument('-i', '--include', action="store_true", default=True, help='include files')
parser.add_argument('-r', '--replace', action="store_true", default=True, help='replace using yaml or json instructions')
parser.add_argument('-m', '--messy', action="store_true", default=False, help='disable clean up of intermediate files')
parser.add_argument('-csl', '--autocsl', action="store_true", default=False, help='automatically create csl yaml heading')
parser.add_argument('-mf', '--move_figures', action="store_true", default=False, help='move all figures to the end and create captions section for submission to journal')
parser.add_argument('-pm', '--pandoc_mermaid', action="store_true", default=False, help='use pandoc-mermaid-filter')
parser.add_argument('-lt', '--latex_template', default=None, help='a latex template to be passed to pandoc') # e.g. a letter format
parser.add_argument('-wt', '--word_template', default=None, help='a word template to be passed to pandoc') # e.g. a letter format
parser.add_argument('-ptp', '--pathtopandoc', default='pandoc', help='specify a particular path to pandoc if desired')
parser.add_argument('-redact', '--redact', action="store_true", default=False, help='redact between <!-- STARTREDACT --> <!-- ENDREDACT --> tags')
parser.add_argument('-s', '--stripcomments', action="store_true", default=False, help='stripcomments in html format')
parser.add_argument('-v', '--verbose', action="store_true", default=False, help='verbose')
parser.add_argument('-x', '--xelatex', action="store_true", default=False, help='use xelatex in pandoc build')
parser.add_argument('-ch', '--chaptermode', action="store_true", default=False, help='use pandoc --top-level-division=chapter')
parser.add_argument('-svg', '--convert_svg', action="store_true", default=False, help='convert svg images to pdf - replaces any pdf files with the same name')
parser.add_argument('-tf', '--tableformat', default="pipe", help='fancy_grid, github, grid, html, jira, latex, latex_booktabs, latex_raw, mediawiki, moinmoin, orgtbl, pipe, plain, presto, psql, rst, simple, textile, tsv, youtrack')
args = parser.parse_args()
#-------------------
if args.filepath:
    args.filepath = os.path.abspath(os.path.expanduser(args.filepath))
else:
    print ("\nNo input file specified. Try a test file using this command:")
    print ("python fixcitations.py -f {} -o md".format(os.path.abspath(config['testfile'])))
    print ("and navigate to this directory to see the output: {}\n".format(os.path.split(config['testfile'])[0]))
    sys.exit()
sourcepath, filename = os.path.split(args.filepath)
outpath = os.path.join(sourcepath, args.outputsubdir)
if outpath != '':
    citefunctions.check_dir(outpath)
filestem = '.'.join(filename.split('.')[:-1])
#-------------------
# read input file
if args.include:
    text = include.parse_includes(args.filepath, tbf=args.tableformat)
    #text = ''.join(lines)
else:
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
if args.autocsl:
    if 'csl' not in workingyaml:
        workingyaml['csl'] = os.path.relpath(os.path.join(config["csldir"], config["csldefault"])) # HARD OVERWRITE instead of merge
# ===> BIB - just read them all and copy into one local version
args.bibfile = os.path.abspath(os.path.expanduser(args.bibfile))
if 'bibliography' in workingyaml.keys():
    print ('using yaml-specified bib: {}'.format(workingyaml['bibliography']))
else:
    workingyaml['bibliography'] = config['default_localbib']  # HARD OVERWRITE instead of merge
print ("Using {} as bibout".format(workingyaml['bibliography']))
# ===> TEMPLATES hierarchy: command argument, yaml
for thistemplatetype in ["latex_template", "word_template"]:
    if thistemplatetype in workingyaml:
        print ("template in yaml")
        if not vars(args)[thistemplatetype]: # this is the same as args.latex_template or args.word_template
            vars(args)[thistemplatetype] = workingyaml[thistemplatetype]
            print ("template assigned from YAML: {}".format(workingyaml[thistemplatetype]))
    if vars(args)[thistemplatetype]:
        if not os.path.exists(vars(args)[thistemplatetype]): # then look in sup dir
            vars(args)[thistemplatetype] = os.path.split(vars(args)[thistemplatetype])[-1] # filename
            cdir = config["{}dir".format(thistemplatetype)]
            print ("***searching for {}".format(cdir))
            if os.path.exists(cdir):
                print ("cdir found. args: ", vars(args)[thistemplatetype])
                if vars(args)[thistemplatetype] in os.listdir(cdir):
                    vars(args)[thistemplatetype] = os.path.join(cdir, vars(args)[thistemplatetype])
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
if args.redact:
    text = include.redact(text)
if args.stripcomments:
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
outbib = bibtexparser.dumps(bib.db)
outbib = citefunctions.make_unicode(outbib)
with open(localbibpath, 'w') as bf:
    bf.write(outbib)
#-----------------
'''
if len(bib.newpmids)>0:
    print ("\nURL for batch import to reference manager [CMD+dbl_click]:")
    print ("https://www.ncbi.nlm.nih.gov/pubmed/?term={}\n".format\
        ('&term='.join(["{}[pmid]".format(x) for x in pmids_to_add_to_inputdb])))
'''
#-----------------
if args.replace:
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
    if args.outputstyle =="md":
        file.write('---\n{}\n---'.format(yaml.dump(workingyaml)).replace("\n\n","\n"))
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
            if args.latex_template:
                pargstring+="--template {}".format(args.latex_template) # e.g. a letter format
                #args.xelatex = True # FORCE xelatex
            if args.pandoc_mermaid:
                pargstring+="--filter pandoc-mermaid"
                #args.xelatex = True # FORCE xelatex
        elif thisformat == "docx":
            if args.word_template:
                pargstring+="--reference-doc {}".format(args.word_template) # e.g. a letter format
        yamlfiles = []
        if os.path.exists(yamlfile):
            yamlfiles.append(yamlfile)
        if args.yaml in os.listdir(config['yamldir']):
            yamlfiles.append(args.yaml)
        citefunctions.callpandoc(os.path.abspath(outputfile),
                    sourcepath,
                    '.{}'.format(thisformat),
                    pargs=pargstring,
                    yaml=" ".join(yamlfiles),
                    out_dir=outpath,
                    x=args.xelatex,
                    ch=args.chaptermode,
                    pathtopandoc=args.pathtopandoc
                    )
    if len(args.pandoc_outputs)>0 and not args.messy and not args.overwrite:
        #then clean up because the intermediate files are probably not wanted
        cmd = "rm {}".format(outputfile)
        print (cmd)
        subprocess.call(cmd, shell=True)



















