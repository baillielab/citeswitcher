#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

import os
import re
import sys
import json
import oyaml as yaml
import select
import difflib
import calendar
import requests
import platform
import subprocess
from Bio import Entrez
from Bio import Medline
from lxml import etree
import xml.etree.ElementTree as ET
import latexchars
#-------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath, 'python-bibtexparser-master/'))
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
from bibtexparser.latexenc import latex_to_unicode
#-------------
import bib
import include
#-------------
def getconfig(cfgfile="null"):
    if cfgfile=="null":
        for cfgname in ['config_local.json', "config.json"]:
            cfgfile = os.path.join(scriptpath, cfgname)
            if os.path.exists(cfgfile):
                break
    with open(cfgfile) as json_data_file:
        data = json.load(json_data_file)
    for item in data:
        if type(data[item]) is str:
            if data[item].startswith("{{scriptpath}}"):
                data[item] = os.path.join(scriptpath, data[item].replace("{{scriptpath}}",""))
    try:
        Entrez.email = data['email']
    except:
        pass
    try:
        Entrez.tool = data['toolname']
    except:
        pass
    if "@" not in data['email']:
        print ("No email in config file: {}".format(cfgfile))
        sys.exit()
    return data
#-------------
markdown_labels_to_ignore = [
    '@fig:','@sec:','@tbl:','@eq:',
    '#fig:','#sec:','#tbl:','#eq:',
    '(\d+(\.\d+)?)(pt|em)',
    ]
#-------------
#regexes
squarebrackets = '\[[\s\S]+?\]' #'[^!]\[[^!][\s\S]+?\]'# '\[[\s\S]+?\]'
wordbrackets = '\[<sup>[\s\S]+?<\/sup>\]\(#[\s\S]+?\)' 
wordmultiplebrackets = '<sup>[\s\S]+?\)<\/sup>' 
curlybrackets = '\{[\s\S]+?\}' #'\{[[^!]\s\S]+?\}'
commentedsquarebrackets = '<!--\[[\s\S]+?\]-->'
roundbrackets = '\([\s\S]+?\)'
latexbrackets = '\\\cite\{[\s\S]+?\}'
#-------------
pubmedsearchstrings = ["PMID:", "PubMed:"]
mdsearchstrings = ["@"]
doisearchstrings = [
            "doi: https://doi.org/",
            "doi: http://doi.org/",
            "doi: http://dx.doi.org/",
            "doi: https://dx.doi.org/",
            "doi:https://doi.org/",
            "doi:http://doi.org/",
            "doi:http://dx.doi.org/",
            "doi:https://dx.doi.org/",
            "doi:",
            "https://doi.org/",
            "http://doi.org/",
            "http://dx.doi.org/",
            "https://dx.doi.org/",
        ]
#-------------
no_user_response_count = 0
#-------------
class cd:
    """Context manager for changing the current working directory"""
    """ by Brian M. Hunt https://stackoverflow.com/users/19212/brian-m-hunt"""
    '''
    use like this:
    with cd(scriptpath):
        ...nested code
    '''
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
#-------------
def fix_permissions(this_path):
    os.system("/bin/chmod 755 %s"%(this_path))

def check_dir(this_dir):
    if not os.path.isdir(this_dir):
        os.mkdir(this_dir)
    fix_permissions(this_dir)

def getext(filepath):
    return os.path.split(filepath)[-1].split('.')[-1]

def newext(filepath, thisext):
    return filepath[:filepath.rfind('.')] + thisext

def uncomment_images(thistext):
    commented_image_regex = '<!--!\[.*?\]\(.*?\).*?-->'
    ci_found = re.findall(commented_image_regex, thistext)
    print ("found: ", len(ci_found), len(thistext))
    for ci in ci_found:
        try:
            thistext = thistext.replace(ci, ci[4:-3])
        except:
            pass
    return thistext

def move_figures(thistext, dropfilepath=True):
    #figureformat = '!\[.*?\]\(.*?\).*?\n'
    figureformat = '!\[[\s\S]+?\]\([\s\S]+?\).*?\n'
    blank_svg = os.path.abspath(os.path.join(scriptpath, "sup-files/no_image.svg"))
    figures_found = list(set(re.findall(figureformat, thistext)))
    lines = [x for x in thistext.split("\n")]
    lines.reverse()
    for x in lines:
        # remove references line. why?
        if len(x.strip()) > 0:
            if x.startswith("#") and "Reference" in x:
                rline = x+"\n"
                thistext = thistext.replace(rline,"")
                print ("Removing this line:", rline)
            break
    if len(figures_found)>0:
        thistext += "# Figure Legends\n\n"
    for i,fig in enumerate(figures_found):
        thistext = thistext.replace(fig,"")
        if dropfilepath:
            startpos = fig.rfind("](")
            fig = fig[:startpos] + re.sub('\]\(.*?\)', ']({})'.format(blank_svg), fig[startpos:], count=0)
        legend = "{}\n\n".format(fig)
        thistext += legend
        thistext += "\n\n"
    return thistext

def find_svgs(thistext):
    imagelink = '!\[.*?\]\(.*?\.svg\)'
    svgs_found = re.findall(imagelink, thistext)
    return svgs_found

def replace_svgs(thistext, thispath):
    '''inkscape -D --export-filename=recruitmentbyweek.pdf recruitmentbyweek.svg '''
    #imagelink = '!\[[\S]*?\]\([\S]*?.svg\)'
    svgcalls = find_svgs(thistext)
    for call in list(set(svgcalls)):
        filename = call.rsplit(']', 1)[1]
        filename = re.findall('\([\S]*?.svg\)',filename)[0][1:-1]
        subdir, filename = os.path.split(filename)
        outpath = os.path.abspath(os.path.join(thispath, subdir))
        pdf_filename = filename.replace(".svg",".pdf")
        callpdf = call.replace(filename, pdf_filename)
        cmd = "inkscape -D --export-filename={} {}".format(pdf_filename, filename)
        print ("Replacing {} with pdf: {}".format(filename, pdf_filename))
        print (outpath)
        print (cmd)
        with cd(outpath):
            subprocess.call(cmd, shell=True)
        thistext = thistext.replace(call, callpdf)
    return thistext

def make_output(thispath, pathtopandoc="pandoc", outputformats=[], localbibonly=False):
    # run fix citations in messy mode
    extra_args = "-svg -flc "
    # -x indicates xelatex mode. Handles special characters. Crashes sid.
    # -flc indicates force lowercase citations
    if localbibonly:
        extra_args += (" -l")
    cmd = '{} {} {} -f {} -m {} -ptp {} '.format(
        sys.executable,
        os.path.join(scriptpath, "fixcitations.py"),
        extra_args,
        thispath,
        " ".join(["-p "+x.replace(".","") for x in outputformats]),
        pathtopandoc
        )
    with cd(scriptpath):
        print ("CWD:", os.getcwd())
        print (cmd)
        subprocess.call(cmd, shell=True)

def callpandoc(f, pandocworkingpath, out_ext, out_dir='', pargs="", yaml="", x=False, ch=False, pathtopandoc="pandoc"):
    # crossref must come before citeproc
    # entire yaml file is pre-pended to the main file
    cmd = '{} --resource-path={} --markdown-headings=atx --filter {}-crossref --citeproc {} {} {} -o {} '.format(
            pathtopandoc,
            pandocworkingpath,
            pathtopandoc,
            pargs,
            yaml,
            f,
            os.path.join(out_dir, os.path.split(newext(f, out_ext))[-1])
        )
    if x:
        cmd += " --pdf-engine=xelatex "
    if ch:
        cmd += " --top-level-division=chapter "
    if out_ext in ['.md','.txt']:
        #cmd += "  -t markdown-citations -t markdown-strict "
        cmd += "  -t markdown-citations "
    else:
        if out_ext not in ['.html']:
            cmd += " -s " # STANDALONE OUTPUT FOR EVERYTHING APART FROM MD/TXT/HTML OUT
    with cd(pandocworkingpath):
        print ("CWD:", os.getcwd())
        print (cmd)
        subprocess.call(cmd, shell=True)

def read_bib_files(bibfiles):
    bfs = ""
    for bibfile in bibfiles:
        if os.path.exists(bibfile):
            # read bibtex file
            try:
                size = os.path.getsize(bibfile)
                #print ("reading bib file ({}) {}".format(size, bibfile))
            except:
                print ("bib file not found at {}".format(bibfile))
                sys.exit()
            if size > 0:
                with open(bibfile) as bf:
                    bfs += bf.read()
            else:
                print ("bib file empty: {}".format(bibfile))
        else:
            print ("File does not exist: {}".format(bibfile))
    try:
        return bibtexparser.bparser.BibTexParser(common_strings=True, homogenize_fields=True, interpolate_strings=False).parse(bfs, partial=False)
    except:
        return BibDatabase()

def sort_db(thisdb, sortby="year"):
    sorter = {}
    for this_entry in thisdb.entries:
        s = this_entry[sortby]
        try:
            s = int(s)
        except:
            pass
        sorter[this_entry['ID']] = s
    theseids = [key for key, value in sorted(iter(sorter.items()), key=lambda k_v: (k_v[1],k_v[0]), reverse=True)]
    thisdb.entries = [thisdb.entries_dict[thisid] for thisid in theseids]

#-------------
def findreplace(inputtext, frdict):
    for f in frdict:
        inputtext = inputtext.replace(f, str(frdict[f]))
    return inputtext

#-------------
def getyaml(text):
    h,r = readheader(text)
    if h.startswith("---"):
        h = h[4:-3] # strip yaml identifiers
    yml = yaml.load(h, Loader=yaml.Loader)
    if yml:
        return yml
    else:
        return {}

def mergeyaml(priority_yaml, extra_yaml):
    # keep the contents of the first yaml if there is a conflict
    if not(priority_yaml):
        priority_yaml = {}
    if extra_yaml:
        for item in extra_yaml:
            try:
                priority_yaml[item]
            except:
                priority_yaml[item] = extra_yaml[item]
    return priority_yaml

def readheader(filecontents):
    '''
        Read a valid markdown yaml header
        Input is full text of file
        Returns list of header items + full text of remainder
    '''
    t = filecontents.strip()
    t = t.replace('\r','\n')
    h = ""
    remainder = filecontents
    lines = [x for x in t.split('\n')] # don't strip because indentation matters
    if lines[0].strip() == '---':
        h1 = re.findall( '---[\s\S]+?---',filecontents)
        h2 = re.findall( '---[\s\S]+?\.\.\.',filecontents)
        if len(h1)>0 and len(h2)>0:
            #print ("both yaml header formats match! Taking the shorter one")
            if len(h1[0]) < len(h2[0]):
                h=h1[0]
                #print ("Choosing ---/---\n", h)
            else:
                h=h2[0]
                #print ("Choosing ---/...\n", h)
        elif len(h1)>0:
            h = h1[0]
        elif len(h2)>0:
            h = h2[0]
        remainder = filecontents.replace(h,'')
    return h, remainder
#-------------

def flatten(thislist):
    listcount = [1 for x in thislist if type(x)==type([])]
    if sum(listcount)==0:
        return thislist
    else:
        return [item for sublist in thislist for item in sublist] #flatten list

nbspace = re.compile(u"\N{NO-BREAK SPACE}", re.IGNORECASE)
def make_unicode(inputstring):
    if type(inputstring) != str:
        inputstring =  inputstring.decode('utf-8')
    inputstring = nbspace.sub(" ", inputstring)
    inputstring = inputstring.replace(u"\u2003", " ") # nonbreaking space
    inputstring = inputstring.replace(u"\u0391", "$\alpha$") # alpha
    inputstring = inputstring.replace(u"\u03B2", "$\beta$") # beta
    inputstring = inputstring.replace(u"\u1D737", "$\beta$") # beta (maths)
    inputstring = inputstring.replace(u"\u2009", " ") # "thin space"
    inputstring = inputstring.replace(u"\ufeff", "") # sometimes appears at start of file
    return inputstring

def get_parenthesised(thistext,parentheses=[squarebrackets, commentedsquarebrackets, curlybrackets]):
    if type(parentheses)!=type([]):
        print ("Error, parentheses must be in a list")
    outlist = []
    for p in parentheses:
        try:
            outlist += re.findall(p,thistext)
        except:
            continue
    # now replace nested parentheses with the internally nested ones
    nested_out = []
    for item in outlist:
        o = item[1:]
        qout = []
        for p in parentheses:
            qout += re.findall(p, o)
        if len(qout) > 0:
            nested_out += qout
        else:
            nested_out.append(item)
    for thislabel in markdown_labels_to_ignore:
        for x in nested_out:
            if re.match(thislabel, x[1:]):
                if not x[1:].startswith(thislabel):
                    print ("regex picked up this match: {} {} which was missed by startswith".format(x[1:], thislabel))
        nested_out = [x for x in nested_out if not re.match(thislabel, x[1:])]
    return nested_out

def remove_parentheses(thistext):
    '''
        remove matched parentheses at both ends
        and remove newlines
    '''
    thistext = thistext.strip()
    for opener in ('(','[','{','\\cite{','<sup>'):
        if thistext.startswith(opener):
            thistext = thistext[len(opener):]
    for thisreg in [
            r'[\s\S]+?<\/sup>\]\(#ref-',
            r'\[[\s\S]+?\]\(#ref-',
            ]:
        thistext = re.sub(thisreg, "", thistext, count=1)
    for closer in ('}',']',')</sup>',')'):
        if thistext.endswith(closer):
            thistext = thistext[:-len(closer)]
    thistext = thistext.replace('\n',' ')
    thistext = thistext.replace('\r',' ')
    return thistext

def split_by_delimiter(this_string, delimiters=[";",","]):
    '''
        return flattened list of non-empty, stripped strings, split by delimiters
    '''
    thislist = [this_string]
    for d in delimiters:
        thislist = [e.split(d) for e in thislist]
        thislist = flatten(thislist)
    thislist = [e.strip() for e in thislist if e!=""]
    return thislist

# ==== get citation blocks functions ====

def get_mixed_citation_blocks(inputtext):
    '''
    return citation blocks for PMID, DOI and MD formats, and new input text with these removed
    '''
    confirmed_blocks = []
    for theseparetheses in [squarebrackets, curlybrackets]:
        for b in get_parenthesised(inputtext, [theseparetheses]):
            for crossreflabel in markdown_labels_to_ignore:
                if remove_parentheses(b).startswith(crossreflabel):
                    continue
            for stem in pubmedsearchstrings+mdsearchstrings+doisearchstrings:
                if stem.lower() in b.lower(): # case-insensitive match
                    confirmed_blocks.append(b)
                    # remove this block so that nested blocks are only counted once
                    inputtext = inputtext.replace(b, '---citationblockremoved---')
                    break
    return confirmed_blocks, inputtext

def get_latex_citation_blocks(inputtext):
    '''
    return latex-formatted citation formats, and new input text with these removed
    '''
    confirmed_blocks = []
    for b in get_parenthesised(inputtext, [latexbrackets]):
        confirmed_blocks.append(b)
        inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks, inputtext

def get_word_citation_blocks(inputtext):
    '''
    return word-formatted citation formats, and new input text with these removed
    '''
    confirmed_blocks = []
    for b in get_parenthesised(inputtext, [wordbrackets]):
        confirmed_blocks.append(b)
        inputtext = inputtext.replace(b, '---citationblockremoved---')
    for b in get_parenthesised(inputtext, [wordmultiplebrackets]):
        numrefs = len(re.findall(r'\[[\s\S]+?\]\(#[\s\S]+?\)',b))
        mw = b[5:-6].split(",")
        if len(mw)>1 and len(mw)==numrefs: # there must be more than one reference, and every single one must be in correct format
            confirmed_blocks.append(b)
            inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks, inputtext

def clean_searchstring(this_string):
    dirtychars = ['“','”','`','‘','’']
    for d in dirtychars:
        this_string = this_string.replace(d,'')
    return this_string.strip()

def get_searchstring_from_wholecitation(wc):
    '''
        return the longest sentence that doesn't look like a list of authors or a citation
    '''
    authorsep = '[A-Z, a-z, .], [A-Z]' # regex for [capital, comma, space, capital]
    lendict = {x:len(x) for x in remove_parentheses(wc).split('.')}
    for x in sorted(iter(lendict.items()), key=lambda k_v: (k_v[1],k_v[0]), reverse=True):
        if len(x[0]) > 0:
            authoriness = float(len(re.findall(authorsep, x[0]))/len(x[0]))
            if authoriness < 0.01:
                if ":" in x[0] and ";" in x[0]: # probably not a title
                    continue
                else:
                    return clean_searchstring(x[0])

def get_wholereference_citation_blocks(inputtext):
    '''
        finds wholecitations AND doi
        return a list of citation blocks, and new input text with these removed
    '''
    confirmed_blocks = []
    for theseparetheses in [squarebrackets, curlybrackets]:
        for b in get_parenthesised(inputtext, [theseparetheses]):
            b = b.replace('\n',' ')
            if "." in b or ":" in b:
                new_entry = parse_wholecitation_block(b)
                if new_entry is not None:
                    confirmed_blocks.append(str(b))
                    inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks, inputtext

def parse_citation_block(thisblock):
    results = {
        'pmids':[],
        'ids':[],
        'notfound':[],
    }
    theseids = split_by_delimiter(remove_parentheses(thisblock))
    for thisid in theseids:
        thisid = clean_id(thisid)
        if thisblock.startswith("\\cite"): #latex formatting
            results['ids'].append(thisid)
        elif thisblock.startswith("[<sup>"): #word formatting single citation
            results['ids'].append(thisid)
        elif thisblock.startswith("<sup>"): #word formatting multipe citations
            results['ids'].append(thisid)
        elif thisid.startswith('PMID:') and len(thisid)>5:
            results['pmids'].append(thisid.replace("PMID:","").strip())
        elif thisid.startswith('@') and len(thisid)>1:
            results['ids'].append(thisid.replace("@","").strip())
        elif thisid.startswith('DOI:') and len(thisid)>4:
            # DOI is not an output format. Translate to an output format here (PMID or MD)
            new_entry = findcitation(thisid.replace("DOI:","").strip(), 'doi')
            if new_entry:
                results['ids'].append(new_entry['ID'])
            else:
                results['notfound'].append(thisid)
        else:
            results['notfound'].append(thisid)
    return results

def casereplace(thistext, catchall, cleanversion):
    '''
    replace <catchall> text with <cleanversion> if it occurs at the beginning of the text
    '''
    catchall = re.compile("^"+re.escape(catchall.strip()), re.IGNORECASE)
    return catchall.sub(cleanversion, thistext)

def clean_id(thisid):
    thisid = thisid.replace(" ","").strip()
    for stem in pubmedsearchstrings:
        thisid = casereplace(thisid, stem, 'PMID:')
    for stem in mdsearchstrings:
        thisid = casereplace(thisid, stem, '@')
    for stem in doisearchstrings:
        thisid = casereplace(thisid, stem, 'DOI:')
    return thisid

askedalready = {}
def parse_wholecitation_block(thisblock):
    ''' take a block of text, and return two lists of ids '''
    try:
        return askedalready[thisblock]
    except:
        outids = []
        # try a doi search first in case there's a doi in there
        d = parse_citation_block(thisblock)
        if len(d['ids'])>0:
            print ("untested: returning ids:", d['ids'])
            outids = d['ids']
        else:
            # if that doesn't work, try a title search
            title = get_searchstring_from_wholecitation(thisblock)
            if title is not None:
                found = findcitation(title, 'title', additionalinfo=thisblock)
                if found is not None:
                    outids.append(found['ID'])
        askedalready[thisblock] = outids, []
        return outids, []  # there is only ever one

#------------
def findcitation(info, infotype='pmid', additionalinfo=''):
    '''
        search the bibdat first
        search online in a variety of ways
        if found, use bib.new to add it to global db
        return a single bibtex entry
        or None if not found or in doubt
    '''
    global no_user_response_count
    info = str(info.strip())
    if infotype == 'pmid':
        try:
            print ("searching for {} in bibdat".format(info))
            return bib.pmids[info]
        except:
            print ("{} not found in bib.pmids: ".format(info))
            pass
        pub = p2b([info])
        if len(pub) > 0:
            if pub[0] != 'null' and pub[0] != None:
                print ("PMID:{} found online".format(info))
                print (pub[0])
                bib.new(pub[0])
                return pub[0]
        print ("PMID:{} NOT FOUND ON PUBMED".format(info))
        return None
    elif infotype == 'doi':
        msg=""
        try:
            return bib.dois[info]
        except:
            msg += ("DOI not in bib file: {}".format(info))
        pub = search_pubmed(info, "doi")
        if len(pub) == 1:
            pubent = p2b(pub[0])
            if len(pubent) > 0:
                if pubent[0] != 'null' and pubent[0] != None:
                    print (msg + "but found in pubmed: {}".format(pubent[0]))
                    bib.new(pubent[0])
                    return pubent[0]
        print (msg + ", nor in pubmed")
        return None
    elif infotype == 'title':
        if no_user_response_count > 2:
            print ("No user response to previous 3 queries. Running in silent mode.")
            return None
        # should add code to search bibdat for title here ***TODO***
        print ("searching pubmed for this title:\t{}".format(info))
        pub = search_pubmed(info, "title")
        if len(pub) == 1:
            pmid = pub[0]
            pubent = p2b(pmid)
            if len(pubent) > 0:
                if pubent[0] != 'null' and pubent[0] != None:
                    question = "--------------\n\
New citation (PMID:{}) found in Pubmed. Please check that input is the same as the found citation for this reference block: \n\
\n{}\n\n{}\n{}\n\n{}\n\n\
Enter y/n within 10 seconds".format(
                        pmid,
                        additionalinfo,
                        "{:>12}:    {}".format('Input Title', info),
                        "{:>12}:    {}".format('Found Title', pubent[0]['Title']),
                        '\n'.join( ["{:>12}:    {}".format(x,pubent[0][x]) for x in pubent[0] if x != "Title"])
                        )
                    #q = input (question)
                    print (question)
                    i,o,e = select.select([sys.stdin],[],[],10) # 10 second timeout
                    if i==[] and o==[] and e==[]:
                        no_user_response_count += 1
                    if i:
                        q = sys.stdin.readline().strip()
                        q = q.strip().upper()
                        if q == "Y":
                            print ('--confirmed--')
                            bib.new(pubent[0])
                            return pubent[0]
        return None


def id2pmid(theseids, notpmidlist=[]):
    '''
        input is a list of ids
    '''
    pmidlist = []
    for thisid in theseids:
        # try to find this id in full_bibdat
        try:
            bib.full_bibdat.entries_dict[thisid]
        except:
            bestmatchingkey = find_similar_keys(thisid, bib.full_bibdat.entries_dict)
            print(("biblatex id not found in biblatex file: {}. Best match in database is {}.".format(thisid, bestmatchingkey)))
            continue
        # if it is found, try to get the pmid
        if 'PMID' in bib.full_bibdat.entries_dict[thisid]:
            pmidlist.append(bib.full_bibdat.entries_dict[thisid]['PMID'])
        elif 'pmid' in bib.full_bibdat.entries_dict[thisid]:
            pmidlist.append(bib.full_bibdat.entries_dict[thisid]['pmid'])
        else:
            print ("pmid not found in bib file: {}. Searching online...".format(thisid))
            try:
                new_entry = findcitation(bib.full_bibdat.entries_dict[thisid]['doi'], 'doi')
            except:
                notpmidlist.append(thisid) # hack
                new_entry = bib.full_bibdat.entries_dict[thisid]
            if new_entry is None:
                new_entry = findcitation(bib.full_bibdat.entries_dict[thisid]['title'], 'title')
                if new_entry is None:
                    notpmidlist.append(thisid)
                    continue
            try:
                pmidlist.append(new_entry['pmid'])
                bib.pmids[new_entry['pmid']] = new_entry
            except:
                notpmidlist.append(thisid)
    return pmidlist, notpmidlist

def pmid2id(thesepmids, others):
    outids = []
    missing_ids = others
    for pmid in thesepmids:
        try:
            outids.append(bib.pmids[pmid]['ID'])
            bib.cite([bib.pmids[pmid]['ID']])
        except:
            new_entry = findcitation(pmid, 'pmid')
            if new_entry is not None:
                outids.append(new_entry['ID'])
                bib.new(new_entry)
            else:
                missing_ids.append(pmid)
    return outids, missing_ids

#------------
def format_inline(thisid):
    au = bib.db.entries_dict[thisid]['Author'].split(" and ")
    if len(au)>1:
        au = au[0] + " et al"
    else:
        au = au[0]
    formatted_citation = "{}. {} {};{}:{}".format(
        au,
        bib.db.entries_dict[thisid]['Journal'].capitalize(),
        bib.db.entries_dict[thisid]['Year'],
        bib.db.entries_dict[thisid]['Volume'],
        bib.db.entries_dict[thisid]['Pages'],
        )
    return formatted_citation

#------------
def pmidout(pmidlist, notpmidlist):
    # make a blockstring
    blockstring = ''
    if len(pmidlist) > 0:
        # add to the outputdatabase if possible but since PMID is the output format, it doesn't matter if we can't find it
        for x in pmidlist:
            try:
                bib.pmids[x]
            except:
                continue
            bib.cite([bib.pmids[x]['ID']])
        blockstring = '[' + ', '.join(["PMID:{}".format(x) for x in pmidlist]) + ']'
        if len(notpmidlist) > 0:
            blockstring += '[*' + ', '.join(notpmidlist) + ']'
    else:
        blockstring = 'null'
    return blockstring

def mdout(theseids, thesemissing=[], outputstyle="md", flc=False):
    if len(theseids) == 0:
        return 'null'
    if flc:
        theseids = [x.lower() for x in theseids]
    # add to the outputdatabase
    bib.cite(theseids)
    # make a blockstring
    blockstring = ''
    if outputstyle == "md":
        blockstring += "[{}]".format('; '.join(["@{}".format(x) for x in theseids]))
        if len(thesemissing) > 0:
            blockstring += "[***{}]".format(', '.join(thesemissing))
    elif outputstyle == "tex":
        blockstring += "\\cite{%s}"%(', '.join(theseids))
        if len(thesemissing) > 0:
            blockstring += "[**{}]".format(', '.join(thesemissing))
    elif outputstyle == "inline":
        blockstring += "({})".format(', '.join([format_inline(x) for x in theseids]))
        if len(thesemissing) > 0:
            blockstring += "[***{}]".format(', '.join(thesemissing))
    return blockstring

def replace_blocks(thistext, outputstyle="md", use_whole=False, flc=False):
    # pmid first as they are the most likely to have errors
    workingtext = thistext
    p, workingtext = get_mixed_citation_blocks(workingtext)
    l, workingtext = [x for x in get_latex_citation_blocks(workingtext) if x not in p]
    w, workingtext = [x for x in get_word_citation_blocks(workingtext) if x not in p+l]
    if use_whole:
        r, workingtext = [x for x in get_wholereference_citation_blocks(workingtext) if x not in p+l+w]
    else:
        r=[]
    print ("Number blocks using pmid or doi or md:{}".format(len(p)))
    print ("Number blocks using latex:{}".format(len(l)))
    print ("Number blocks using wholeref:{}".format(len(r)))
    replacedict = {}
    # there are slightly different procedures for each reference type:
    for b in p+l+w:
        citedhere = parse_citation_block(b)
        if outputstyle == 'md' or outputstyle=='tex' or outputstyle=='inline':
            theseids, notfound = pmid2id(citedhere['pmids'], citedhere['notfound']) # ids added to bib
            theseids = list(set(theseids+citedhere['ids']))
            replacedict[b] = mdout(theseids, notfound, outputstyle, flc=flc)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(citedhere['ids'], citedhere['notfound']) # ids added to bib
            pm = list(set(pm+citedhere['pmids']))
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in r:
        theseids, theseothers = parse_wholecitation_block(b)
        if outputstyle == 'md' or outputstyle=='tex' or outputstyle=='inline':
            replacedict[b] = mdout(theseids, theseothers, outputstyle, flc=flc)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(theseids) # ids added to bib
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in replacedict:
        if replacedict[b] == 'null':
            print ("{:>80} ... left alone".format(b))
            continue
        print ("{:>80} ==> {}".format(b, replacedict[b]))
        thistext = thistext.replace(b, replacedict[b])
    return thistext

#-----------------

def find_similar_keys(this_string, thisdict):
    '''
        search the keys of thisdict for a string similar to this_string
    '''
    topscore = -1
    bestmatch = 'none'
    for i, key in enumerate(thisdict.keys()):
        sim = difflib.SequenceMatcher(None, this_string, key).ratio()
        if sim>topscore:
            bestmatch = key
            topscore = sim
    return bestmatch


#------ PUBMED FUCTIONS -------

def search_pubmed(search_string, restrictfields=""):
    '''return a list of pmids for articles found by a search string'''
    if search_string is None:
        return []
    if len(search_string.strip())==0:
        return []
    max_returns = 500
    days = "all"
    try:
        handle = Entrez.esearch(db="pubmed", term=search_string, retmax=max_returns, rettype="medline", field=restrictfields)
        r = Entrez.read(handle)
    except:
        print ("pubmed search failure")
        return []
    return r['IdList']

# --------------------

def p2b(pmidlist):
    ''' by Nick Loman '''

    if type(pmidlist) != list:
        pmidlist = [str(pmidlist)]

    ## Fetch XML data from Entrez.
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    url = '{}?db=pubmed&id={}&rettype=abstract'.format(efetch, ','.join(pmidlist))
    try:
        r = requests.get(url)
    except:
        return []
    ##print(r.text) # to examine the returned xml
    ## Loop over the PubMed IDs and parse the XML using https://docs.python.org/2/library/xml.etree.elementtree.html
    bibout = []
    par = etree.XMLParser(encoding='utf-8', recover=True)
    root = ET.fromstring(r.text, parser=par)
    if root:
        for PubmedArticle in root.iter('PubmedArticle'):
            PMID = PubmedArticle.find('./MedlineCitation/PMID')
            ISSN = PubmedArticle.find('./MedlineCitation/Article/Journal/ISSN')
            Volume = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/Volume')
            Issue = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/Issue')
            Year = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
            Month = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Month')
            Title = PubmedArticle.find('./MedlineCitation/Article/Journal/Title')
            ArticleTitle = PubmedArticle.find('./MedlineCitation/Article/ArticleTitle')
            MedlinePgn = PubmedArticle.find('./MedlineCitation/Article/Pagination/MedlinePgn')
            Abstract = PubmedArticle.find('./MedlineCitation/Article/Abstract/AbstractText')
            # jkb additions
            PMCID = None
            DOI = None
            theseids = PubmedArticle.findall('./PubmedData/ArticleIdList/ArticleId')
            for thisid in theseids:
                if thisid.attrib['IdType'] == 'pmc':
                    PMCID = thisid
                elif thisid.attrib['IdType'] == 'doi':
                    DOI = thisid
            # format author list
            authors = []
            for Author in PubmedArticle.iter('Author'):
                try:
                    LastName = Author.find('LastName').text
                    ForeName = Author.find('ForeName').text
                except AttributeError:  # e.g. CollectiveName
                    continue
                authors.append('{}, {}'.format(LastName, ForeName))
            ## Use InvestigatorList instead of AuthorList
            if len(authors) == 0:
                ## './MedlineCitation/Article/Journal/InvestigatorList'
                for Investigator in PubmedArticle.iter('Investigator'):
                    try:
                        LastName = Investigator.find('LastName').text
                        ForeName = Investigator.find('ForeName').text
                    except AttributeError:  # e.g. CollectiveName
                        continue
                    authors.append('{}, {}'.format(LastName, ForeName))
            if Year is None:
                _ = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate')
                Year = _.text[:4]
                Month = '{:02d}'.format(list(calendar.month_abbr).index(_.text[5:8]))
            else:
                Year = Year.text
                if Month is not None:
                    Month = Month.text
            '''
            try:
                for _ in (PMID.text, Volume.text, Title.text, ArticleTitle.text, MedlinePgn.text, Abstract.text, ''.join(authors)):
                    if _ is None:
                        continue
                    assert '{' not in _, _
                    assert '}' not in _, _
            except AttributeError:
                pass
            '''

            # make the bibtex formatted output.
            bib = {}
            if len(authors)>0:
                authorname = authors[0].split(',')[0]
            else:
                authorname = ''
            try:
                titlewords = [x for x in ArticleTitle.text.split(' ') if len(x)>3]
            except:
                print ("PUBMED ERROR - no article title for PMID:{}: {}".format(PMID.text, ArticleTitle.text))
                continue
            if len(titlewords)>2:
                titlestring = ''.join(titlewords[:3])
            elif len(titlewords)>0:
                titlestring = ''.join(titlewords)
            else:
                titlestring = ''
            if len(authorname+titlestring)==0:
                titlestring = "PMID{}_".format(PMID.text)
            new_id = '{}{}{}'.format(authorname, titlestring, Year).lower()
            new_id = re.sub(r'\W+', '', new_id)
            bib["ID"] = latexchars.replace_accents(new_id)
            bib["Author"] = ' and '.join(authors)
            bib["Title"] = ArticleTitle.text
            bib["Journal"] = Title.text
            bib["Year"] = Year
            if Volume is not None:
                bib["Volume"] = Volume.text
            if Issue is not None:
                bib["Number"] = Issue.text
            if MedlinePgn is not None:
                bib["Pages"] = MedlinePgn.text
            if Month is not None:
                bib["Month"] = Month
            # bib[""] = (' Abstract={{{}}},'.format(Abstract.text))
            if PMCID is not None:
                bib["pmcid"] = PMCID.text
            if DOI is not None:
                bib["doi"] = DOI.text
            if ISSN is not None:
                bib["ISSN"] = ISSN.text
            bib["pmid"] = PMID.text
            # always return clean latex
            bib = {d:latex_to_unicode(bib[d]) for d in bib.keys()}
            bibout.append(bib)
    return bibout

#---------------

def getlink(entry):
    '''
        return DOI link,  PMC link, or PMID link in that order
    '''
    try:
        url = "http://dx.doi.org/{}".format(entry['doi'])
    except:
        try:
            url = "https://www.ncbi.nlm.nih.gov/pmc/articles/{}/".format(entry['pmcid'])
        except:
            url = "http://www.ncbi.nlm.nih.gov/pubmed/{}".format(entry['pmid'])
    return url





