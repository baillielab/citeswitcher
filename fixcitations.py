#!/usr/bin/env python
# encoding: utf-8

import argparse
import calendar
import copy
import difflib
import hashlib
import io
import json
import os
import re
import requests
import select
import sys
import xml.etree.ElementTree as ET

scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath, 'dependencies/'))
import oyaml as yaml
import Entrez
import Medline
import latexchars

sys.path.append(os.path.join(scriptpath, 'dependencies/python-bibtexparser-master/'))
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
from bibtexparser.latexenc import latex_to_unicode
from bibtexparser.latexenc import string_to_latex

# patch for yaml error in pythonista on ipad
import collections
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", DeprecationWarning)
    if not hasattr(collections, 'Hashable'):
        import collections.abc
        collections.Hashable = collections.abc.Hashable

#---
default_localbibname = "cs.bib"
default_global_bibfile = '_bibfiles/lib.pmid.bib'
#---


additionaldicts = []
def init():
    global full_bibdat
    full_bibdat = BibDatabase()
    global db # cited
    db = BibDatabase()
    global pmids
    pmids = {}
    global dois
    dois = {}
    additionaldicts.append((pmids, "pmid"))
    additionaldicts.append((dois, "doi"))
init()

def id_to_lower():
    for entry in full_bibdat.entries:
        entry["ID"] = entry["ID"].lower()

def make_alt_dicts():
    for entry in full_bibdat.entries:
        for thisdict, thislabel in additionaldicts:
            try:
                entry[thislabel]
            except:
                continue
            try:
                thisdict[entry[thislabel]]
                if args.verbose:
                    print("duplicate {} in biblatex database:{}".format(thislabel, entry[thislabel]))
            except:
                pass
            thisdict[entry[thislabel]] = entry

def new(entry):
    '''
        add new reference to full_bibdat.entries, full_bibdat.entries_dict and all additionaldicts
    '''
    if entry is not None:
        try:
            entry['ENTRYTYPE']
        except:
            entry['ENTRYTYPE'] = 'article'
        already_in = False
        try:
            full_bibdat.entries_dict[entry['ID']] # existing full_bibdat entry ALWAYS takes precedence
            already_in = True
        except:
            entry = cleanbib(entry)
        if already_in:
            # add any additional fields from online by merging dictionaries
            entry = {**full_bibdat.entries_dict[entry['ID']], **entry}
            full_bibdat.entries_dict[entry['ID']] = entry
        else:
            print ("adding entry to bib:\n {}".format(entry['ID']))
            full_bibdat.entries_dict[entry['ID']] = entry
            full_bibdat.entries = [entry] + full_bibdat.entries
        for thisdict, thislabel in additionaldicts:
            try:
                entry[thislabel]
            except:
                continue
            try:
                thisdict[entry[thislabel]]
            except:
                thisdict[entry[thislabel]] = entry

def cite(theseids):
    if args.verbose:
        print ("cite function has been asked to find:\n",theseids)
    fails = []
    for thisid in theseids:
        try:
            db.entries_dict[thisid]
        except:
            try:
                full_bibdat.entries_dict[thisid]
            except:
                print ("id not found in full bibtex file: ", thisid)
                fails.append(thisid)
                continue
            db.entries = [full_bibdat.entries_dict[thisid]] + db.entries
            db.entries_dict[thisid] = full_bibdat.entries_dict[thisid]
    return fails

def cleanbib(bibtex_entry):
    return {d:string_to_latex(bibtex_entry[d]) for d in bibtex_entry.keys()}

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

markdown_labels_to_ignore = [
    '@fig:','@sec:','@tbl:','@eq:',
    '#fig:','#sec:','#tbl:','#eq:',
    r'(\d+(\.\d+)?)(pt|em)',
    ]

squarebrackets = r'\[[\s\S]+?\]' #'[^!]\[[^!][\s\S]+?\]'# '\[[\s\S]+?\]'
wordbrackets = r'\[<sup>[\s\S]+?<\/sup>\]\(#[\s\S]+?\)' 
wordmultiplebrackets = r'<sup>[\s\S]+?\)<\/sup>' 
curlybrackets = r'\{[\s\S]+?\}' #'\{[[^!]\s\S]+?\}'
commentedsquarebrackets = r'<!--\[[\s\S]+?\]-->'
roundbrackets = r'\([\s\S]+?\)'
latexbrackets = r'\\\cite\{[\s\S]+?\}'
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
no_user_response_count = 0

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

def fix_permissions(this_path):
    os.system("/bin/chmod 755 %s"%(this_path))

def check_dir(this_dir):
    if not os.path.isdir(this_dir):
        os.mkdir(this_dir)
    fix_permissions(this_dir)

def read_bib_files(localbibfile, globalbibfile=None):
    bfs = ""
    original_contents = ""
    bibfiles = [localbibfile]
    if globalbibfile:
        bibfiles.append(globalbibfile)
    for bibfile in bibfiles:
        if os.path.exists(bibfile):
            size = os.path.getsize(bibfile)
            if size > 0:
                try:
                    with open(bibfile, encoding="utf-8") as bf:
                        content = bf.read()
                except UnicodeDecodeError:
                    with open(bibfile, encoding="latin1") as bf:
                        content = bf.read()
            else:
                print("Bib file empty: {}".format(bibfile))
                content = ""
        else:
            print("File does not exist: {}".format(bibfile))
            continue  # Skip to the next file
        bfs += content
        if globalbibfile and bibfile == globalbibfile:
            original_contents += content
    try:
        bib_database = bibtexparser.bparser.BibTexParser(
            common_strings=True,
            homogenize_fields=True,
            interpolate_strings=False
        ).parse(bfs, partial=False)
        return bib_database, original_contents
    except Exception as e:
        print(f"Error parsing BibTeX files: {e}")
        return BibDatabase(), original_contents

def serialize_bib_database(bib_database):
    writer = BibTexWriter()
    writer.order_entries_by = None  # Keep the order of entries
    return writer.write(bib_database)

def hash_content(content):
    return hashlib.md5(content.encode('utf-8')).hexdigest()

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

def findreplace(inputtext, frdict):
    for f in frdict:
        inputtext = inputtext.replace(f, str(frdict[f]))
    return inputtext

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
            if priority_yaml[item] is None:
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
        h1 = re.findall(r'---[\s\S]+?---',filecontents)
        h2 = re.findall(r'---[\s\S]+?\.\.\.',filecontents)
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

def findcitation(info, infotype='pmid', additionalinfo=''):
    '''
        search the bibdat first
        search online in a variety of ways
        if found, use new to add it to global db
        return a single bibtex entry
        or None if not found or in doubt
    '''
    global no_user_response_count
    info = str(info.strip())
    if infotype == 'pmid':
        try:
            print ("searching for {} in bibdat".format(info))
            return pmids[info]
        except:
            print ("{} not found in pmids: ".format(info))
            pass
        pub = p2b([info])
        if len(pub) > 0:
            if pub[0] != 'null' and pub[0] != None:
                print ("PMID:{} found online".format(info))
                print (pub[0])
                new(pub[0])
                return pub[0]
        print ("PMID:{} NOT FOUND ON PUBMED".format(info))
        return None
    elif infotype == 'doi':
        msg=""
        try:
            return dois[info]
        except:
            msg += ("DOI not in bib file: {}".format(info))
        pub = search_pubmed(info, "doi")
        if len(pub) == 1:
            pubent = p2b(pub[0])
            if len(pubent) > 0:
                if pubent[0] != 'null' and pubent[0] != None:
                    print (msg + "but found in pubmed: {}".format(pubent[0]))
                    new(pubent[0])
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
                            new(pubent[0])
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
            full_bibdat.entries_dict[thisid]
        except:
            bestmatchingkey = find_similar_keys(thisid, full_bibdat.entries_dict)
            print(("biblatex id not found in biblatex file: {}. Best match in database is {}.".format(thisid, bestmatchingkey)))
            continue
        # if it is found, try to get the pmid
        if 'PMID' in full_bibdat.entries_dict[thisid]:
            pmidlist.append(full_bibdat.entries_dict[thisid]['PMID'])
        elif 'pmid' in full_bibdat.entries_dict[thisid]:
            pmidlist.append(full_bibdat.entries_dict[thisid]['pmid'])
        else:
            print ("pmid not found in bib file: {}. Searching online...".format(thisid))
            try:
                new_entry = findcitation(full_bibdat.entries_dict[thisid]['doi'], 'doi')
            except:
                notpmidlist.append(thisid) # hack
                new_entry = full_bibdat.entries_dict[thisid]
            if new_entry is None:
                new_entry = findcitation(full_bibdat.entries_dict[thisid]['title'], 'title')
                if new_entry is None:
                    notpmidlist.append(thisid)
                    continue
            try:
                pmidlist.append(new_entry['pmid'])
                pmids[new_entry['pmid']] = new_entry
            except:
                notpmidlist.append(thisid)
    return pmidlist, notpmidlist

def pmid2id(thesepmids, others):
    outids = []
    missing_ids = others
    for pmid in thesepmids:
        try:
            outids.append(pmids[pmid]['ID'])
            cite([pmids[pmid]['ID']])
        except:
            new_entry = findcitation(pmid, 'pmid')
            if new_entry is not None:
                outids.append(new_entry['ID'])
                new(new_entry)
            else:
                missing_ids.append(pmid)
    return outids, missing_ids

def format_inline(thisid):
    au = db.entries_dict[thisid]['Author'].split(" and ")
    if len(au)>1:
        au = au[0] + " et al"
    else:
        au = au[0]
    formatted_citation = "{}. {} {};{}:{}".format(
        au,
        db.entries_dict[thisid]['Journal'].capitalize(),
        db.entries_dict[thisid]['Year'],
        db.entries_dict[thisid]['Volume'],
        db.entries_dict[thisid]['Pages'],
        )
    return formatted_citation

def pmidout(pmidlist, notpmidlist):
    # make a blockstring
    blockstring = ''
    if len(pmidlist) > 0:
        # add to the outputdatabase if possible but since PMID is the output format, it doesn't matter if we can't find it
        for x in pmidlist:
            try:
                pmids[x]
            except:
                continue
            cite([pmids[x]['ID']])
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
    cite(theseids)
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

def search_pubmed(search_string, restrictfields=""):
    """
    Return a list of PMIDs for articles found by a search string.

    Args:
        search_string (str): The search term for querying PubMed.
        restrictfields (str): Optional field restriction for the search (e.g., 'title').

    Returns:
        list: A list of PMIDs for matching articles.
    """
    if not search_string or not search_string.strip():
        print("Empty search string provided.")
        return []

    max_returns = 500
    print(f"Searching PubMed for: {search_string}")

    # Modify search string for DOI searches
    if restrictfields.lower() == 'doi':
        search_string = f"{search_string}[AID]"
        restrictfields = ""  # Clear restrictfields

    # Build the query parameters
    params = {
        'db': 'pubmed',
        'term': search_string,
        'retmax': max_returns,
        'retmode': 'xml',
        'email': 'your_email@example.com',  # Replace with your actual email
        # 'api_key': 'your_api_key',        # Uncomment and set if you have an API key
    }
    # Do not include 'field' parameter if empty
    if restrictfields:
        params['field'] = restrictfields

    # Make the GET request to the ESearch API
    try:
        response = requests.get(
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
            params=params
        )
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"PubMed search failure: {e}")
        return []

    # Parse the XML response
    try:
        root = ET.fromstring(response.text)
    except ET.ParseError as e:
        print(f"Error parsing PubMed response: {e}")
        return []

    # Extract PMIDs from the response
    idlist = root.find('IdList')
    pmid_list = [id_elem.text for id_elem in idlist.findall('Id')] if idlist is not None else []

    foundcount = len(pmid_list)
    print(f"{foundcount} records found.")

    return pmid_list

def p2b(pmidlist):
    ''' by Nick Loman '''

    if type(pmidlist) != list:
        pmidlist = [str(pmidlist)]

    # Fetch XML data from Entrez.
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    url = '{}?db=pubmed&id={}&rettype=abstract'.format(efetch, ','.join(pmidlist))
    try:
        r = requests.get(url)
        r.raise_for_status()
    except requests.exceptions.RequestException:
        return []

    # Parse the XML using xml.etree.ElementTree
    bibout = []
    try:
        root = ET.fromstring(r.text)
    except ET.ParseError:
        return []

    for PubmedArticle in root.findall('.//PubmedArticle'):
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

        # Additional IDs
        PMCID = None
        DOI = None
        theseids = PubmedArticle.findall('./PubmedData/ArticleIdList/ArticleId')
        for thisid in theseids:
            if thisid.attrib.get('IdType') == 'pmc':
                PMCID = thisid
            elif thisid.attrib.get('IdType') == 'doi':
                DOI = thisid

        # Format author list
        authors = []
        for Author in PubmedArticle.findall('.//Author'):
            LastName = Author.find('LastName')
            ForeName = Author.find('ForeName')
            if LastName is not None and ForeName is not None:
                authors.append('{}, {}'.format(LastName.text, ForeName.text))

        # Use InvestigatorList if no authors found
        if not authors:
            for Investigator in PubmedArticle.findall('.//Investigator'):
                LastName = Investigator.find('LastName')
                ForeName = Investigator.find('ForeName')
                if LastName is not None and ForeName is not None:
                    authors.append('{}, {}'.format(LastName.text, ForeName.text))

        # Handle missing Year and Month
        if Year is None:
            MedlineDate = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate')
            if MedlineDate is not None and MedlineDate.text:
                Year = MedlineDate.text[:4]
                try:
                    month_abbr = MedlineDate.text[5:8]
                    Month = '{:02d}'.format(list(calendar.month_abbr).index(month_abbr))
                except ValueError:
                    Month = None
            else:
                Year = 'Unknown'
                Month = None
        else:
            Year = Year.text
            if Month is not None:
                Month = Month.text

        # Create a unique ID
        bib = {}
        authorname = authors[0].split(',')[0] if authors else ''
        try:
            titlewords = [x for x in ArticleTitle.text.split(' ') if len(x) > 3]
        except AttributeError:
            print("PUBMED ERROR - no article title for PMID:{}.".format(PMID.text if PMID else 'Unknown'))
            continue

        if len(titlewords) > 2:
            titlestring = ''.join(titlewords[:3])
        elif titlewords:
            titlestring = ''.join(titlewords)
        else:
            titlestring = ''

        if not (authorname + titlestring):
            titlestring = "PMID{}_".format(PMID.text if PMID else 'Unknown')

        new_id = '{}{}{}'.format(authorname, titlestring, Year).lower()
        new_id = re.sub(r'\W+', '', new_id)
        try:
            bib["ID"] = latexchars.replace_accents(new_id)
        except:
            bib["ID"] = new_id
        # Populate bibtex fields
        bib["Author"] = ' and '.join(authors)
        bib["Title"] = ArticleTitle.text if ArticleTitle is not None else ''
        bib["Journal"] = Title.text if Title is not None else ''
        bib["Year"] = Year
        if Volume is not None:
            bib["Volume"] = Volume.text
        if Issue is not None:
            bib["Number"] = Issue.text
        if MedlinePgn is not None:
            bib["Pages"] = MedlinePgn.text
        if Month is not None:
            bib["Month"] = Month
        if PMCID is not None:
            bib["pmcid"] = PMCID.text
        if DOI is not None:
            bib["doi"] = DOI.text
        if ISSN is not None:
            bib["ISSN"] = ISSN.text
        if PMID is not None:
            bib["pmid"] = PMID.text

        # Append to output list
        bibout.append(bib)

    return bibout


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

def main(
    filepath,
    bibfile,
    yaml_file,
    wholereference,
    outputstyle,
    overwrite,
    localbibonly,
    force_lowercase_citations,
    verbose
):
    # Determine source path and filenames
    sourcepath, filename = os.path.split(filepath)
    outpath = os.path.join(sourcepath)
    if outpath != '':
        check_dir(outpath)
    filestem = '.'.join(filename.split('.')[:-1])

    # Read input file
    with io.open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    # READ FROM INPUT FILE FIRST
    infileyaml = getyaml(text)
    workingyaml = copy.copy(infileyaml)
    # THEN READ FROM FILENAME.YAML
    yamlfile = os.path.join(sourcepath, filestem + ".yaml")
    if os.path.exists(yamlfile):
        with open(yamlfile) as f:
            workingyaml = mergeyaml(workingyaml, getyaml(f.read()))
    # THEN READ FROM SPECIFIED yaml_file
    if os.path.exists(yaml_file):
        with open(yaml_file) as f:
            workingyaml = mergeyaml(workingyaml, getyaml(f.read()))
    if verbose:
        print("Filepath:", filepath)
    input_file_extension = filepath.split('.')[-1]
    if filepath.endswith((".md", ".txt", ".qmd")):
        if outputstyle == 'null':
            outputstyle = 'md'
    elif filepath.endswith(".tex"):
        if outputstyle == 'null':
            outputstyle = 'tex'
    if overwrite:
        print("Overwriting original file")
        citelabel = "."
    elif outputstyle in ('pubmed', 'pmid'):
        print("Outputstyle = pmid")
        citelabel = ".citepmid."
    elif outputstyle in ('markdown', 'md'):
        print("Outputstyle = md")
        citelabel = ".citemd."
    elif outputstyle in ('latex', 'tex'):
        print("Outputstyle = tex")
        citelabel = ".citetex."
    elif outputstyle == 'inline':
        print("Outputstyle = inline")
        citelabel = "."
    else:
        citelabel = "."

    outputfile = os.path.join(outpath, filestem + citelabel + input_file_extension)

    text = make_unicode(text)
    text = readheader(text)[1]
    if verbose:
        print("Read input file:", filepath)



    # BIB - read them all and copy into one local version
    bibfile = os.path.abspath(os.path.expanduser(bibfile))
    if 'bibliography' in workingyaml.keys():
        print('Using YAML-specified bib:', workingyaml['bibliography'])
    else:
        workingyaml['bibliography'] = default_localbibname  # HARD OVERWRITE
    print("Using {} as bibout".format(workingyaml['bibliography']))
    
    # Prepare bibliography
    localbibpath = os.path.join(sourcepath, workingyaml['bibliography'])
    db, _ = read_bib_files(localbibpath)
    if localbibonly:
        print("\n*** Reading local bib file only: {} ***\n".format(localbibpath))
        full_bibdat, _ = read_bib_files(localbibpath)
    else:
        if verbose:
            print("Reading bibfiles:", bibfile, localbibpath)
        full_bibdat, original_bib_content = read_bib_files(localbibpath, bibfile)
    
    if force_lowercase_citations:
        print("Forcing lowercase citations")
        id_to_lower()
    make_alt_dicts()
    text = replace_blocks(text, outputstyle, use_whole=wholereference, flc=force_lowercase_citations)

    # Save bibliography for the current manuscript
    print('\nSaving bibliography for this file here:', localbibpath)
    outbib = bibtexparser.dumps(db)
    outbib = make_unicode(outbib)
    with open(localbibpath, "w", encoding="utf-8") as bf:
        bf.write(outbib)

    new_bib_content = serialize_bib_database(full_bibdat)

    # Save new global bibliography 
    if hash_content(original_bib_content) != hash_content(new_bib_content):
        print('\nSaving global bibliography here:', bibfile)
        bibfile_dir = os.path.dirname(bibfile)
        if not os.path.exists(bibfile_dir):
            os.makedirs(bibfile_dir, exist_ok=True)
        with open(bibfile, "w", encoding="utf-8") as bf:
            bf.write(new_bib_content)
        print("Global bibliography updated.")
    else:
        print("Global bibliography unchanged. Skipping write.")

    # Save new text file
    with io.open(outputfile, 'w', encoding='utf-8') as file:
        if outputstyle == "md":
            file.write('---\n{}\n---'.format(yaml.dump(workingyaml)).replace("\n\n", "\n"))
        file.write(text + "\n\n")
        print("Outputfile:", outputfile)

if __name__ == "__main__":
    config = getconfig()
    parser = argparse.ArgumentParser()
    # Essential arguments
    parser.add_argument("-f", '--filepath', help='Path to the input file', default="../genomicc-manuscript/manuscript.tex")
    # Additional files to specify
    parser.add_argument('-b', '--bibfile', default=default_global_bibfile, help='BibTeX file')
    parser.add_argument('-y', '--yaml', default='_quarto.yml', help='YAML file to use')
    # Other options
    parser.add_argument('-w', '--wholereference', action="store_true", default=False, help='Try to match whole references.')
    parser.add_argument('-o', '--outputstyle', type=str, choices=['md', '.qmd', 'markdown', 'tex', 'latex', 'pubmed', 'pmid', 'inline'], default='null', help='Output references format')
    parser.add_argument('-ow', '--overwrite', action="store_true", default=False, help='Overwrite input file with new version')
    parser.add_argument('-l', '--localbibonly', action="store_true", default=False, help='Use only local BibTeX file')
    parser.add_argument('-flc', '--force_lowercase_citations', action="store_true", default=False, help='Force all citation references into lowercase')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='Verbose output')
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(
        filepath=os.path.abspath(os.path.expanduser(args.filepath)),
        bibfile=args.bibfile,
        yaml_file=args.yaml,
        wholereference=args.wholereference,
        outputstyle=args.outputstyle,
        overwrite=args.overwrite,
        localbibonly=args.localbibonly,
        force_lowercase_citations=args.force_lowercase_citations,
        verbose=args.verbose
    )
