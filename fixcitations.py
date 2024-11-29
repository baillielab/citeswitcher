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
from collections import OrderedDict
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
default_global_bibfile = os.path.join(scriptpath,'_bibfiles/lib.pmid.bib')
#---

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
latexbrackets = r'\\cite\{[\s\S]+?\}'
pubmedsearchstrings = ["PMID:", "PubMed:"]
markdown_citations = r'\[@[\s\S]+?\]'
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

additionaldicts = []
def init():
    global full_bibdat
    full_bibdat = BibDatabase()
    global cited_bibdat # cited
    cited_bibdat = BibDatabase()
    global pmids
    pmids = {}
    global dois
    dois = {}
    additionaldicts.append((pmids, "pmid"))
    additionaldicts.append((dois, "doi"))
init()

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

def remove_duplicates_preserve_order(thislist):
    return list(OrderedDict.fromkeys(thislist))

def fix_permissions(this_path):
    os.system("/bin/chmod 755 %s"%(this_path))

def check_dir(this_dir):
    if not os.path.isdir(this_dir):
        os.mkdir(this_dir)
    fix_permissions(this_dir)

def id_to_lower(bibdat):
    for entry in bibdat.entries:
        entry["ID"] = entry["ID"].lower()
    #bibdat._entries_dict = None
    bibdat.entries = bibdat.entries  # Reassign entries to itself to force cache update
    return bibdat

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

def add_entry_to_bibdatabase(entry, bibdat, additional_dicts=None):
    '''
    Adds a new entry to the given BibDatabase object, checking for duplicates based on ID, DOI, and PMID.
    Merges entries if a duplicate is found, giving precedence to existing data.

    Parameters:
    - entry: The new BibTeX entry as a dictionary.
    - bibdat: The BibDatabase object to which the entry should be added.
    - additional_dicts: A list of tuples containing dictionaries for additional identifiers (e.g., DOI, PMID)
      and their corresponding field names in entries. Example: [(dois_dict, 'doi'), (pmids_dict, 'pmid')]

    Returns:
    - None. The function modifies bibdat in place.
    '''

    if entry is None:
        return  # Nothing to add
    if 'ENTRYTYPE' not in entry: # Ensure the entry has an ENTRYTYPE
        entry['ENTRYTYPE'] = 'article'
    existing_entry = None
    if entry['ID'] in bibdat.entries_dict: # First, try to find an existing entry by ID
        existing_entry = bibdat.entries_dict[entry['ID']]
    else: # If not found by ID, try to find by DOI or PMID
        if additional_dicts:
            for id_dict, id_field in additional_dicts:
                if id_field in entry:
                    id_value = entry[id_field]
                    if id_value in id_dict:
                        existing_entry = id_dict[id_value]
                        break  # Found existing entry by DOI or PMID
    if existing_entry:
        #existing_entry = merge_entries(existing_entry, entry)
        merged_entry = merge_entries(existing_entry, entry)
        index = bibdat.entries.index(existing_entry)
        bibdat.entries[index] = merged_entry
        if args.verbose:
            print("Merged two entries in a bib database under this ID:\n {}".format(existing_entry['ID']))
    else:
        #entry = cleanbib(entry) # commented out because \_ underscore escape is annoying
        bibdat.entries.append(entry)
    bibdat.entries = bibdat.entries  # Reassign entries to itself to force cache update
    if additional_dicts:
        for id_dict, id_field in additional_dicts:
            if id_field in entry:
                id_value = entry[id_field]
                id_dict[id_value] = entry

def cite(theseids):
    """
    Adds the specified citation IDs to the local bibliography database (cited_bibdat).
    Ensures that entries in cited_bibdat are identical to those in full_bibdat (global bibliography).
    If an entry is not found in full_bibdat, it attempts to find and add it via online search.
    """
    global full_bibdat, cited_bibdat
    if args.verbose:
        print("cite function has been asked to handle:\n", theseids)
    fails = []
    for thisid in theseids:
        # Ensure the citation ID is in lowercase if force_lowercase_citations is True
        if args.force_lowercase_citations:
            thisid = thisid.lower()
        # Check if the entry is already in cited_bibdat
        if thisid in cited_bibdat.entries_dict:
            if args.verbose:
                print("\t\tfound", thisid, "in cited_bibdat")
            continue  # Entry is already in local bibliography
        # Check if the entry exists in full_bibdat
        try:
            full_bibdat.entries_dict[thisid]
            if args.verbose:
                print("\t\tfound", thisid, "in full_bibdat")
        except:
            new_entry = findcitation(thisid, 'id')
            if new_entry:
                add_entry_to_bibdatabase(new_entry, full_bibdat, additional_dicts=additionaldicts)
                add_entry_to_bibdatabase(new_entry, cited_bibdat, additional_dicts=None)
            else:
                print("Cite function unable to find entry for ID:", thisid)
                fails.append(thisid)
            continue
        add_entry_to_bibdatabase(full_bibdat.entries_dict[thisid], cited_bibdat, additional_dicts=None)
    return fails

def cleanbib(bibtex_entry):
    def safe_string_to_latex(value):
        if value is None:
            return ''  # or return None, depending on your needs
        return string_to_latex(str(value))  # Convert to string first
    
    return {d: safe_string_to_latex(bibtex_entry[d]) for d in bibtex_entry.keys()}

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

def merge_entries(entry1, entry2):
    """
    Merge two BibTeX entries, giving precedence to the first entry's data.
    Returns the merged entry.
    entry1 takes precedence
    """
    merged_entry = entry1.copy()
    for key, value in entry2.items():
        if key not in merged_entry or not merged_entry[key]:
            merged_entry[key] = value
    return merged_entry

def merge_bibdat_duplicates(bib_database1, bib_database2=BibDatabase()):
    entries_ordered = OrderedDict()
    id_changes = {}
    for entry in bib_database1.entries + bib_database2.entries:
        entry_id = entry['ID']
        doi = entry.get('doi')
        pmid = entry.get('pmid')
        existing_entry = entries_ordered.get(entry_id)
        if existing_entry: # Duplicate by ID
            merged_entry = merge_entries(existing_entry, entry) # existing_entry takes precedence
            entries_ordered[entry_id] = merged_entry
        else: # Check for duplicate by DOI or PMID
            duplicate_found = False
            for e_id, e in entries_ordered.items():
                if (doi and e.get('doi') == doi) or (pmid and e.get('pmid') == pmid):
                    print (f"Duplicate found: {entry_id} will be merged into {e_id}")
                    merged_entry = merge_entries(e, entry)
                    entries_ordered[e_id] = merged_entry
                    if entry_id != e_id:
                        id_changes[entry_id] = e_id
                    duplicate_found = True
                    break
            if not duplicate_found:
                entries_ordered[entry_id] = entry
    print (entries_ordered['eqtlgen2021'])
    merged_bib_database = BibDatabase()
    for entry in list(entries_ordered.values()):
        add_entry_to_bibdatabase(entry, merged_bib_database)
    print (merged_bib_database.entries_dict['eqtlgen2021'])
    return merged_bib_database, id_changes

def parse_bib_contents(bibfilecontents):
    try:
        bib_database = bibtexparser.bparser.BibTexParser(
            common_strings=True,
            homogenize_fields=True,
            interpolate_strings=False
        ).parse(bibfilecontents, partial=False)
        return bib_database
    except Exception as e:
        print(f"Error parsing parse_bib_contents: {e}")
        return BibDatabase()

def read_bib_files(localbibfile, globalbibfile=None):
    original_contents = ""
    # Read and parse the local bib file
    parser = BibTexParser(common_strings=True)
    local_bib_database = BibDatabase()
    if os.path.exists(localbibfile):
        size = os.path.getsize(localbibfile)
        if size > 0:
            try:
                with open(localbibfile, encoding="utf-8") as bf:
                    content = bf.read()
            except UnicodeDecodeError:
                with open(localbibfile, encoding="latin1") as bf:
                    content = bf.read()
            local_bib_database = parse_bib_contents(content)
            try:
                local_bib_database = parse_bib_contents(content)
            except Exception as e:
                print(f"Error parsing BibTeX file {localbibfile}: {e}")
        else:
            print("Bib file empty: {}".format(localbibfile))
            original_contents = ""
    else:
        print("File does not exist: {}".format(localbibfile))
        original_contents = ""
    
    print (local_bib_database.entries_dict['eqtlgen2021'])
    local_bib_database, id_changes_local = merge_bibdat_duplicates(local_bib_database)
    local_bib_database.entries = local_bib_database.entries
    print (local_bib_database.entries_dict['eqtlgen2021'])

    # Read and parse the global bib file
    global_bib_database = BibDatabase()
    if globalbibfile and os.path.exists(globalbibfile):
        size = os.path.getsize(globalbibfile)
        if size > 0:
            try:
                with open(globalbibfile, encoding="utf-8") as bf:
                    content = bf.read()
            except UnicodeDecodeError:
                with open(globalbibfile, encoding="latin1") as bf:
                    content = bf.read()
            original_contents = content
            try:
                global_bib_database = parse_bib_contents(content)
            except Exception as e:
                print(f"Error parsing BibTeX file {globalbibfile}: {e}")
        else:
            print("Bib file empty: {}".format(globalbibfile))
            original_contents = ""
        combined_bib_database, id_changes_global = merge_bibdat_duplicates(local_bib_database, global_bib_database)
    elif globalbibfile:
        print("File does not exist: {}".format(globalbibfile))
        original_contents = ""
    else:
        combined_bib_database = local_bib_database
        id_changes_global = {}
    print ("\n")

    print (combined_bib_database.entries_dict['eqtlgen2021'])

    id_changes = {**id_changes_local, **id_changes_global}
    return combined_bib_database, original_contents, id_changes

def serialize_bib_database(bib_database):
    writer = BibTexWriter()
    writer.order_entries_by = None  # Keep the order of entries
    return writer.write(bib_database)

def hash_content(content):
    return hashlib.md5(content.encode('utf-8')).hexdigest()

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

def make_unicode(inputstring):
    if type(inputstring) != str:
        inputstring =  inputstring.decode('utf-8')
    nbspace = re.compile(u"\N{NO-BREAK SPACE}", re.IGNORECASE)
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
        nested_out = [x for x in nested_out if not re.match(thislabel, x[1:])]
        if args.verbose:
            for x in nested_out:
                if re.match(thislabel, x[1:]):
                    if not x[1:].startswith(thislabel):
                        print ("regex picked up this match: {} {} which was missed by startswith".format(x[1:], thislabel))
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

def get_wholereference_citation_blocks(inputtext, force_lowercase=False):
    '''
        finds wholecitations AND doi
        return a list of citation blocks, and new input text with these removed
    '''
    confirmed_blocks = []
    for theseparetheses in [squarebrackets, curlybrackets]:
        for b in get_parenthesised(inputtext, [theseparetheses]):
            b = b.replace('\n',' ')
            if "." in b or ":" in b:
                new_entry = parse_wholecitation_block(b, force_lowercase)
                if new_entry is not None:
                    confirmed_blocks.append(str(b))
                    inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks, inputtext

def parse_citation_block(thisblock, force_lowercase=False):
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
    if force_lowercase:
        results['ids'] = [x.lower() for x in results['ids']]
    return results

def clean_id(thisid):

    def casereplace(thistext, catchall, cleanversion):
        '''
        replace <catchall> text with <cleanversion> if it occurs at the beginning of the text
        '''
        catchall = re.compile("^"+re.escape(catchall.strip()), re.IGNORECASE)
        return catchall.sub(cleanversion, thistext)

    thisid = thisid.replace(" ","").strip()
    for stem in pubmedsearchstrings:
        thisid = casereplace(thisid, stem, 'PMID:')
    for stem in mdsearchstrings:
        thisid = casereplace(thisid, stem, '@')
    for stem in doisearchstrings:
        thisid = casereplace(thisid, stem, 'DOI:')
    return thisid

askedalready = {}
def parse_wholecitation_block(thisblock, force_lowercase=False):
    ''' take a block of text, and return two lists of ids '''
    try:
        return askedalready[thisblock]
    except:
        outids = []
        # try a doi search first in case there's a doi in there
        d = parse_citation_block(thisblock, force_lowercase)
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
        if force_lowercase:
            outids = [x.lower() for x in outids]
        return outids, []  # there is only ever one

def findcitation(info, infotype='pmid', additionalinfo='', force_search=False):
    '''
        search the bibdat first
        search online in a variety of ways
        if found, add it to global db
        return a single bibtex entry
        or None if not found or in doubt
    '''
    global no_user_response_count
    info = str(info.strip())
    if infotype == 'pmid':
        if not force_search:
            try:
                return pmids[info]
            except:
                print ("{} not found in pmids: ".format(info))
                pass
        pub = p2b([info])
        if len(pub) > 0:
            if pub[0] != 'null' and pub[0] != None:
                print ("PMID:{} found online".format(info))
                print (pub[0])
                add_entry_to_bibdatabase(pub[0], full_bibdat, additional_dicts=additionaldicts)
                return pub[0]
        print ("PMID:{} NOT FOUND ON PUBMED".format(info))
        return None
    elif infotype == 'doi':
        msg=""
        if not force_search:
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
                    add_entry_to_bibdatabase(pubent[0], full_bibdat, additional_dicts=additionaldicts)
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
                            add_entry_to_bibdatabase(pubent[0], full_bibdat, additional_dicts=additionaldicts)
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
            print(("id2pmid: bibtex id not found in bibtex files: {}. Best match in database is {}.".format(thisid, bestmatchingkey)))
            continue
        # if it is found, try to get the pmid
        if 'PMID' in full_bibdat.entries_dict[thisid]:
            pmidlist.append(full_bibdat.entries_dict[thisid]['PMID'])
        elif 'pmid' in full_bibdat.entries_dict[thisid]:
            pmidlist.append(full_bibdat.entries_dict[thisid]['pmid'])
        else:
            print ("PMID not found in bib file: {}. Searching online...".format(thisid))
            if 'doi' in full_bibdat.entries_dict[thisid]:
                if args.verbose: 
                    print (f"searching pubmed for DOI: {full_bibdat.entries_dict[thisid]['doi']}")
                new_entry = findcitation(full_bibdat.entries_dict[thisid]['doi'], 'doi', '', force_search=True)
            elif 'title' in full_bibdat.entries_dict[thisid]:
                if args.verbose: 
                    print (f"searching pubmed for Title: {full_bibdat.entries_dict[thisid]['title']}")
                new_entry = findcitation(full_bibdat.entries_dict[thisid]['title'], 'title', additionalinfo='', force_search=True)
            if new_entry is None:
                notpmidlist.append(thisid)
                continue
            if 'pmid' in new_entry:
                pmidlist.append(new_entry['pmid'])
                pmids[new_entry['pmid']] = new_entry
            else:
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
                add_entry_to_bibdatabase(new_entry, full_bibdat, additional_dicts=additionaldicts)
            else:
                missing_ids.append(pmid)
    return outids, missing_ids

def format_inline(thisid):
    au = cited_bibdat.entries_dict[thisid]['Author'].split(" and ")
    if len(au)>1:
        au = au[0] + " et al"
    else:
        au = au[0]
    formatted_citation = "{}. {} {};{}:{}".format(
        au,
        cited_bibdat.entries_dict[thisid]['Journal'].capitalize(),
        cited_bibdat.entries_dict[thisid]['Year'],
        cited_bibdat.entries_dict[thisid]['Volume'],
        cited_bibdat.entries_dict[thisid]['Pages'],
        )
    return formatted_citation

def pmidout(pmidlist, notpmidlist):
    '''
    make a blockstring to replace PMIDs in a citation list
    '''
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
    '''
    make a blockstring to replace IDs in a citation list
    '''
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

def replace_ids_in_text(text, id_changes):
    """
    Replace occurrences of old IDs with new IDs in the text, using the global regex patterns.
    Returns both the updated text and a report of changes made.
    
    Returns:
        tuple: (updated_text, report) where report is a string describing all changes made
    """
    patterns = [
        markdown_citations,
        latexbrackets
    ]
    combined_pattern = '(' + '|'.join(patterns) + ')'
    report_lines = []    
    def replace_ids_in_match(match):
        citation_block = match.group(0)
        original_block = citation_block  # Keep for debugging
        for old_id, new_id in id_changes.items():
            if old_id != new_id:  # Use word boundaries to match whole IDs only
                citation_block = re.sub(r'\b' + re.escape(old_id) + r'\b', new_id, citation_block)
        if citation_block != original_block:
            newline = "{:>80} ==> {}".format(original_block, citation_block)
            if newline not in report_lines:
                report_lines.append(newline)
        return citation_block
    
    updated_text = re.sub(combined_pattern, replace_ids_in_match, text)
    report = "\n".join(report_lines) if report_lines else "No replacements made"
    
    return updated_text, report

def replace_blocks(thistext, outputstyle="md", use_whole=False, flc=False):
    workingtext = thistext
    p, workingtext = get_mixed_citation_blocks(workingtext) # pmid first as they are the most likely to have errors
    l, workingtext = [x for x in get_latex_citation_blocks(workingtext) if x not in p]
    w, workingtext = [x for x in get_word_citation_blocks(workingtext) if x not in p+l]
    if use_whole:
        r, workingtext = [x for x in get_wholereference_citation_blocks(workingtext, flc) if x not in p+l+w]
    else:
        r=[]
    print ("\nNumber blocks using pmid or doi or md:{}".format(len(p)))
    print ("Number blocks using latex:{}".format(len(l)))
    print ("Number blocks using wholeref:{}".format(len(r)))
    replacedict = {}
    # there are slightly different procedures for each reference type:
    for b in p+l+w:
        citedhere = parse_citation_block(b, flc)
        if outputstyle == 'md' or outputstyle=='tex' or outputstyle=='inline':
            theseids, notfound = pmid2id(citedhere['pmids'], citedhere['notfound']) # ids added to bib
            theseids = remove_duplicates_preserve_order(theseids+citedhere['ids'])
            replacedict[b] = mdout(theseids, notfound, outputstyle, flc=flc)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(citedhere['ids'], citedhere['notfound']) # ids added to bib
            pm = remove_duplicates_preserve_order(pm+citedhere['pmids'])
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in r:
        theseids, theseothers = parse_wholecitation_block(b, flc)
        if outputstyle == 'md' or outputstyle=='tex' or outputstyle=='inline':
            replacedict[b] = mdout(theseids, theseothers, outputstyle, flc=flc)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(theseids) # ids added to bib
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    if len(replacedict)>0:
        print ("\n Citations replaced in text:")
    for b in replacedict:
        if replacedict[b] == 'null' or b==replacedict[b]:
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

def print_all_ids(bibdat):
    print ("=======\nPrinting all ids:\n")
    for entry in bibdat.entries:
        print(entry['ID'])

def main(
    filepath,
    globalbibfile,
    yaml_file,
    wholereference,
    outputstyle,
    safemode,
    localbibonly,
    force_lowercase_citations,
    verbose
):
    global full_bibdat, cited_bibdat
    # Determine source path and filenames
    sourcepath, filename = os.path.split(filepath)
    outpath = os.path.join(sourcepath)
    if outpath != '':
        check_dir(outpath)
    filestem = '.'.join(filename.split('.')[:-1])

    # Read input file
    with io.open(filepath, "r", encoding="utf-8") as f:
        text = f.read()
    text = make_unicode(text)
    if verbose:
        print("Read input file:", filepath)

    # READ YAML FROM INPUT FILE FIRST
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
    
    # GET OUTPUTSTYLE
    input_file_extension = filepath.split('.')[-1]
    if filepath.endswith((".md", ".txt", ".qmd")):
        if outputstyle == 'null':
            outputstyle = 'md'
    elif filepath.endswith(".tex"):
        if outputstyle == 'null':
            outputstyle = 'tex'
    if safemode:
        if outputstyle in ('pubmed', 'pmid'):
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
        print("Overwriting original file")
        citelabel = "."

    # BIB - read them all and copy into one local version
    if 'bibliography' in workingyaml.keys():
        print('Using YAML-specified bib:', workingyaml['bibliography'])
    else:
        workingyaml['bibliography'] = default_localbibname  # HARD OVERWRITE
    localbibpath = os.path.join(sourcepath, workingyaml['bibliography'])
    print(f"Using {localbibpath} as bibout")
    original_bib_content = None

    # Read bib files and get ID changes
    if localbibonly:
        full_bibdat, original_bib_content, id_changes = read_bib_files(localbibpath)
    else:
        globalbibfile = os.path.abspath(os.path.expanduser(globalbibfile))
        full_bibdat, original_bib_content, id_changes = read_bib_files(localbibpath, globalbibfile)
    full_bibdat.entries = full_bibdat.entries  # Force cache update

    # Force lowercase citations if required
    if force_lowercase_citations:
        print("Forcing lowercase citations")
        full_bibdat = id_to_lower(full_bibdat)
        cited_bibdat = id_to_lower(cited_bibdat)

    # Rebuild additional dictionaries
    make_alt_dicts()

    # Replace IDs in text based on id_changes
    if id_changes:
        text, id_change_report = replace_ids_in_text(text, id_changes)

    if len(id_change_report)>0:
        print ("\n ID changes in the following citations:")
        print (id_change_report)

    text = readheader(text)[1]
    text = replace_blocks(text, outputstyle, use_whole=wholereference, flc=force_lowercase_citations)

    if args.verbose:
        print(f"Number of entries before duplicate removal (cited_bibdat): {len(cited_bibdat.entries)}")
    cited_bibdat, id_changes_cited = merge_bibdat_duplicates(cited_bibdat)
    if args.verbose:
        print(f"Number of entries after duplicate removal (cited_bibdat): {len(cited_bibdat.entries)}")
    
    if args.verbose:
        print(f"Number of entries before duplicate removal (full_bibdat): {len(full_bibdat.entries)}")
    full_bibdat, id_changes_full = merge_bibdat_duplicates(full_bibdat)
    if args.verbose:
        print(f"Number of entries after duplicate removal (full_bibdat): {len(full_bibdat.entries)}")
    
    if len(id_changes_cited)>0:
        for item in id_changes_cited:
            print ("Additional id change: {}  ==> {}".format(item, id_changes_cited[item]))
    if args.verbose:
        if len(id_changes_full)>0:
            for item in id_changes_full:
                print ("Additional id change: {}  ==> {}".format(item, id_changes_full[item]))

    bibdir, bibfilename = os.path.split(localbibpath)
    bibstem = '.'.join(bibfilename.split('.')[:-1])
    localbibpath = os.path.join(bibdir, bibstem + citelabel + "bib")
    print('\nSaving bibliography for this file here:', localbibpath)
    outbib = bibtexparser.dumps(cited_bibdat)
    outbib = make_unicode(outbib)
    with open(localbibpath, "w", encoding="utf-8") as bf:
        bf.write(outbib)

    new_bib_content = serialize_bib_database(full_bibdat)
    # Save new global bibliography 
    if original_bib_content is not None:
        if hash_content(original_bib_content) != hash_content(new_bib_content):
            if not args.safemode:
                print('\nSaving updated global bibliography here:', globalbibfile)
                os.makedirs(os.path.dirname(globalbibfile), exist_ok=True)
                with open(globalbibfile, "w", encoding="utf-8") as bf:
                    bf.write(new_bib_content)
        else:
            print("Global bibliography unchanged. Skipping write.")

    # Save new text file
    outputfile = os.path.join(outpath, filestem + citelabel + input_file_extension)
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
    parser.add_argument('-gb', '--globalbibfile', default=default_global_bibfile, help='BibTeX file')
    parser.add_argument('-y', '--yaml', default='_quarto.yml', help='YAML file to use')
    # Other options
    parser.add_argument('-w', '--wholereference', action="store_true", default=False, help='Try to match whole references.')
    parser.add_argument('-o', '--outputstyle', type=str, choices=['md', '.qmd', 'markdown', 'tex', 'latex', 'pubmed', 'pmid', 'inline'], default='null', help='Output references format')
    parser.add_argument('-s', '--safemode', action="store_true", default=False, help='Overwrite input file with new version')
    parser.add_argument('-l', '--localbibonly', action="store_true", default=False, help='Use only local BibTeX file')
    parser.add_argument('-flc', '--force_lowercase_citations', action="store_true", default=False, help='Force all citation references into lowercase')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='Verbose output')
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(
        filepath=os.path.abspath(os.path.expanduser(args.filepath)),
        globalbibfile=args.globalbibfile,
        yaml_file=args.yaml,
        wholereference=args.wholereference,
        outputstyle=args.outputstyle,
        safemode=args.safemode,
        localbibonly=args.localbibonly,
        force_lowercase_citations=args.force_lowercase_citations,
        verbose=args.verbose
    )


