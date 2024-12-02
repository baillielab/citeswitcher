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

class BibData:
    def __init__(self):
        self.entries = []
        self.entries_dict = {}   # Keyed by ID for quick access
        self.pmids_dict = {}     # Keyed by PMID for quick access
        self.dois_dict = {}      # Keyed by DOI for quick access
        self.id_changes = {}     # Maps old IDs to new IDs
        self.comments = []       # to prevent bibtexparser warnings
        self.preambles = []      # to prevent bibtexparser warnings
        self.strings = {}        # to prevent bibtexparser warnings

    def add_entry(self, entry):
        # Ensure the entry has an 'ENTRYTYPE', default to 'article'
        if 'ENTRYTYPE' not in entry or not entry['ENTRYTYPE']:
            entry['ENTRYTYPE'] = 'article'

        entry_id = entry['ID']
        # Resolve entry_id if it has been changed
        resolved_id = self.id_changes.get(entry_id, entry_id)
        existing_entry = self.entries_dict.get(resolved_id)

        if existing_entry:
            # Merge entries if duplicate IDs are found
            merged_entry, changed_ids = self.merge_entries(existing_entry, entry)
            index = self.entries.index(existing_entry)
            self.entries[index] = merged_entry
            self.entries_dict[merged_entry['ID']] = merged_entry

            # Update id_changes with any new ID changes
            self.id_changes.update(changed_ids)
        else:
            # If the entry ID has changed, record the change
            if resolved_id != entry_id:
                self.id_changes[entry_id] = resolved_id
                entry['ID'] = resolved_id

            self.entries.append(entry)
            self.entries_dict[resolved_id] = entry

        # Index the entry by PMID
        pmid = entry.get('pmid') or entry.get('PMID')
        if pmid:
            self.pmids_dict[pmid] = entry

        # Index the entry by DOI
        doi = entry.get('doi')
        if doi:
            self.dois_dict[doi] = entry

    def merge_entries(self, entry1, entry2):
        """
        Merge two entries, giving precedence to entry1.
        Records any changes in IDs and returns them.
        """
        merged_entry = entry1.copy()
        changed_ids = {}

        for key, value in entry2.items():
            if key not in merged_entry or not merged_entry[key]:
                merged_entry[key] = value

        # If IDs are different, map old ID to new ID
        id1 = entry1['ID']
        id2 = entry2['ID']
        if id1 != id2:
            # Prefer the ID from entry1 (existing entry)
            self.id_changes[id2] = id1
            changed_ids[id2] = id1
            # Update the entries_dict to reflect the ID change
            if id2 in self.entries_dict:
                del self.entries_dict[id2]
            self.entries_dict[id1] = merged_entry

        return merged_entry, changed_ids

    def id_to_lower(self):
        new_entries = []
        new_entries_dict = {}
        new_pmids_dict = {}
        new_dois_dict = {}
        new_id_changes = {}

        for entry in self.entries:
            old_id = entry['ID']
            new_id = old_id.lower()

            if new_id != old_id:
                entry['ID'] = new_id
                # Record the ID change
                self.id_changes[old_id] = new_id
                new_id_changes[old_id] = new_id

            # Ensure 'ENTRYTYPE' is present
            if 'ENTRYTYPE' not in entry or not entry['ENTRYTYPE']:
                entry['ENTRYTYPE'] = 'article'

            new_entries.append(entry)
            new_entries_dict[new_id] = entry

            # Update PMID index
            pmid = entry.get('pmid') or entry.get('PMID')
            if pmid:
                new_pmids_dict[pmid] = entry

            # Update DOI index
            doi = entry.get('doi')
            if doi:
                new_dois_dict[doi] = entry

        self.entries = new_entries
        self.entries_dict = new_entries_dict
        self.pmids_dict = new_pmids_dict
        self.dois_dict = new_dois_dict
        self.id_changes.update(new_id_changes)

    def get_entry_by_id(self, entry_id):
        # Resolve entry_id if it has been changed
        resolved_id = self.id_changes.get(entry_id, entry_id)
        return self.entries_dict.get(resolved_id)

    def get_entry_by_pmid(self, pmid):
        return self.pmids_dict.get(pmid)

    def get_entry_by_doi(self, doi):
        return self.dois_dict.get(doi)

    def get_entries(self):
        return self.entries

    def get_id_changes(self):
        return self.id_changes


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
        if args.force_lowercase_citations:
            thisid = thisid.lower()
        if cited_bibdat.get_entry_by_id(thisid):
            if args.verbose:
                print("\t\tfound", thisid, "in cited_bibdat")
            continue  # Entry is already in local bibliography
        entry = full_bibdat.get_entry_by_id(thisid)
        if entry:
            if args.verbose:
                print("\t\tfound", thisid, "in full_bibdat")
            cited_bibdat.add_entry(entry)
        else:
            new_entry = findcitation(thisid, 'id')
            if new_entry:
                full_bibdat.add_entry(new_entry)
                cited_bibdat.add_entry(new_entry)
            else:
                print("Cite function unable to find entry for ID:", thisid)
                fails.append(thisid)
    return fails

def merge_bibdat_duplicates(bib_data1, bib_data2=None):
    """
    Merges two BibData instances, handling duplicates based on IDs, DOIs, and PMIDs.
    Modifies bib_data1 in place and records id_changes in bibdata1.get_id_changes()

    Parameters:
    - bib_data1: The primary BibData instance.  This one takes precedence.
    - bib_data2: The secondary BibData instance to merge into bib_data1.

    Returns:
    - merged_bib_data: A new BibData instance containing merged entries.
    - id_changes: A dictionary mapping old IDs: new IDs if changes occurred.
    """

    entries_ordered = OrderedDict()
    if bib_data2 is None:
        bib_data2 = BibData()
    print (f"Merging two bib_dats: {len(bib_data1.entries)} and {len(bib_data2.entries)} entries")

    # Combine entries from both BibData instances
    combined_entries = bib_data1.get_entries() + bib_data2.get_entries()

    for entry in combined_entries:
        entry_id = entry['ID']
        doi = entry.get('doi')
        pmid = entry.get('pmid') or entry.get('PMID')
        existing_entry = entries_ordered.get(entry_id)

        if existing_entry:
            # Duplicate found by ID
            bib_data1.merge_entries(existing_entry, entry)
            entries_ordered[entry_id] = bib_data1.get_entry_by_id(entry_id)
        else:
            # Check for duplicate by DOI or PMID
            duplicate_found = False
            for e_id, e in entries_ordered.items():
                e_doi = e.get('doi')
                e_pmid = e.get('pmid') or e.get('PMID')
                if (doi and e_doi == doi) or (pmid and e_pmid == pmid):
                    print(f"Duplicate found: {entry_id} will be merged into {e_id} len(entries_ordered) = {len(entries_ordered)}")
                    bib_data1.merge_entries(e, entry)
                    entries_ordered[e_id] = bib_data1.get_entry_by_id(e_id)
                    duplicate_found = True
                    break
            if not duplicate_found:
                entries_ordered[entry_id] = entry
                bib_data1.add_entry(entry) # add all entries to bib_data1
    print (f"Merged length = {len(bib_data1.entries)} entries")

def parse_bib_contents(bibfilecontents):
    try:
        parser = BibTexParser(common_strings=True)
        bib_database = parser.parse(bibfilecontents)
        bib_data = BibData()
        for entry in bib_database.entries:
            bib_data.add_entry(entry)
        return bib_data
    except Exception as e:
        print(f"Error parsing parse_bib_contents: {e}")
        return BibData()

def read_bib_file(bibfilepath):
    original_contents = ""
    bibdata = BibData()
    if os.path.exists(bibfilepath):
        size = os.path.getsize(bibfilepath)
        if size > 0:
            try:
                with open(bibfilepath, encoding="utf-8") as bf:
                    content = bf.read()
            except UnicodeDecodeError:
                with open(bibfilepath, encoding="latin1") as bf:
                    content = bf.read()
            original_contents = copy.copy(content)
            try:
                bibdata = parse_bib_contents(content)
            except Exception as e:
                print(f"Error parsing BibTeX file {bibfilepath}: {e}")
        else:
            print("Bib file empty: {}".format(bibfilepath))
            original_contents = ""
    else:
        print("File does not exist: {}".format(bibfilepath))
        original_contents = ""
    # Force lowercase citations if required
    if args.force_lowercase_citations:
        print("Forcing lowercase citations")
        bibdata.id_to_lower()
    return bibdata, original_contents

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
    """
    Searches for a citation in full_bibdat or online.
    """
    global no_user_response_count, full_bibdat
    info = info.strip()
    if infotype == 'pmid':
        if not force_search:
            entry = full_bibdat.get_entry_by_pmid(info)
            if entry:
                return entry
            else:
                print(f"PMID {info} not found in full_bibdat.")
        # Search online
        pub_entries = p2b([info])
        if pub_entries:
            new_entry = pub_entries[0]
            if new_entry:
                full_bibdat.add_entry(new_entry)
                return new_entry
        print(f"PMID {info} not found online.")
        return None
    elif infotype == 'doi':
        if not force_search:
            entry = full_bibdat.get_entry_by_doi(info)
            if entry:
                return entry
            else:
                print(f"DOI {info} not found in full_bibdat.")
        # Search online
        pmids = search_pubmed(info, "doi")
        if pmids:
            pub_entries = p2b(pmids)
            if pub_entries:
                new_entry = pub_entries[0]
                if new_entry:
                    full_bibdat.add_entry(new_entry)
                    return new_entry
        print(f"DOI {info} not found online.")
        return None
    elif infotype == 'title':
        if no_user_response_count > 2:
            print("No user response to previous 3 queries. Running in silent mode.")
            return None
        print(f"Searching PubMed for title: {info}")
        pmids = search_pubmed(info, "title")
        if len(pmids) == 1:
            pmid = pmids[0]
            pub_entries = p2b([pmid])
            if pub_entries:
                new_entry = pub_entries[0]
                question = "--------------\n\
New citation (PMID:{}) found in Pubmed. Please check that input is the same as the found citation for this reference block: \n\
\n{}\n\n{}\n{}\n\n{}\n\n\
Enter y/n within 10 seconds".format(
                    pmid,
                    additionalinfo,
                    "{:>12}:    {}".format('Input Title', info),
                    "{:>12}:    {}".format('Found Title', pub_entries[0]['Title']),
                    '\n'.join( ["{:>12}:    {}".format(x,pub_entries[0][x]) for x in pub_entries[0] if x != "Title"])
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
                        full_bibdat.add_entry(new_entry)
                        return new_entry
        print(f"Title {info} not found online.")

def id2pmid(theseids, notpmidlist=[]):
    """
    Converts a list of IDs to PMIDs.
    """
    global full_bibdat
    pmidlist = []
    for thisid in theseids:
        # Try to find this ID in full_bibdat
        entry = full_bibdat.get_entry_by_id(thisid)
        if entry:
            pmid = entry.get('pmid') or entry.get('PMID')
            if pmid:
                pmidlist.append(pmid)
            else:
                print(f"PMID not found for ID {thisid}. Searching online...")
                # Attempt to find the entry online using DOI or title
                new_entry = None
                doi = entry.get('doi')
                title = entry.get('title') or entry.get('Title')
                if doi:
                    new_entry = findcitation(doi, 'doi', force_search=True)
                elif title:
                    new_entry = findcitation(title, 'title', force_search=True)
                if new_entry and (new_entry.get('pmid') or new_entry.get('PMID')):
                    pmid = new_entry.get('pmid') or new_entry.get('PMID')
                    pmidlist.append(pmid)
                    full_bibdat.add_entry(new_entry)
                else:
                    notpmidlist.append(thisid)
        else:
            print(f"ID {thisid} not found in full_bibdat.")
            notpmidlist.append(thisid)
    return pmidlist, notpmidlist

def pmid2id(thesepmids, others):
    """
    Converts a list of PMIDs to IDs.
    """
    global full_bibdat
    outids = []
    missing_ids = others
    for pmid in thesepmids:
        entry = full_bibdat.get_entry_by_pmid(pmid)
        if entry:
            outids.append(entry['ID'])
        else:
            print(f"PMID {pmid} not found. Searching online...")
            new_entry = findcitation(pmid, 'pmid')
            if new_entry:
                outids.append(new_entry['ID'])
                full_bibdat.add_entry(new_entry)
            else:
                missing_ids.append(pmid)
    return outids, missing_ids

def format_inline(thisid):
    global cited_bibdat
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

def replace_ids_in_text(text):
    """
    Replace occurrences of old IDs with new IDs in the text, using the global regex patterns.
    Returns both the updated text and a report of changes made.
    
    Returns:
        tuple: (updated_text, report) where report is a string describing all changes made
    """
    global full_bibdat
    patterns = [
        markdown_citations,
        latexbrackets
    ]
    combined_pattern = '(' + '|'.join(patterns) + ')'
    report_lines = []    
    def replace_ids_in_match(match):
        citation_block = match.group(0)
        original_block = citation_block  # Keep for debugging
        for old_id, new_id in full_bibdat.get_id_changes().items():
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
            cite(theseids)
            replacedict[b] = mdout(theseids, notfound, outputstyle, flc=flc)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(citedhere['ids'], citedhere['notfound']) # ids added to bib
            pm = remove_duplicates_preserve_order(pm+citedhere['pmids'])
            cite(citedhere['ids'])
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in r:
        theseids, theseothers = parse_wholecitation_block(b, flc)
        cite(theseids)
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
    full_bibdat = BibData()
    cited_bibdat = BibData()

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
    original_fullbib_content = None

    # Read bib files and get ID changes
    if localbibonly:
        local_bibdat, _ = read_bib_file(localbibpath)
        full_bibdat = copy.copy(local_bibdat)
        merge_bibdat_duplicates(full_bibdat)
    else:
        local_bibdat, _ = read_bib_file(localbibpath)
        globalbibfile = os.path.abspath(os.path.expanduser(globalbibfile))
        full_bibdat, original_fullbib_content = read_bib_file(globalbibfile)
        merge_bibdat_duplicates(full_bibdat, local_bibdat)

    text, id_change_report = replace_ids_in_text(text)
    if len(id_change_report)>0:
        print ("\n ID changes in the following citations:")
        print (id_change_report)

    text = readheader(text)[1]
    text = replace_blocks(text, outputstyle, use_whole=wholereference, flc=force_lowercase_citations)
    

    # save local cs.bib file
    # keep any uncited items in localbibdat because the user might want them. But remove duplicates. User can handle this manually. 
    merge_bibdat_duplicates(cited_bibdat, local_bibdat)
    bibdir, bibfilename = os.path.split(localbibpath)
    bibstem = '.'.join(bibfilename.split('.')[:-1])
    localbibpath = os.path.join(bibdir, bibstem + citelabel + "bib")
    print('\nSaving bibliography for this file here:', localbibpath)
    outbib = bibtexparser.dumps(cited_bibdat)
    outbib = make_unicode(outbib)
    with open(localbibpath, "w", encoding="utf-8") as bf:
        bf.write(outbib)

    # Save new global bibliography 
    new_bib_content = serialize_bib_database(full_bibdat)
    if original_fullbib_content is not None:
        if hash_content(original_fullbib_content) != hash_content(new_bib_content):
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


