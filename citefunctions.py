#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

import os
import re
import sys
import json
import difflib
import calendar
import requests
import platform
import subprocess
from Bio import Entrez
from Bio import Medline
import xml.etree.ElementTree as ET
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
from bibtexparser.latexenc import latex_to_unicode
import latexchars
#-------------
import bib
#-------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#-------------
def getconfig(cfgfile):
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
    return data
#-------------
markdown_labels_to_ignore = [
    '@fig:','@sec:','@tbl:','@eq:',
    '#fig:','#sec:','#tbl:','#eq:'
    ]
#-------------
#regexes
squarebrackets = '\[[\s\S]+?\]'
curlybrackets = '\{[\s\S]+?\}'
roundbrackets = '\([\s\S]+?\)'
latexbrackets = '\\\cite\{[\s\S]+?\}'
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

def callpandoc(f, out_ext, out_dir='', args="", yaml=""):
    cmd = 'pandoc --pdf-engine=xelatex {} --filter pandoc-crossref --filter pandoc-citeproc {} {} -o {}'.format(args, yaml, f, os.path.join(out_dir, newext(f, out_ext)))
    print (cmd)
    subprocess.call(cmd, shell=True)

def read_bib_files(bibfiles):
    bfs = ""
    bibdb = BibDatabase()
    for bibfile in bibfiles:
        if os.path.exists(bibfile):
            # read bibtex file
            try:
                size = os.path.getsize(bibfile)
            except:
                print ("bib file not found at {}".format(bibfile))
                sys.exit()
            if size > 0:
                with open(bibfile) as bf:
                    bfs += bf.read()
            else:
                print ("bib file empty: {}".format(bibfile))
    bibdb = bibtexparser.bparser.BibTexParser(common_strings=True, homogenize_fields=True).parse(bfs)
    return bibdb

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
        inputtext = inputtext.replace(f, frdict[f])
    return inputtext

#-------------
def getyaml(filepath):
    with open(filepath) as f:
        y = readheader(f.read())
    yml = {}
    for line in y[0]:
        line = line.split(": ")
        if len(line)>1:
            yml[line[0]]=line[1]
    return yml


def readheader(filecontents):
    '''
        Read a valid markdown header
        Input is full text of file
        Returns list of header items + full text of remainder
    '''
    t = filecontents.strip()
    t = t.replace('\r','\n')
    lines = [x for x in t.split('\n')]
    header = []
    remainder = filecontents
    if lines[0]=='---':
        if '...' in lines:
            header = lines[1:lines.index('...')]
            remainder = '\n'.join(lines[lines.index('...')+1:])
        elif '---' in lines[1:]:
            header = lines[1:lines.index('---')]
            remainder = '\n'.join(lines[lines[1:].index('---')+1:])
    return header, remainder

def addheader(filecontents, bibtexfile, cslfilepath='null'):
    '''
        Add components to markdown header. If no header exists, add one.
    '''
    filecontents = filecontents.strip()
    header, remainder = readheader(filecontents)
    headerkeys = [x.split(':')[0] for x in header]
    '''
    if 'title' not in headerkeys:
        header.append('title: {}'.format(os.path.split(bibtexfile)[0].replace('.bib','')))
    if 'date' not in headerkeys:
        header.append('date: \\today')
    '''
    if 'csl' not in headerkeys and cslfilepath != 'null':
        header.append('csl: {}'.format(cslfilepath))
    if 'bibliography' not in headerkeys:
        header.append('bibliography: {}'.format(bibtexfile))
    else:
        header[header.index([y for y in header if y.startswith("bibliography")][0])] = 'bibliography: {}'.format(bibtexfile)
    return '---\n{}\n...\n{}'.format('\n'.join(header), remainder)

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
    inputstring = inputstring.replace(u"\ufeff", "") # sometimes appears at start of file
    return inputstring

def get_parethesised(thistext,parentheses=[squarebrackets, curlybrackets]):
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
        nested_out = [x for x in nested_out if not x[1:].startswith(thislabel)]
    return nested_out

def remove_parentheses(thistext):
    '''
        remove matched parentheses at both ends
        and remove newlines
    '''
    thistext = thistext.strip()
    for opener in ('(','[','{','\\cite{'):
        if thistext.startswith(opener):
            thistext = thistext[len(opener):]
    for closer in (')','}',']'):
        if thistext.endswith(closer):
            thistext = thistext[:-1]
    thistext = thistext.replace('\n',' ')
    thistext = thistext.replace('\r',' ')
    return thistext

def split_by_delimiter(this_string, delimiters=[";",",","[","]"," ","\n","\t","\r"]):
    '''
        return flattened list of non-empty, stripped strings, split by delimiters
    '''
    thislist = this_string.split()
    for d in delimiters:
        thislist = [e.split(d) for e in thislist]
        thislist = flatten(thislist)
    thislist = [e.strip() for e in thislist if e!=""]
    return thislist

# get citation blocks functions

def get_pmid_citation_blocks(inputtext):
    confirmed_blocks = []
    for theseparetheses in [squarebrackets, curlybrackets, roundbrackets]:
        for b in get_parethesised(inputtext, [theseparetheses]):
            if "PMID" in b or "pmid" in b:
                confirmed_blocks.append(b)
                # remove this block so that nested blocks
                # are only counted once
                inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks, inputtext

def get_latex_citation_blocks(inputtext):
    confirmed_blocks = []
    for b in get_parethesised(inputtext, [latexbrackets]):
        confirmed_blocks.append(b)
        inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks, inputtext

def get_md_citation_blocks(inputtext):
    confirmed_blocks = []
    for b in get_parethesised(inputtext, [squarebrackets]):
        if "@" in b:
            confirmed = True
            for crossreflabel in markdown_labels_to_ignore:
                if remove_parentheses(b).startswith(crossreflabel):
                    confirmed = False
            if confirmed:
                confirmed_blocks.append(b)
                inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks, inputtext

def get_doi_citation_blocks(inputtext):
    confirmed_blocks = []
    for theseparetheses in [squarebrackets, curlybrackets, roundbrackets]:
        for b in get_parethesised(inputtext, [theseparetheses]):
            if "doi:" in b or "DOI:" in b or "https://doi.org/" in b:
                confirmed_blocks.append(b)
                inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks, inputtext

def get_searchstring_from_wholecitation(wc):
    '''
        return the longest sentence that doesn't look like a list of authors or a citation
    '''
    authorsep = '[A-Z], [A-Z]' # regex for [capital, comma, space, capital]
    lendict = {x:len(x) for x in remove_parentheses(wc).split('.')}
    for x in sorted(iter(lendict.items()), key=lambda k_v: (k_v[1],k_v[0]), reverse=True):
        if len(x[0]) > 0:
            authoriness = float(len(re.findall(authorsep, x[0]))/len(x[0]))
            if authoriness < 0.02:
                if ":" in x[0] and ";" in x[0]: # probably not a title
                    continue
                else:
                    return x[0]

def get_wholereference_citation_blocks(inputtext):
    '''
        finds wholecitations AND doi
        return a list of citation blocks
    '''
    confirmed_blocks = []
    askedalready = {}
    for theseparetheses in [squarebrackets, curlybrackets]:
        for b in get_parethesised(inputtext, [theseparetheses]):
            if "." in b or ":" in b:
                try:
                    new_entry = askedalready[b.replace('\n',' ')]
                except:
                    new_entry = parse_wholecitation_block(b)
                if new_entry is not None:
                    confirmed_blocks.append(str(b))
                    inputtext = inputtext.replace(b, '---citationblockremoved---')
                    askedalready[b.replace('\n',' ')] = new_entry
    return confirmed_blocks, inputtext

# parse citation blocks functions

def parse_id_block(thisblock):
    ''' take a block of text, and return two lists of ids '''
    thisblock = remove_parentheses(thisblock)
    pth = split_by_delimiter(thisblock)
    pth = [x.replace("@","") for x in pth if x.startswith('@')]
    theseids = list(set(pth))
    bib.cite(theseids)
    return theseids, []

def parse_pmid_block(thisblock):
    ''' take a block of text, and return two lists of ids '''
    thisblock = remove_parentheses(thisblock)
    # better to use regex for this
    for p in ["PMID ", "PMID: ", "pmid:", "pmid :", "pmid: "]:
        thisblock = thisblock.replace(p, "PMID:")
    thisblock = split_by_delimiter(thisblock)
    pmidstyle = []
    notpmid = []
    for x in thisblock:
        if x.startswith('PMID'):
            # then assume this is a correctly formatted PMID
            x = x.replace('PMID:','')
            x = x.replace('PMID','')
            x = x.strip()
            if len(x)>0:
                pmidstyle.append(x)
        else:
            try:
                found = get_article_from_pmid(x)
                if found['PMID'] == x:
                    pmidstyle.append(x)
                else:
                    raise ValueError('found more than one')
            except:
                print ("failed to get pmid:{}".format(x))
                notpmid.append(x)
    return list(set(pmidstyle)), list(set(notpmid))

def parse_doi_block(thisblock):
    ''' take a block of text, and return two lists of ids '''
    doistyle = []
    notdoi = []
    d = remove_parentheses(thisblock).replace("DOI:","doi:")
    d = d.replace("https://doi.org/","doi:")
    d = d.replace("doi: ", "doi:")
    d = d.replace("https://doi.org/", "")
    d = d.replace("http://doi.org/", "")
    d = d.replace(",", " ")
    dlist = [x.replace("doi:","") for x in d.split(' ') if x.startswith('doi:')]
    for x in dlist:
        found = findcitation(x, 'doi')
        if found is not None:
            doistyle.append(found['ID'])
        else:
            notdoi.append(x)
    return doistyle, notdoi

def parse_wholecitation_block(thisblock):
    ''' take a block of text, and return two lists of ids '''
    outids = []
    # try a doi search first in case there's a doi in there
    d = parse_doi_block(thisblock)[0]
    if len(d)>0:
        outids = d
    else:
        # if that doesn't work, try a title search
        title = get_searchstring_from_wholecitation(thisblock)
        print ("searching:", title)
        found = findcitation(title, 'title', additionalinfo=thisblock)
        print ("title:", title, found)
        if found is not None:
            outids.append(found['ID'])
    return outids, []  # there is only ever one

#------------

def findcitation(info, infotype='pmid', additionalinfo='', ):
    '''
        search the bibdat first
        search online in a variety of ways
        if found, use bib.new to add it to global db
        return a single bibtex entry
        or None if not found or in doubt
    '''
    if infotype == 'pmid':
        try:
            # print ("searching for {} in bibdat".format(info))
            return bib.pmids[info]
        except:
            # print ("{} not found in bib.pmids: ".format(info))
            # print (bib.pmids.keys())
            pass
        pub = p2b([str(info)])
        if len(pub) > 0:
            if pub[0] != 'null' and pub[0] != None:
                print ("PMID:{} found online".format(info))
                bib.new(pub[0])
                return pub[0]
        print ("PMID:{} NOT FOUND ON PUBMED".format(info))
        return None
    elif infotype == 'doi':
        try:
            return bib.dois[info]
        except:
            pass
        print ("this is a doi: {}".format(info))
        pub = search_pubmed(info, "doi")
        if len(pub) == 1:
            pubent = p2b(pub[0])
            if len(pubent) > 0:
                if pubent[0] != 'null' and pubent[0] != None:
                    bib.new(pubent[0])
                    return pubent[0]
        print ("doi search failure: {}\nPubmed return: {}".format(info, pub))
        return None
    elif infotype == 'title':
        if additionalinfo == '':
            additionalinfo = title
        pub = search_pubmed(info, "title")
        if len(pub) == 1:
            pmid = pub[0]
            pubent = p2b(pmid)
            if len(pubent) > 0:
                if pubent[0] != 'null' and pubent[0] != None:
                    q = input ("--------------\n\
Reference block (PMID:{}) found. \
Please check that input:\n\n{}\n\n\
Is the same as the found citation:\n{}\n\n\
Enter y/n\
                        ".format(pmid, additionalinfo, '\n'.join( ["{:>12}:    {}".format(x,pubent[0][x]) for x in pubent[0]]) ))
                    q = q.strip().upper()
                    if q == "Y":
                        print ('--confirmed--')
                        bib.new(pubent[0])
                        return pubent[0]
        return None


def id2pmid(theseids):
    '''
        input is a list of ids
    '''
    pmidlist = []
    notpmidlist = []
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
    # add to the outputdatabase
    bib.cite([bib.pmids[x]['ID'] for x in pmidlist])
    # make a blockstring
    blockstring = ''
    if len(pmidlist) > 0:
        blockstring = '[' + ', '.join(["PMID:{}".format(x) for x in pmidlist]) + ']'
        if len(notpmidlist) > 0:
            blockstring += '[*' + ', '.join(notpmidlist) + ']'
    else:
        blockstring = 'null'
    return blockstring

def mdout(theseids, thesemissing=[], outputstyle="md"):
    # add to the outputdatabase
    bib.cite(theseids)
    # make a blockstring
    blockstring = ''
    if len(theseids) == 0:
        blockstring = 'null'
    elif outputstyle == "md":
        blockstring += "[{}]".format('; '.join(["@{}".format(x) for x in theseids]))
        if len(thesemissing) > 0:
            blockstring += "[***{}]".format(', '.join(thesemissing))
    elif outputstyle == "tex":
        bib.cite(theseids)
        blockstring += "\\cite\{{}\}".format(', '.join(theseids))
        if len(thesemissing) > 0:
            blockstring += "[**{}]".format(', '.join(thesemissing))
    elif outputstyle == "inline":
        blockstring += "({})".format(', '.join([format_inline(x) for x in theseids]))
        if len(thesemissing) > 0:
            blockstring += "[***{}]".format(', '.join(thesemissing))
    return blockstring

def replace_blocks(thistext, outputstyle="md"):
    # pmid first as they are the most likely to have errors
    workingtext = thistext
    p, workingtext = get_pmid_citation_blocks(workingtext)
    l, workingtext = [x for x in get_latex_citation_blocks(workingtext) if x not in p]
    m, workingtext = [x for x in get_md_citation_blocks(workingtext) if x not in p+l]
    d, workingtext = [x for x in get_doi_citation_blocks(workingtext) if x not in p+l+m]
    r, workingtext = [x for x in get_wholereference_citation_blocks(workingtext) if x not in p+l+m+d]
    print ("Number blocks using pmid:{}".format(len(p)))
    print ("Number blocks using latex:{}".format(len(l)))
    print ("Number blocks using markdown:{}".format(len(m)))
    print ("Number blocks using doi:{}".format(len(d)))
    print ("Number blocks using wholeref:{}".format(len(r)))
    replacedict = {}
    for b in p:
        thesepmids, theseothers = parse_pmid_block(b)
        theseids, notfound = pmid2id(thesepmids, theseothers) # ids added to bib
        if outputstyle == 'md' or outputstyle=='tex' or outputstyle=='inline':
            replacedict[b] = mdout(theseids, notfound, outputstyle)
        else:
            continue
    for b in l:
        theseids, notfound = parse_id_block(b)  # ids added to bib
        if outputstyle=='md' or outputstyle=='inline':
            replacedict[b] = mdout(theseids, notfound, outputstyle)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(theseids) # ids added to bib
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in m:
        theseids, notfound = parse_id_block(b)  # ids added to bib
        if outputstyle=='tex' or outputstyle=='inline':
            replacedict[b] = mdout(theseids, notfound, outputstyle)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(theseids) # ids added to bib
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in d:
        theseids, notfound = parse_doi_block(b)  # ids added to bib
        if outputstyle == 'md' or outputstyle=='tex' or outputstyle=='inline':
            replacedict[b] = mdout(theseids, notfound, outputstyle)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(theseids) # ids added to bib
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in r: 
        theseids, theseothers = parse_wholecitation_block(b)
        if outputstyle == 'md' or outputstyle=='tex' or outputstyle=='inline':
            replacedict[b] = mdout(theseids, theseothers, outputstyle)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(theseids) # ids added to bib
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in replacedict:
        if replacedict[b] == 'null':
            print ("{:>70} left alone".format(b))
            continue
        print ("{:>70} ==> {}".format(b, replacedict[b]))
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
    handle = Entrez.esearch(db="pubmed", term=search_string, retmax=max_returns, rettype="medline", field=restrictfields)
    r = Entrez.read(handle)
    return r['IdList']

def get_links_from_pmid(pmid):
    return Medline.parse(Entrez.elink(dbfrom="pubmed",id=pmid,cmd="prlinks"))

def get_article_from_pmid(pmid):
    output = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
    details = Medline.parse(output)
    thisarticle = next(details)
    return thisarticle

def get_doi_from_pmid(pmid):
    articledetails = get_article_from_pmid(pmid)
    doi_list = [string.strip(x.replace('[doi]','')) for x in articledetails['AID'] if x.endswith('[doi]')]
    if len(doi_list)>0:
        return doi_list[0]
    else:
        raise ValueError('This pmid ({}) was not found.'.format(pmid))


# --------------------
# by Nick Loman

def p2b(pmidlist):
    ''' by Nick Loman '''

    if type(pmidlist) != list:
        pmidlist = [pmidlist]

    ## Fetch XML data from Entrez.
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    url = '{}?db=pubmed&id={}&rettype=abstract'.format(efetch, ','.join(pmidlist))
    r = requests.get(url)
    ##print(r.text) # to examine the returned xml
    ## Loop over the PubMed IDs and parse the XML using https://docs.python.org/2/library/xml.etree.elementtree.html
    bibout = []
    root = ET.fromstring(r.text)
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
        try:
            for _ in (PMID.text, Volume.text, Title.text, ArticleTitle.text, MedlinePgn.text, Abstract.text, ''.join(authors)):
                if _ is None:
                    continue
                assert '{' not in _, _
                assert '}' not in _, _
        except AttributeError:
            pass

        # make the bibtex formatted output.
        bib = {}
        if len(authors)>0:
            authorname = authors[0].split(',')[0]
        else:
            authorname = ''
        titlewords = [x for x in ArticleTitle.text.split(' ') if len(x)>3]
        if len(titlewords)>2:
            titlestring = ''.join(titlewords[:3])
        elif len(titlewords)>0:
            titlestring = ''.join(titlewords)
        else:
            titlestring = ''
        if len(authorname+titlestring)==0:
            titlestring = "PMID{}_".format(PMID.text)
        new_id = '{}{}{}'.format(authorname, titlestring, Year)
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





