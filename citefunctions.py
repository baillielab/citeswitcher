#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8


import platform
import re
import string
import json
import sys
import os
import subprocess
import requests
import difflib
from Bio import Entrez
from Bio import Medline
import json
import xml.etree.ElementTree as ET
import calendar


'''
def read_blocks():
    should read all blocks, latex, md or pmid, and assign to categories:
    1. formatted latex block with recognised ids
    2. formatted md block with recognised ids
    3. formatted pmid block with recognised ids
    4. block with unrecognised ids - this will be rewritten with [output] and [unmatched]

'''

#-------------
def getconfig():
    cfgfile = os.path.join(os.path.dirname(__file__), 'config.json')
    with open(cfgfile) as json_data_file:
        data = json.load(json_data_file)
    return data
#-------------
config = getconfig()
Entrez.email=config['email']
#-------------
def get_home_dir():
    if platform.system()=="Linux":
        return config['linux_home']
    elif platform.system()=="Darwin":
        return config['mac_home']
    else: 
        print("*** unrecognised operating system (citefunctions.get_home_dir()): {} ***".format(platform.system()))
        sys.exit()

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

def callpandoc(f, out_ext, out_dir='', args="--toc "):
    cmd = 'pandoc {} --filter pandoc-crossref --filter pandoc-citeproc {} -o {}'.format(args, f, os.path.join(out_dir, newext(f, out_ext)))
    print (cmd)
    subprocess.call(cmd, shell=True)

#-------------
def findreplace(inputtext, frdict):
    for f in frd:
        inputtext = inputtext.replace(f, frdict[f])
    return inputtext

#-------------
def readheader(filecontents):
    '''
        Read a valid markdown header
    '''
    t = filecontents.strip()
    t = t.replace('\r','\n')
    lines = t.split('\n')
    lines = [x.strip() for x in lines]
    header = []
    remainder = filecontents
    if lines[0]=='---' and '...' in lines:
        header = lines[1:lines.index('...')]
        remainder = '\n'.join(lines[lines.index('...')+1:])
    return header, remainder

def addheader(filecontents, bibtexfile, cslfilepath='null'):
    '''
        Add components to markdown header. If no header exists, add one.
    '''
    filecontents = filecontents.strip()
    header, remainder = readheader(filecontents)
    headerkeys = [x.split(':')[0] for x in header]
    if 'title' not in headerkeys:
        header.append('title: {}'.format(os.path.split(bibtexfile)[0].replace('.bib','')))
    if 'date' not in headerkeys:
        header.append('date: \\today')
    if 'bibliography' not in headerkeys:
        header.append('bibliography: {}'.format(bibtexfile))
    if 'csl' not in headerkeys and cslfilepath != 'null':
        header.append('csl: {}'.format(cslfilepath))
    return '---\n{}\n...\n{}'.format('\n'.join(header), remainder)

#-------------
def make_dictionaries(bib_db):
    pmids = {}
    ids = {}
    for entry in bib_db:
        # biblatex id dict
        try:
            ids[entry['ID']]
            print("duplicate ID in biblatex database:{}".format(entry["ID"]))
            if 'pmid' in entry:
                # replace this entry with a new one that has a PMID
                ids[entry['ID']]=entry
        except:
            ids[entry['ID']]=entry
        # pmid id dict
        try:
            entry["pmid"]
        except:
            continue
        try:
            pmids[entry["pmid"]]
            print("duplicate PMID in biblatex database:{}".format(entry["pmid"]))
        except:
            pass
        pmids[entry["pmid"]] = entry
    return pmids, ids

#-------------

def flatten(thislist):
    listcount = [1 for x in thislist if type(x)==type([])]
    if sum(listcount)==0:
        return thislist
    else:
        return [item for sublist in thislist for item in sublist] #flatten list

def make_unicode(input):
    if type(input) != str:
        input =  input.decode('utf-8')
        return input
    else:
        return input

def get_parethesised(thistext,parentheses=['\[.+?\]', '\{.+?\}']):
    if type(parentheses)!=type([]):
        print ("Error, parentheses must be in a list")
    outlist = []
    for p in parentheses:
        try:
            outlist += re.findall(p,thistext)
        except:
            continue
    return outlist

def remove_parentheses(thistext):
    thistext = thistext.strip()
    for opener in ('(','[','{','\\cite{'):
        if thistext.startswith(opener):
            thistext = thistext[len(opener):]
    for closer in (')','}',']'):
        if thistext.endswith(closer):
            thistext = thistext[:-1]
    return thistext

def get_latex_citation_blocks(inputtext):
    return get_parethesised(inputtext, ['\\\cite\{.+?\}'])

def get_md_citation_blocks(inputtext):
    confirmed_blocks = []
    for b in get_parethesised(inputtext, ['\[.+?\]']):
        if "@" in b:
            confirmed = True
            for crossreflabel in ['@fig:','@sec:','@tbl:','@eq:']:
                if remove_parentheses(b).startswith(crossreflabel):
                    confirmed = False
            if confirmed:
                confirmed_blocks.append(b)
    return confirmed_blocks

def get_pmid_citation_blocks(inputtext):
    confirmed_blocks = []
    for theseparetheses in ['\[.+?\]', '\{.+?\}', '\(.+?\)']:
        for b in get_parethesised(inputtext, [theseparetheses]):
            if "PMID" in b or "pmid" in b:
                print (b)
                confirmed_blocks.append(b)
                # remove this block so that nested blocks
                # are only counted once
                inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks

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

def get_all_id_citations(inputtext, style='latex'):
    pth = get_latex_citation_blocks(inputtext)
    pth += get_md_citation_blocks(inputtext)
    citationlist = []
    for block in pth:
        citationlist += parse_id_block(block)[0]
    return list(set(citationlist))

def get_all_pmid_citations(inputtext):
    pth = get_pmid_citation_blocks(inputtext)
    citationlist = []
    for block in pth:
        citationlist += parse_pmid_block(block)[0]
    return list(set(citationlist))

def parse_id_block(thisblock):
    thisblock = remove_parentheses(thisblock)
    pth = split_by_delimiter(thisblock)
    return list(set(pth)), []

def parse_pmid_block(thisblock):
    thisblock = remove_parentheses(thisblock)
    # better to use regex for this
    for p in ["PMID ", "PMID: ", "pmid:", "pmid :", "pmid: "]:
        thisblock = thisblock.replace(p, "PMID:")
    thisblock = split_by_delimiter(thisblock)
    pmidstyle = []
    notpmid = []
    for x in thisblock:
        if x.startswith('PMID'):
            # then assume this is a correctly formatted PMID
            x = x.replace('PMID:','')
            x = x.replace('PMID','')
            x = x.strip()
            if len(x)>0:
                pmidstyle.append(x)
        else:
            print ("failed to get pmid:{}".format(x))
            notpmid.append(x)
    return list(set(pmidstyle)), list(set(notpmid))

#------------

def id2pmid(theseids, id_db):
    pmidlist = []
    notpmidlist = []
    for thisid in theseids:
        if thisid.startswith("PMID:"):
            pmidlist.append(thisid)
        try:
            pmidlist.append(ids[thisid]['PMID'])
        except:
            notpmidlist.append(thisid)
    return pmidlist, notpmidlist

def pmid2id(thesepmids, others, pmid_db, id_db):
    outids = []
    missing_ids = others
    for pmid in thesepmids:
        try:
            outids.append(pmid_db[pmid]['ID'])
        except:
            try: 
                id_db[pmid] # if no error, this is an id, not a pmid
                outids.append(pmid)
            except:
                missing_ids.append(pmid)
    return outids, missing_ids

def pmidout(pmidlist, notpmidlist):
    blockstring = ''
    if len(pmidlist) > 0:
        blockstring = '[' + ', '.join(["PMID:{}".format(x) for x in pmidlist]) + ']'
    if len(notpmidlist) > 0:
        blockstring += '[*' + ', '.join(notpmidlist) + ']'
    return blockstring

def texout(theseids, thesemissing=[]):
    blockstring = ''
    if len(theseids) > 0:
        blockstring += "\\cite\{{}\}".format(', '.join(theseids))
    if len(thesemissing) > 0:
        blockstring += "[**{}]".format(', '.join(thesemissing))
    return blockstring

def mdout(theseids, thesemissing=[]):
    blockstring = ''
    if len(theseids) > 0:
        blockstring += "[{}]".format(', '.join(["@{}".format(x) for x in theseids]))
    if len(thesemissing) > 0:
        blockstring += "[***{}]".format(', '.join(thesemissing))
    return blockstring

def replace_blocks(thistext, pmid_db, id_db, outputstyle="md"):
    # pmid first as they are the most likely to have errors
    p = get_pmid_citation_blocks(thistext)
    l = [x for x in get_latex_citation_blocks(thistext) if x not in p]
    m = [x for x in get_md_citation_blocks(thistext) if x not in p and x not in l]
    replacedict = {}
    for b in m:
        theseids = parse_id_block(b)[0]
        if outputstyle=='tex':
            replacedict[b] = texout(theseids)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(theseids, id_db)
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in l:
        theseids = parse_id_block(b)[0]
        if outputstyle=='md':
            replacedict[b] = mdout(theseids)
        elif outputstyle=='pmid':
            pm, notpm = id2pmid(theseids, id_db)
            replacedict[b] = pmidout(pm, notpm)
        else:
            continue
    for b in p:
        thesepmids, theseothers = parse_pmid_block(b)
        theseids, notfound = pmid2id(thesepmids, theseothers, pmid_db, id_db)
        if outputstyle == 'md':
            replacedict[b] = mdout(theseids, notfound)
        elif outputstyle=='tex':
            replacedict[b] = texout(theseids, notfound)
        else:
            continue
    for b in replacedict:
        print ("{:>70} ==> {}".format(b, replacedict[b]))
        thistext = thistext.replace(b, replacedict[b])
    return thistext

#-----------------

def bibadd(thisdb, thisentry):
    if thisentry is not None:
        try:
            thisentry['ENTRYTYPE']
        except:
            thisentry['ENTRYTYPE'] = 'article'
        try:
            thisdb[thisentry]
        except:
            thisdb.entries.append(thisentry)

def find_similar_keys(this_string, thisdict):
    '''
        search the keys of thisdict for a string similar to this_string
    '''
    topscore=-1
    bestmatch='none'
    for i, key in enumerate(thisdict.keys()):
        sim = difflib.SequenceMatcher(None, this_string, key).ratio()
        if sim>topscore:
            bestmatch = key
            topscore = sim
    return bestmatch


#------ PUBMED FUCTIONS -------

def id_translator(thisid):
    '''
    search pubmed API; return dict of {'pmid':PMID, 'pmcid':PMCID, 'doi':DOI, 'error':'errormessage'} where available
    '''
    # PMCID must start with PMC
    searchstring = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool={}&email={}&ids={}&format=json".format(config['toolname'], Entrez.email, thisid)
    r = requests.get(searchstring)
    r = json.loads(r.text)
    dopubmed = False
    if 'records' in r:
        if 'status' not in r['records'][0]: # then this is not an error
            return r['records'][0]
        elif r['records'][0]['status']=='error':
            dopubmed = True
    elif 'status' in r:
        if r['status']=='error':
            dopubmed = True
    else:
        return {'error':'fail: parsing API output: {}'.format(r)}
    # if the API{'error': fail:ed to find this one, search pubmed for it instead}
    if dopubmed:
        try:
            r2 = search_pubmed(thisid)
        except:
            return {'error':'fail: pubmed search'}
        if len(r2)==1: # only one record found
            pmid = r2[0]
            try:
                doi = get_doi_from_pmid(pmid)
            except:
                return {'error':"fail: doi"}
            return {'pmid':pmid,'doi':doi}
        else:
            return {'error':"fail: not unique pmid {}".format(thisid)}
    else:
        return {'error':"fail: not found, not trying pubmed"}
    
def search_pubmed(search_string):
    '''return a list of pmids for articles found by a search string'''
    max_returns = 500
    days = "all"
    handle = Entrez.esearch(db="pubmed", term=search_string, retmax=max_returns, rettype="medline")
    r = Entrez.read(handle)
    return r['IdList']

def get_links_from_pmid(pmid):
    return Entrez.elink(dbfrom="pubmed",id=pmid,cmd="prlinks")

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

def p2b(thispmid):
    ''' by Nick Loman '''

    ## Fetch XML data from Entrez.
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    pmids = [str(thispmid)]
    url = '{}?db=pubmed&id={}&rettype=abstract'.format(efetch, ','.join(pmids))
    r = requests.get(url)
    ##print(r.text) # to examine the returned xml

    ## Loop over the PubMed IDs and parse the XML using https://docs.python.org/2/library/xml.etree.elementtree.html
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
        ids = PubmedArticle.findall('./PubmedData/ArticleIdList/ArticleId')
        for thisid in ids:
            if thisid.attrib['IdType'] == 'pmc':
                PMCID = thisid
            elif thisid.attrib['IdType'] == 'doi':
                DOI = thisid
        # format author list
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

        # Print the bibtex formatted output.
        bib = {}
        try:
            bib["ID"] = '{}{}{}'.format(authors[0].split(',')[0], ''.join([x for x in ArticleTitle.text.split(' ') if len(x)>3][:3]), Year)
        except IndexError:
            print ('IndexError', pmids, file=sys.stderr, flush=True)
        except AttributeError:
            print ('AttributeError', pmids, file=sys.stderr, flush=True)
        bib["Author"] = ' AND '.join(authors)
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
        bib["ISSN"] = ISSN.text
        bib["pmid"] = PMID.text
        return bib













