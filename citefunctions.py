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
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
import latexchars
#-------------
def getconfig(cfgfile):
    with open(cfgfile) as json_data_file:
        data = json.load(json_data_file)
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

def callpandoc(f, out_ext, out_dir='', args=""):
    cmd = 'pandoc {} --number-sections --filter pandoc-crossref --filter pandoc-citeproc {} -o {}'.format(args, f, os.path.join(out_dir, newext(f, out_ext)))
    print (cmd)
    subprocess.call(cmd, shell=True)

def read_bib_file(bibfile):
    # read bibtex file
    try:
        size = os.path.getsize(bibfile)
    except:
        print ("bib file not found at {}".format(bibfile))
        sys.exit()
    if size > 0:
        with open(bibfile) as bf:
            bibdb = bibtexparser.bparser.BibTexParser(common_strings=True, homogenize_fields=True).parse_file(bf)
    else:
        bibdb = BibDatabase()
        print ("bib file empty: {}".format(bibfile))
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
    ids = [key for key, value in sorted(iter(sorter.items()), key=lambda k_v: (k_v[1],k_v[0]), reverse=True)]
    thisdb.entries = [thisdb.entries_dict[thisid] for thisid in ids]

#-------------
def findreplace(inputtext, frdict):
    for f in frdict:
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

def make_unicode(inputstring):
    if type(inputstring) != str:
        inputstring =  inputstring.decode('utf-8')
        return inputstring
    else:
        return inputstring

def get_parethesised(thistext,parentheses=['\[.+?\]', '\{.+?\}']):
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
    return nested_out

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
            for crossreflabel in markdown_labels_to_ignore:
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
                confirmed_blocks.append(b)
                # remove this block so that nested blocks
                # are only counted once
                inputtext = inputtext.replace(b, '---citationblockremoved---')
    return confirmed_blocks

def get_wholereference_citation_blocks(inputtext):
    confirmed_blocks = {}
    for theseparetheses in ['\[.+?\]', '\{.+?\}']:
        for b in get_parethesised(inputtext, [theseparetheses]):
            if "." in b or ":" in b:
                if "@" not in b and "PMID" not in b and "pmid" not in b and "#" not in b:
                    # try searching for the whole ref
                    pub = search_pubmed(b)
                    if len(pub) != 1:
                        # try searching just for the title (the longest sentence in the ref)
                        lendict = {x:len(x) for x in remove_parentheses(b).split('.')}
                        title = sorted(iter(lendict.items()), key=lambda k_v: (k_v[1],k_v[0]), reverse=True)[0][0]
                        pub = search_pubmed(title, "title")
                    if len(pub) == 1:
                        pmid = pub[0]
                        p = p2b(pmid)
                        if len(p) > 0:
                            p = p[0]
                            q = input ("--------------\n\
Reference block (PMID:{}) found. \
Please check that input:\n\n{}\n\n\
Is the same as the found citation:\n{}\n\n\
Enter y/n\
                                ".format(pmid, b, '\n'.join( ["{:>12}:    {}".format(x,p[x]) for x in p]) ))
                            q = q.strip().upper()
                            if q == "Y":
                                print ('--confirmed--')
                                confirmed_blocks[pmid] = b
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
    pth = [x.replace("@","") for x in pth if x.startswith('@')]
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

#------------

def id2pmid(theseids, id_db):
    pmidlist = []
    notpmidlist = []
    for thisid in theseids:
        if thisid.startswith("PMID:"):
            pmidlist.append(thisid)
        if 'PMID' in id_db[thisid]:
            pmidlist.append(id_db[thisid]['PMID'])
        elif 'pmid' in id_db[thisid]:
            pmidlist.append(id_db[thisid]['pmid'])
        else:
            print ("pmid not found in bib file: {}. Searching online...".format(thisid))
            pub = search_pubmed(id_db[thisid]['doi'], "doi")
            if len(pub) == 1:
                pmid = pub[0]
                pmidlist.append(pmid)
                id_db[thisid]['PMID'] = pmid
            else:
                pub = search_pubmed(id_db[thisid]['title'], "title")
                if len(pub) == 1:
                    pmid = pub[0]
                    p = p2b(pmid)
                    print ("maybe rescued? Not included as this code not written yet...")
                    print (p)
                notpmidlist.append(thisid)
                continue
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
    else:
        blockstring = 'null'
    return blockstring

def texout(theseids, thesemissing=[]):
    blockstring = ''
    if len(theseids) > 0:
        blockstring += "\\cite\{{}\}".format(', '.join(theseids))
        if len(thesemissing) > 0:
            blockstring += "[**{}]".format(', '.join(thesemissing))
    else:
        blockstring = 'null'
    return blockstring

def mdout(theseids, thesemissing=[]):
    blockstring = ''
    if len(theseids) > 0:
        blockstring += "[{}]".format('; '.join(["@{}".format(x) for x in theseids]))
        if len(thesemissing) > 0:
            blockstring += "[***{}]".format(', '.join(thesemissing))
    else:
        blockstring = 'null'
    return blockstring

def replace_blocks(thistext, pmid_db, id_db, outputstyle="md"):
    # pmid first as they are the most likely to have errors
    p = get_pmid_citation_blocks(thistext)
    l = [x for x in get_latex_citation_blocks(thistext) if x not in p]
    m = [x for x in get_md_citation_blocks(thistext) if x not in p and x not in l]
    wr = get_wholereference_citation_blocks(thistext) # dict because only found articles included
    r = {x:wr[x] for x in wr if wr[x] not in p and wr[x] not in l and wr[x] not in m}
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
    for pmid in r:
        b = r[pmid]
        thesepmids = [pmid]
        theseothers = []
        if outputstyle=='pmid':
            replacedict[b] = pmidout([pmid],[])
        else:
            theseids, notfound = pmid2id(thesepmids, theseothers, pmid_db, id_db)
            if outputstyle == 'md':
                replacedict[b] = mdout(theseids, notfound)
            elif outputstyle=='tex':
                replacedict[b] = texout(theseids, notfound)
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

def bibadd(thisdb, thisentry):
    '''
        add new reference at START of entries
    '''
    if thisentry is not None:
        try:
            thisentry['ENTRYTYPE']
        except:
            thisentry['ENTRYTYPE'] = 'article'
        try:
            thisdb[thisentry]
        except:
            thisentry = latexchars.cleanbib(thisentry)
            thisdb.entries = [thisentry] + thisdb.entries

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

def search_pubmed(search_string, restrictfields=""):
    '''return a list of pmids for articles found by a search string'''
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

def p2b(pmids):
    ''' by Nick Loman '''

    if type(pmids) != list:
        pmids = [pmids]

    ## Fetch XML data from Entrez.
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    url = '{}?db=pubmed&id={}&rettype=abstract'.format(efetch, ','.join(pmids))
    r = requests.get(url)
    ##print(r.text) # to examine the returned xml

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
        new_id = latexchars.replace_accents(new_id, mode="biblatex")
        bib["ID"] = new_id
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
        bib["ISSN"] = ISSN.text
        bib["pmid"] = PMID.text

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





