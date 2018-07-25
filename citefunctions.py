#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8


import platform
import re
import string
import json
import sys
import requests
import difflib
from Bio import Entrez
from Bio import Medline
import json
import xml.etree.ElementTree as ET
import calendar





#-------------
def getconfig():
    with open('config.json') as json_data_file:
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
        print("unrecognised operating system:")
        print(platform.system())
        sys.exit()

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

def get_latex_citation_blocks(inputtext):
    return get_parethesised(inputtext, ['\\cite\{.+?\}'])

def get_md_citation_blocks(inputtext):
    return get_parethesised(inputtext, ['\[.+?\]'])

def get_pmid_citation_blocks(inputtext):
    confirmed_blocks = []
    for b in get_parethesised(inputtext, ['\[.+?\]', '\{.+?\}']):
        if "PMID" in b or "pmid" in b:
            confirmed_blocks.append(b)
    return confirmed_blocks

def split_by_delimiter(thislist, delimiters=[";",",","[","]"]):
    '''
        input is a list of strings
        return flattened list of non-empty, stripped strings, split by delimiters
    '''
    for d in delimiters:
        thislist = [e.split(d) for e in thislist]
        thislist = flatten(thislist)
    thislist = [e.strip() for e in thislist if e!=""]
    return thislist

def get_latex_citations(inputtext):
    pth = get_latex_citation_blocks(inputtext)
    pth = split_by_delimiter(pth)
    pth = [e.replace('cite{','').replace('}','') for e in pth]
    return list(set(pth))
    
def get_md_citations(inputtext):
    pth = get_md_citation_blocks(inputtext)
    pth = split_by_delimiter(pth)
    pth = [x.replace('@','') for x in pth if x.startswith('@')]
    return list(set(pth))

def get_pmid_citations(inputtext):
    pth = get_pmid_citation_blocks(inputtext)
    pth = split_by_delimiter(pth)
    for p in ["PMID: ","PMID :", "pmid:", "pmid :", "pmid: "]:
        pth = [e.replace(p, "PMID:") for e in pth]
    pmidstyle = []
    for x in pth:
        if x.startswith('PMID:'):
            # then assume this is a correctly formatted PMID
            pmidstyle.append(x.replace('PMID:',''))
        else:
            # rescue some integers if they look like PMIDs
            print(x)
            if len(x)>4 and len(x)<10:
                try:
                    int(x)
                except:
                    continue
                pmidstyle.append(x)
    return pmidstyle

def bibadd(thisdb,thisentry):
    if thisentry is not None:
        try:
            thisentry['ENTRYTYPE']
        except:
            thisentry['ENTRYTYPE'] = 'article'
        try:
            thisdb[thisentry]
        except:
            thisdb.entries.append(thisentry)

def convert_separators_to_preferred(thisblock, possible_separators = [',',';',' ','\t'], preferred_separator = ', '):
    for sep in possible_separators:
        thisblock = thisblock.replace(sep, '--|holdingseparator|--')
    blockitems = [x for x in thisblock.split('--|holdingseparator|--') if len(x)>0]
    return preferred_separator.join(blockitems)

def replace_id_with_pmid(thistext, thisid, thispmid):
    ''' 
        return thistext after replacing all occurences of thisid (either a latex- or md-style reference) with thispmid
    '''
    #print 'replacing {} with {}'.format(thisid, thispmid)
    pmidstring = 'PMID:{}'.format(thispmid)
    for block in get_latex_citation_blocks(thistext):
        newblock = block[block.find('{')+1:block.rfind('}')]
        newblock = '['+newblock+']'
        newblock = convert_separators_to_preferred(newblock)
        newblock = newblock.replace(thisid, pmidstring)
        thistext = thistext.replace('\\'+block, newblock)
        thistext = thistext.replace(block, newblock)
    for block in get_md_citation_blocks(thistext):
        newblock = '['+block[1:-1]+']'
        newblock = convert_separators_to_preferred(newblock)
        newblock = newblock.replace('@'+thisid, pmidstring)
        newblock = newblock.replace(thisid, pmidstring)
        thistext = thistext.replace(block, newblock)
    return thistext

def replace_pmid_with_id(thistext, thisid, thispmid, style='tex'):
    pmidstring = 'PMID:{}'.format(thispmid)
    for block in get_pmid_citation_blocks(thistext):
        if block.startswith('{') and block.endswith('}'):
            block = '['+block[1:-1]+']'
        newblock = block[block.find('[')+1:block.rfind(']')]
        if style == 'tex':
            newblock = convert_separators_to_preferred(newblock, preferred_separator=', ')
            newblock = newblock.replace(pmidstring, thisid)
            newblock = '\\cite{' + newblock + '}'
        elif style == 'md':
            if not thisid.startswith('@'):
                thisid = '@'+thisid
            newblock = convert_separators_to_preferred(newblock, preferred_separator='; ')
            newblock = newblock.replace(pmidstring, thisid)
            newblock = '[' + newblock + ']'
        thistext = thistext.replace(block, newblock)
    return thistext

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













