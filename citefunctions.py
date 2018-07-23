#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8


import platform, re, string, json, sys
import requests
import difflib
from Bio import Entrez
from Bio import Medline
import json

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
        print "unrecognised operating system:"
        print platform.system()
        sys.exit()

def flatten(thislist):
    listcount = [1 for x in thislist if type(x)==type([])]
    if sum(listcount)==0:
        return thislist
    else:
        return [item for sublist in thislist for item in sublist] #flatten list

def make_unicode(input):
    if type(input) != unicode:
        input =  input.decode('utf-8')
        return input
    else:
        return input

def get_parethesised(thistext,parentheses=['\[.+?\]', '\{.+?\}']):
    if type(parentheses)!=type([]):
        print ("Error, parentheses must be in a list")
    outlist = []
    for p in parentheses:
        outlist += re.findall(p,thistext)
    return outlist

def get_latex_citation_blocks(inputtext):
    return get_parethesised(inputtext, ['\\cite\{.+?\}'])

def get_md_citation_blocks(inputtext):
    return get_parethesised(inputtext, ['\[.+?\]', '\{.+?\}'])

def get_pmid_citation_blocks(inputtext):
    return get_parethesised(inputtext, ['\[.+?\]', '\{.+?\}'])

def split_by_delimiter(thislist, delimiters=[";",",","[","]"]):
    '''
        input is a list of strings
        return flattened list of non-empty, stripped strings, split by delimiters
    '''
    for d in delimiters:
        thislist = [string.split(e,d) for e in thislist]
        thislist = flatten(thislist)
    thislist = [string.strip(e) for e in thislist if e!=""]
    return thislist

def get_latex_citations(inputtext):
    pth = get_latex_citation_blocks(inputtext)
    pth = split_by_delimiter(pth)
    pth = [e.replace('cite{','').replace('}','') for e in pth]
    return list(set(pth))
    
def get_md_citations(inputtext):
    pth = get_md_citation_blocks(inputtext)
    pth = split_by_delimiter(pth)
    pth = [x.replace('@','') for x in pth if x.startswith(u'@')]
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
            print x
            if len(x)>4 and len(x)<10:
                try:
                    int(x)
                except:
                    continue
                pmidstyle.append(x)
    return pmidstyle

def bibadd(thisdb,thisentry):
    try:
        thisdb[thisentry]
    except:
        thisdb.entries.append(thisentry)

def replace_id_with_pmid(thistext, thisid, thispmid):
    ''' 
        return thistext after replacing all occurences of thisid with thispmid
    '''
    #print 'replacing {} with {}'.format(thisid, thispmid)
    pmidtext = u'PMID:{}'.format(thispmid)
    latexblocks = get_latex_citation_blocks(thistext)
    for block in latexblocks:
        newblock = block[block.find('{')+1:block.rfind('}')]
        newblock = '['+newblock+']'
        newblock = newblock.replace(thisid, pmidtext)
        thistext = thistext.replace('\\'+block, newblock)
        thistext = thistext.replace(block, newblock)
    for block in get_md_citation_blocks(thistext):
        newblock = '['+block[1:-1]+']'
        newblock = newblock.replace(u'@'+thisid, pmidtext)
        newblock = newblock.replace(thisid, pmidtext)
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
    output = Entrez.efetch(db="pubmed", id=pmid, rettype="medline",retmode="text")
    details = Medline.parse(output)
    thisarticle = details.next()
    return thisarticle

def get_doi_from_pmid(pmid):
    articledetails = get_article_from_pmid(pmid)
    doi_list = [string.strip(x.replace('[doi]','')) for x in articledetails['AID'] if x.endswith('[doi]')]
    if len(doi_list)>0:
        return doi_list[0]
    else:
        raise ValueError('This pmid ({}) was not found.'.format(pmid))




















































