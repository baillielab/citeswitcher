﻿#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

import string,os,sys
#-------------------
cwd = os.getcwd()
sys.path.append(cwd)
import citefunctions
config = citefunctions.getconfig(os.path.join(os.path.dirname(__file__), 'config.json'))
#-------------------
alltests = False
#-------------------

def test_translator():
    print(citefunctions.id_translator('10.1093/nar/gks1195'))
    print(citefunctions.id_translator('23193287'))
    print(citefunctions.id_translator('adsfsddd'))

text0 = 'this is some ´ ` ´´®†©¨¥© text with a in [@bibtexid] in it  (‘A2B trial’).'
text1 = 'The safety and efficacy of medication practice \cite{wang} is at least partly determined genetically[PMID:23423423 PMID 23421432; PMID :234231424].  CPIC dosing guidelines[\\cite{relling_cpic_2011}], based on pretreatment testing for multiple single gene-drug interactions, exist for 100 drugs [\\cite{noauthor_pharmgkb_nodate}, #cite{whirl-carrillo_pharmacogenomics_2012}] and there is RCT evidence for a limited number of guidelines [\\cite{price_pharmacogenomic_2013}].  There are identified alleles for both dexmedetomidine [\\cite{posti_polymorphism_2013, yagar_role_2011, kurnik_genetic_2011}] and clonidine [@nurnberger_effect_2003; @yang_association_2010] which affect the metabolism or molecular pathways of the drugs, and which are associated with differential treatment efficacy (DTE) in groups of patients which are identifiable pretreatment.'
text2 = 'In addition to the sedative effect for which illness[@KimEffectsdexmedetomidineratio2014; @Sedationimprovesearly2009], which may mediate'
text3 = 'In addition to direct sedative effects on the central nervous system, alpha_2-agonistsmodulate [@bibtexid; PMID:1234xx1243]. We hypothesise that detectable, modifiable inflammator signals in peripheral blood cause delirium in some patients[@bibtexid] Also \\cite{bibtexid}.'


def test_replacement(thistext):
    bibdat = citefunctions.read_bib_file('test/formatted/eme.bib')
    pmids, ids = citefunctions.make_dictionaries(bibdat.entries)
    print(thistext)
    print(citefunctions.replace_blocks(thistext,pmids, ids, 'pmid'))

def test_mdcite():
    print(citefunctions.get_md_citations(text1))

def test_latexcite():
    print(citefunctions.get_latex_citations(text2))

def test_stringmatch():
    d = {
        "abcde":1,
        "defgs":1,
        "wange":1,
        "ewang":1,
        "frmamd":1,
        "bibblebob":1,
    }
    citefunctions.find_similar_keys('abcde', d)

def test_p2b(thispmid):
    #new_article = citefunctions.get_article_from_pmid(thispmid)
    print(citefunctions.p2b(thispmid))

def test_readhead(thisfile):
    f=open(thisfile)
    text = f.read()
    f.close()
    print (citefunctions.readheader(text))

def test_addhead(thisfile):
    f=open(thisfile)
    text = f.read()
    f.close()
    newfile =  (citefunctions.addheader(text, 'test/eme.bib', 'sup-files/csl/nature.csl'))
    print ("new:")
    print (newfile)

def test_replaceblocks():
    pmids={}
    ids={}
    print (citefunctions.replace_blocks(text1, {}, {}, 'md'))

def test_par():
    thistext = 'asdf [ asf aas d as [asdf ]asf asdf asdf] ads fa sdf'
    print (thistext)
    print ("EEE:", citefunctions.get_parethesised(thistext))

test_par()
#test_replaceblocks()
#test_addhead('test/eme.txt')
#test_readhead('test/eme.txt')
#test_p2b('29494619')
#test_p2b(29494619)
#test_p2b('xxxx')
#test_translator()
#test_replacement(text2)
#test_mdcite()
#test_latexcite()
#test_translator()
#test_translator()
#test_translator()
#test_stringmatch()