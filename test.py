#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

import string,os,sys
#-------------------
cwd = os.getcwd()
sys.path.append(cwd)
import citefunctions
#-------------------
alltests = False
#-------------------

def test_translator():
    print citefunctions.id_translator('10.1093/nar/gks1195')
    print citefunctions.id_translator('23193287')
    print citefunctions.id_translator('adsfsddd')

text = u'this is some ´ ` ´´®†©¨¥© text with a in [@bibtexid] in it  (‘A2B trial’).'
text1 = u'The safety and efficacy of medication practice is at least partly determined genetically.  CPIC dosing guidelines[\\cite{relling_cpic_2011}], based on pretreatment testing for multiple single gene-drug interactions, exist for 100 drugs [\\cite{noauthor_pharmgkb_nodate}, #cite{whirl-carrillo_pharmacogenomics_2012}] and there is RCT evidence for a limited number of guidelines [\\cite{price_pharmacogenomic_2013}].  There are identified alleles for both dexmedetomidine [\\cite{posti_polymorphism_2013, yagar_role_2011, kurnik_genetic_2011}] and clonidine [\\cite{nurnberger_effect_2003, yang_association_2010}] which affect the metabolism or molecular pathways of the drugs, and which are associated with differential treatment efficacy (DTE) in groups of patients which are identifiable pretreatment.  Estimates of this DTE are imprecise.  Realising the benefit of genomic testing relies on reducing cost of t'
text2 = u'In addition to the sedative effect for which they are usually prescribed, alpha_2-agonists have extensive immune-modulating effects in animal models of critical illness[@KimEffectsdexmedetomidineratio2014; @Sedationimprovesearly2009], which may mediate'

def test_replacement():
    btid = u'bibtexid'
    pmid = u'1234xx1243'
    print text1
    text1 = citefunctions.replace_id_with_pmid(text1, 'relling_cpic_2011', pmid)
    print text1



def test_mdcite():
    print citefunctions.get_md_citations(text1)

def test_latexcite():
    print citefunctions.get_latex_citations(text2)

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




#test_translator()
#test_replacement()
#test_mdcite()
#test_latexcite()
#test_translator()
#test_translator()
#test_translator()
test_stringmatch()