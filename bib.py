#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

from bibtexparser.bibdatabase import BibDatabase
from bibtexparser.latexenc import string_to_latex

def init():
    global full_bibdat
    full_bibdat = BibDatabase()
    global db # cited 
    db = BibDatabase()
    global pmids
    pmids = {}
    global dois
    dois = {}

def make_alt_dicts():
    thesedics = [
            (pmids, "pmid"),
            (dois, "doi"),
        ]
    for entry in full_bibdat.entries:
        for thisdict, thislabel in thesedics:
            try:
                entry[thislabel]
            except:
                continue
            try:
                thisdict[entry[thislabel]]
                print("duplicate {} in biblatex database:{}".format(thislabel, entry[thislabel]))
            except:
                pass
            thisdict[entry[thislabel]] = entry

def new(entry):
    '''
        add new reference to pmids, db.entries, db.entries_dict
    '''
    if entry is not None:
        try:
            entry['ENTRYTYPE']
        except:
            entry['ENTRYTYPE'] = 'article'
        try:
            db.entries_dict[entry['ID']]
        except:
            entry = cleanbib(entry)
            db.entries = [entry] + db.entries
            db.entries_dict[entry['ID']] = entry
            try:
                entry["pmid"]
            except:
                return 0
            try:
                pmids[entry["pmid"]]
            except:
                pmids[entry["pmid"]] = entry

def supplement(theseids):
    #print ("supp:", theseids)
    for thisid in theseids:
        try:
            full_bibdat.entries_dict[thisid]
        except:
            continue
        db.entries = [full_bibdat.entries_dict[thisid]] + db.entries
        db.entries_dict[thisid] = full_bibdat.entries_dict[thisid]

def cleanbib(bibtex_entry):
    return bibtex_entry # this function disabled as it may not be necessary now because p2b() always returns clean latex
    #return {d:string_to_latex(bibtex_entry[d]) for d in bibtex_entry.keys()}





