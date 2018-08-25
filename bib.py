#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

from bibtexparser.bibdatabase import BibDatabase
import latexchars

def init():
    global full_bibdat
    full_bibdat = BibDatabase()
    global db # cited 
    db = BibDatabase()
    global pmids
    pmids = {}

def make_dictionaries():
    for entry in full_bibdat.entries:
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

def new(entry):
    '''
        add new reference to pmids, db.entries_dict and bib db
    '''
    if entry is not None:
        try:
            entry['ENTRYTYPE']
        except:
            entry['ENTRYTYPE'] = 'article'
        try:
            db.entries_dict[entry['ID']]
        except:
            entry = latexchars.cleanbib(entry)
            db.entries = [entry] + db.entries
            print ("testing", db.entries_dict[entry['ID']])
            try:
                entry["pmid"]
            except:
                return 0
            try:
                pmids[entry["pmid"]]
            except:
                pmids[entry["pmid"]] = entry


def supplement(theseids):
    print (theseids)
    for thisid in theseids:
        try:
            db.entries_dict[thisid]
        except:
            continue
        db.entries = [ids[thisid]] + db.entries




