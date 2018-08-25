#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

from bibtexparser.bibdatabase import BibDatabase
import latexchars

def init():
    global db
    db = BibDatabase()
    global pmids
    pmids = {}
    global ids
    ids = {}
    global newcitations
    newcitations = []

def make_dictionaries(bib_db):
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

def new(entry):
    '''
        add new reference to pmids, ids and bib db
    '''
    if entry is not None:
        try:
            entry['ENTRYTYPE']
        except:
            entry['ENTRYTYPE'] = 'article'
        try:
            ids[entry['ID']]
        except:
            entry = latexchars.cleanbib(entry)
            db.entries = [entry] + db.entries
            ids[entry['ID']] = entry
            newcitations.append(entry['ID'])
            try:
                entry["pmid"]
            except:
                return 0
            try:
                pmids[entry["pmid"]]
            except:
                pmids[entry["pmid"]] = entry


def supplement(theseids):
    for thisid in theseids:
        try:
            ids[thisid]
        except:
            continue
        db.entries = [ids[thisid]] + db.entries




