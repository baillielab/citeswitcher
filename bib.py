#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

from bibtexparser.bibdatabase import BibDatabase
from bibtexparser.latexenc import string_to_latex

additionaldicts = []
verbose = False

def init():
    global full_bibdat
    full_bibdat = BibDatabase()
    global db # cited
    db = BibDatabase()
    global pmids
    pmids = {}
    global dois
    dois = {}
    additionaldicts.append((pmids, "pmid"))
    additionaldicts.append((dois, "doi"))

def make_alt_dicts():
    for entry in full_bibdat.entries:
        for thisdict, thislabel in additionaldicts:
            try:
                entry[thislabel]
            except:
                continue
            try:
                thisdict[entry[thislabel]]
                if verbose:
                    print("duplicate {} in biblatex database:{}".format(thislabel, entry[thislabel]))
            except:
                pass
            thisdict[entry[thislabel]] = entry

def new(entry):
    '''
        add new reference to full_bibdat.entries, full_bibdat.entries_dict and all additionaldicts
    '''
    if entry is not None:
        try:
            entry['ENTRYTYPE']
        except:
            entry['ENTRYTYPE'] = 'article'
        try:
            full_bibdat.entries_dict[entry['ID']] # existing full_bibdat entry ALWAYS takes precedence
            # add any additional fields from online by merging dictionaries
            entry = {**full_bibdat.entries_dict[entry['ID']], **entry}
        except:
            pass
        entry = cleanbib(entry)
        full_bibdat.entries = [entry] + full_bibdat.entries
        full_bibdat.entries_dict[entry['ID']] = entry
        for thisdict, thislabel in additionaldicts:
            try:
                entry[thislabel]
            except:
                continue
            try:
                thisdict[entry[thislabel]]
            except:
                thisdict[entry[thislabel]] = entry

def cite(theseids):
    for thisid in theseids:
        try:
            db.entries_dict[thisid]
        except:
            db.entries = [full_bibdat.entries_dict[thisid]] + db.entries
            db.entries_dict[thisid] = full_bibdat.entries_dict[thisid]

def cleanbib(bibtex_entry):
    return {d:string_to_latex(bibtex_entry[d]) for d in bibtex_entry.keys()}

















