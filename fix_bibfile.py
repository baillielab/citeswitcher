#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

'''

python fix_citations.py -f test/eme1.md -o pmid
'''

#-------------------
import string,os,sys
import io
import bibtexparser
import subprocess
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
#-------------------
cwd = os.getcwd()
sys.path.append(cwd)
import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
# - essential
parser.add_argument('-f', '--filepath',    help='filepath', default=config['default_bibfile'])
parser.add_argument('-u', '--updatebibfile',    help='bibfile', default=config['default_updatebibfile'])
args = parser.parse_args()
#-------------------
# read bibtex file
with open(args.bibfile) as bibtex_file:
    bibdat = bibtexparser.bparser.BibTexParser(common_strings=True, homogenize_fields=True).parse_file(bibtex_file)
print ("read bib file")
#-------------------
# make dictionaries
pmids = {}
ids = {}
for entry in bibdat.entries:
    try:
        ids[entry['ID']]
        print(("duplicate ID in biblatex database:", entry["ID"]))
        if 'pmid' in entry:
            # replace this entry with a new one that has a PMID
            ids[entry['ID']]=entry
    except:
        ids[entry['ID']]=entry
    try:
        entry["pmid"]
    except:
        continue
    try:
        pmids[entry["pmid"]]
        print(("duplicate PMID in biblatex database:", entry["pmid"]))
    except:
        pass
    pmids[entry["pmid"]] = entry
#-------------------
# find all citations
idcitations = citefunctions.get_md_citations(text)
idcitations += citefunctions.get_latex_citations(text)
pmidcitations = citefunctions.get_pmid_citations(text)
print(idcitations)
print(pmidcitations)
#-------------------
# make new output database
db = BibDatabase()
for thisid in idcitations:
    print(thisid)
    try:
        ids[thisid]
    except:
        bestmatchingkey = citefunctions.find_similar_keys(thisid, ids)
        print(("bibtex id not found in biblatex file: {}. Best match in database is {}.".format(thisid, bestmatchingkey)))
        continue
    citefunctions.bibadd(db,ids[thisid])
# search through all the cited papers in this file and add them to the db 
update = BibDatabase()
for thispmid in pmidcitations:
    try:
        pmids[thispmid]
    except:
        b = citefunctions.p2b(thispmid)
        if b != 'null' and b != None:
            pmids[thispmid] = b
            citefunctions.bibadd(update,pmids[thispmid])
            print ("PMID:{} found online".format(thispmid))
        else:
            print ("PMID:{} NOT FOUND ON PUBMED".format(thispmid))
            continue
    citefunctions.bibadd(db,pmids[thispmid])
#-----------------
# update database with all the PMIDs we can find
for entry in db.entries:
    downloaded_id_list={'error':'no id to search for'}
    if 'pmid' not in entry:
        for idtype in ['doi', 'pmcid']:
            if idtype in entry:
                downloaded_id_list = citefunctions.id_translator(entry[idtype])
                if 'pmid' in downloaded_id_list:
                    entry['pmid'] = downloaded_id_list['pmid']
                    break
    if 'pmid' not in entry:
        # if no success with doi or pmcid, search pubmed
        searchstring = "{}".format(entry['title']).replace('{','').replace('}','')
        found = citefunctions.search_pubmed(searchstring)
        if len(found) == 1:
            entry['pmid'] = found[0]
    if 'pmid' in entry:
        if args.outputstyle == 'pubmed' or args.outputstyle == 'pmid':
            text = citefunctions.replace_id_with_pmid(text, thisid = entry['ID'], thispmid = entry['pmid'])
        elif args.outputstyle == 'markdown' or args.outputstyle == 'md':
            text = citefunctions.replace_pmid_with_id(text, thisid = entry['ID'], thispmid = entry['pmid'], style='md')
        elif args.outputstyle == 'latex' or args.outputstyle == 'tex':  
            text = citefunctions.replace_pmid_with_id(text, thisid = entry['ID'], thispmid = entry['pmid'], style='tex')
    else:
        print("no PMID found online for {}. {}".format(entry['ID'], downloaded_id_list['error']))

#-----------------
# save remote bibliography
with open(bibout, 'w') as bibtex_file:
    bibtexparser.dump(db, bibtex_file)
# save update bibliography
with open(args.updatebibfile, 'w') as bibtex_file:
    bibtexparser.dump(update, bibtex_file)
#-----------------
# save new text file
text = citefunctions.addheader(text, bibtex_file, args.cslfile)
with io.open(outputfile, 'w', encoding='utf-8') as file:
    file.write(text)
#-----------------
# run pandoc
if args.pandoc:
    cmd = "python run_pandoc.py -f {}".format(outputfile)
    subprocess.call(cmd, shell=True)
























