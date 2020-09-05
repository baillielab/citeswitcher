#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

'''
    read through a bib file and add pmids for every citation that has one
    store the work in a json file so it isn't repeated next time we run
    can be run as a cron job to keep pmids up to date


    pyparsing error is corrected by downloading latest bibtexparser
https://stackoverflow.com/questions/49321998/bibtexparser-pyparsing-parseexception-expected-end-of-text
'''

#-------------------
import os
import io
import sys
import json
import subprocess
#-------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath, 'python-bibtexparser-master/'))
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
#-------------------
import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bibfile',    help='bibfile', default=config['zotero_bibfile'])
parser.add_argument('-o', '--outbib',    help='outbib', default=config['default_bibfile'])
args = parser.parse_args()
#-------------------
bibfile = os.path.expanduser(args.bibfile)
outbib = os.path.expanduser(args.outbib)
#-------------------
speedup_store = outbib.replace('.bib','.json')
bibdat = citefunctions.read_bib_files([bibfile])
#-------------------
if os.path.exists(speedup_store):
    with open(speedup_store) as f:
        try:
            already = json.load(f)
        except:
            already = {}
else:
    already = {}
#-------------------
c=0
d=0
e=0
for thisentry in bibdat.entries:
    try:
        # always look in the existing bib db first
        thisentry['pmid']
        c+=1
    except:
        try:
            thisentry['doi']
            d+=1
        except:
            #print ("No doi: {}".format(thisentry['ID']))
            continue
        try:
            already[thisentry['doi']]
        except:
            # only search if this hasn't already been done
            pub = citefunctions.search_pubmed(thisentry['doi'], "doi")
            if len(pub) == 1:
                pmid = pub[0]
                thisentry['pmid'] = pmid
                e+=1
                already[thisentry['doi']] = pmid
            else:
                already[thisentry['doi']] = 'failed'
                print ("failed {} {} {}".format(thisentry['ID'], thisentry['doi'], pub))
print ("\nPMIDs already:{}, DOIs:{}, PMIDs found:{}".format(c,d,e))
#-----------------
# save details to speed up next run
with open(speedup_store,'w') as o:
    json.dump(already,o)
#-----------------
# save update bibliography
with open(outbib, 'w') as bf:
    bibtexparser.dump(bibdat, bf)
#-----------------
























