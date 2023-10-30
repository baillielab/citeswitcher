#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
#version 0.1

'''
read through a dir containing md files
find doi, pmid, in yaml
put doi, pmid in yaml if missing
put date in yaml if missing
get the front page of a pdf file if there is one in zotero 
'''

#-------------------
import os
import io
import sys
import oyaml as yaml
import json
import copy
import datetime
import subprocess
from PyPDF2 import PdfWriter,PdfReader
#-------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath, 'python-bibtexparser-master/'))
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
#-------------------
import bib
bib.init()
import include
import wordcount
import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir', default="./", help='working directory')
parser.add_argument('-i', '--imagedir', default="../img/paperpics/", help='image dir')
parser.add_argument('-v', '--verbose',    action="store_true", default=False,    help='increases verbosity')
args = parser.parse_args()
#-------------------
config['default_bibfile'] = os.path.abspath(os.path.expanduser(config['default_bibfile']))
#-------------------
print (config['default_bibfile'])
bib.full_bibdat = citefunctions.read_bib_files([config['default_bibfile']])
bib.make_alt_dicts()
#-------------------

def get_date(thisid):
    try:
        bib.full_bibdat.entries_dict[thisid]
    except:
        print (thisid, "not found in bib database")
        return ""
    # sort weird uppercase bib labels
    if "Year" in bib.full_bibdat.entries_dict[thisid]:
        bib.full_bibdat.entries_dict[thisid] =  {k.lower(): v for k, v in bib.full_bibdat.entries_dict[thisid].items()}
    # DATE
    # first handle a bizarre string format introduced sometimes by bparser:
    if "month" in bib.full_bibdat.entries_dict[thisid]:
        if isinstance(bib.full_bibdat.entries_dict[thisid]['month'], bibtexparser.bibdatabase.BibDataStringExpression):
            bib.full_bibdat.entries_dict[thisid]['month'] = bib.full_bibdat.entries_dict[thisid]['month'].get_value()
        thisdate = '{} {}'.format(bib.full_bibdat.entries_dict[thisid]['month'][:3], bib.full_bibdat.entries_dict[thisid]['year'])
        outdate = datetime.datetime.strptime(thisdate, '%b %Y')
    else:
        outdate = datetime.datetime.strptime("{}".format(bib.full_bibdat.entries_dict[thisid]['year']), '%Y')
    return outdate

def get_front_page(thisid, imagedir): 
    try:
        bib.full_bibdat.entries_dict[thisid]
    except:
        print (thisid, "not found in bib database")
        return False
    pdf = os.path.join(imagedir,"{}.pdf".format(thisid))
    pdf_page = os.path.join(imagedir,"{}.page.pdf".format(thisid))
    png_page = os.path.join(imagedir,"{}.png".format(thisid))
    filefound = False
    if not os.path.exists(pdf):
        if "file" in bib.full_bibdat.entries_dict[thisid]:
            files = bib.full_bibdat.entries_dict[thisid]["file"].split(";")
            for f in files:
                if f.endswith(".pdf"):
                    if os.path.exists(f):
                        cmd = 'cp "{}" "{}"'.format(f, pdf)
                        subprocess.call(cmd, shell=True)
                        filefound = True
    if not filefound:
        print ("No PDF file found for {}".format(thisid))
        return False
    else:
        if not os.path.exists(pdf_page):
            inputpdf = PdfReader(open(pdf, "rb"))
            output = PdfWriter()
            output.add_page(inputpdf.pages[0])
            with open(pdf_page, "wb") as outputStream:
                output.write(outputStream)
        if not os.path.exists(png_page):
            cmd = 'gm convert -density 100 "{}" "{}"'.format(pdf_page, png_page)
            subprocess.call(cmd, shell=True)
        return True

def checkdic(key, thisdict):
    '''
    return True if the dict contains a string with len>0 at this key
    '''
    if key in thisdict:
        if thisdict[key] is not None:
            if isinstance(thisdict[key], str):
                if len(thisdict[key].strip())>0:
                    return True
            else:
                return True
    return False

#-------------------
args.dir = os.path.abspath(args.dir)
ignorelist = ["__README.md"]
files = [x for x in os.listdir(args.dir) if not x.startswith(".") and not x in ignorelist]

#files = ["generalisable_stratification.md"]

for i,thisfile in enumerate(files):
    filepath = os.path.join(args.dir, thisfile)
    if os.path.isdir(filepath):
        continue
    with io.open(filepath, "r", encoding="utf-8") as f:
        text = f.read()
    text = citefunctions.make_unicode(text)
    h,r = citefunctions.readheader(text)
    y = citefunctions.getyaml(h)
    thisid = None
    newyaml = {}
    if "id" in y:
        if checkdic("id", y):
            thisid = y["id"]
            try:
                thisone = bib.full_bibdat.entries_dict[thisid]
            except:
                print ("FAILED to find {}\n... in this bibfile: {}".format(thisid, config["default_bibfile"]))
    if thisid is None and "pmid" in y:
        if checkdic("pmid", y):
            thisone = citefunctions.findcitation(str(y["pmid"]), "pmid")
            thisid = thisone["ID"]
    if thisid is None and "doi" in y:
        if checkdic("doi", y):
            thisone = citefunctions.findcitation(str(y["doi"]), "doi")
            thisid = thisone["ID"]
    if thisid is not None:
        newyaml["id"] = thisid
        if "pmid" in thisone:
            newyaml["pmid"] = thisone["pmid"]
        if "doi" in thisone:
            newyaml["doi"] = thisone["doi"]
        if "journal" in thisone:
            newyaml["journal"] = thisone["journal"]
        if not checkdic("date", y):
            newdate = get_date(thisid)
            datestring = ('{:%Y-%m-%d}'.format(newdate))
            newyaml["date"] = datestring
        # make image file
        if not checkdic("paperpic",y):
            imagedirectory = os.path.abspath(os.path.join(args.dir, args.imagedir))
            if not os.path.exists(imagedirectory):
                print ("ERROR - image directory not found here: {}".format(imagedirectory))
            if get_front_page(thisid, imagedirectory):
                newyaml["paperpic"] = "{}.png".format(thisid)

    if len(newyaml)>0:
        y = citefunctions.mergeyaml(y,newyaml)
        r = r.replace("\n\n","\n")
        newcontents = "---\n{}---\n\n{}".format(yaml.dump(y).replace("\n\n","\n"),r)
        with open(filepath,"w") as o:
            o.write(newcontents)























