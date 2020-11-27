import os, sys
import difflib
import numpy as np
import pandas as pd
from collections import OrderedDict

# input format
# <name> \t <comma-separated special labels> \t <affiliation1> \t <affiliation2> ... ...

#-----------------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--authorsource', default="", help='xlsx or csv')
parser.add_argument('-o', '--outputdir', default="autofiles", help='dir')
parser.add_argument('-c', '--contfile', default="contributions.md", help='filename')
parser.add_argument('-l', '--labelcol', default="labels", help='column header for labels')
parser.add_argument('-s', '--similarity', default=0.7, type=float, help='Similarity metric for catching duplicate affiliations. 0-1, 1 being a perfect match.')
args = parser.parse_args()
#-----------------------------

args.authorsource = os.path.abspath(os.path.expanduser(args.authorsource))

if args.outputdir == "autofiles": #default
    outputdir = os.path.join(os.path.split(args.authorsource)[0],"autofiles")
else:
    outputdir = args.outputdir

symbols = [
    "&Dagger;",
    "\*",
    "&dagger;",
    "&sect;",
    "&Vert;",
    "&para;",
    "a",
    "b",
    "c",
    "d",
    "e",
    "f",
    "g",
    "h",
    ]

def check_affiliation(thisaffiliation):
    return len(thisaffiliation.strip())>0

def get_initials(name, surnamelen=1):
    these_initials = ""
    namewords = name.split(" ")
    for i,word in enumerate(namewords):
        thislen = 1
        if i==len(namewords)-1:
            thislen = surnamelen
        if "-" in word:
            a = "-".join([x[:thislen] for x in word.split("-")])
        else:
            a = word[:thislen]
        these_initials+=a
    return these_initials

# READ AUTHOR SOURCE FILE
if args.authorsource.endswith(".csv"):
    df = pd.read_csv(args.authorsource, sep='\t')
elif args.authorsource.endswith(".xlsx"):
    df = pd.read_excel(args.authorsource, dtype=str)
df = df.replace(np.nan, '', regex=True)
df = df.applymap(str) # make sure everything is a string.
print (df)
affiliationcols = [x for x in df.columns if x.startswith("affiliation")]
nonaffiliationcols = [x for x in df.columns if not x.startswith("affiliation")]
for affcol in affiliationcols:
    df[affcol] = df[affcol].str.rstrip('.') # remove trailing '.'
sourcelabels = [x for x in df["labels"].dropna().unique() if len(x)>0]
if len(sourcelabels)>len(symbols):
    print ("Not enough symbols: {} needed, only {} specified".format(len(sourcelabels), len(symbols)))
labeldict = {x:symbols[i] for i,x in enumerate(sourcelabels)}
affiliations = []
for index, row in df.iterrows():
    for aff in row.drop(nonaffiliationcols):
        if aff not in affiliations:
            if check_affiliation(aff):
                affiliations.append(aff)

duplicates = []
for i,a in enumerate(affiliations):
    m = difflib.get_close_matches(a, affiliations, 6, args.similarity)[1:]
    if len(m)>0:
        if a not in duplicates:
            print ("\n*** Possible duplicate affiliation:\n{}\n{}".format(a,"\n".join(m)))
        duplicates += [a] + m

dn = df[["Name"]+affiliationcols]
dn = dn.sort_values(by="Name")
dn = dn.drop_duplicates()
dn.to_csv(os.path.join(outputdir,"nonmatching_affiliations.csv"))

# MAKE A SINGLE LIST OF CONTRIBUTORS FOR EACH SECTION
authsecs = [x for x in df["author_section"].dropna().unique() if len(x)>0]
for author_section in authsecs:
    these_authors = []
    these_affiliations = []
    these_labels = []
    outputfilename = ''.join(e for e in author_section if e.isalnum()).lower() + ".md"
    s = df.loc[df['author_section'] == author_section]
    subs = {x:[] for x in s["author_subsection"].dropna().unique()} #NB includes empty
    for index, row in s.iterrows():
        thisname = row["Name"]
        superscript = []
        if row["labels"] != "":
            lab = sorted(list(set([x for x in row["labels"].split(',') if len(x)>0])))
            superscript += [labeldict[x] for x in lab]
            these_labels += lab
        for aff in row.drop(nonaffiliationcols):
            if check_affiliation(aff):
                superscript.append(str(affiliations.index(aff)+1))
                these_affiliations.append(aff)
        thisname = "{}".format(thisname)
        if len(superscript)>0:
            thisname += "^{}^".format(",".join(superscript))
        #print(row["affiliation1"])
        these_authors.append(thisname)
        subs[row["author_subsection"]].append(thisname)

    # get unique values in order
    these_labels = list(OrderedDict.fromkeys(these_labels))
    these_affiliations = list(OrderedDict.fromkeys(these_affiliations))
    #these_affiliations = list(OrderedDict.fromkeys(affiliations)) # one global affiliation record

    with open(os.path.join(outputdir, outputfilename),"w") as o:
        #o.write("#### {}\n\n".format(author_section))
        for sub in subs:
            if len(sub.strip())>0:
                o.write("##### {}\n\n".format(sub))
            o.write(",\n".join(subs[sub]))
            o.write(".\n\n")
        for label in these_labels:
            o.write("{} - {} \n\n".format(labeldict[label], label))
        o.write("\n\n")
        for i,aff in enumerate(these_affiliations):
            o.write("^{}^{}  \n".format(affiliations.index(aff)+1, aff))
        o.write("\n\n")

# Write author contributions statement
main = df.loc[df["author_section"] == "Main author list"]
names = main["Name"].to_list()
initials = {}
already = []
for name in names:
    z = 1
    ini = get_initials(name, z)
    while z<4 and ini in already:
        z+=1
        ini = get_initials(name, z)
    if ini in already:
        print ("duplicate not resolved! {} {}".format(name, ini))
    already.append(ini)
    initials[name] = ini
cont = {}
for index, row in main.iterrows():
    ini = initials[row["Name"]]
    contributions = [x.replace(".","").strip() for x in row["contribution"].split(",")]
    for c in contributions:
        if len(c)==0:
            continue
        try:
            cont[c]
        except:
            cont[c]=[]
        cont[c].append(ini)
with open(os.path.join(outputdir, args.contfile),"w") as o:
    for c in cont:
        if len(c)>1:
            o.write("{} and {} contributed to {}.\n".format(
                ', '.join(cont[c][:-1]),
                cont[c][-1],
                c
                ))
        else:
            o.write("{} contributed to {}.\n".format(
                cont[c][0],
                c
                ))


print (len(df["Name"].unique()), "unique names")




















