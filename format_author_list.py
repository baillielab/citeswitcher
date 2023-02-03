import os, sys
import string
import select
import difflib
import numpy as np
import pandas as pd
from collections import OrderedDict

# input format
# <name> \t <tab-separated special labels> \t <affiliation1> \t <affiliation2> ... ...
# author_section and author_subsection determine which list this author's details will be added to. For a single list, these should all be the same. 
# "main" has a special meaning and should be used to make sure a contributions file is generated.


#-----------------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--authorsource', default="", help='xlsx or csv')
parser.add_argument('-n', '--numstartaffiliation', default=0, type=int, help='number to start counting from')
parser.add_argument('-o', '--outputdir', default="autofiles", help='dir')
parser.add_argument('-m', '--small_affiliations',    action="store_false", default=True,    help='use scriptsize for affiliations')
parser.add_argument('-fa', '--fix_affiliations',    action="store_true", default=False,    help='manual fix of duplicate affiliations')
parser.add_argument('-c', '--contfile', default="contributions.md", help='filename')
parser.add_argument('-l', '--labelcol', default="labels", help='column header for labels')
parser.add_argument('-s', '--similarity', default=0.7, type=float, help='Similarity metric for catching duplicate affiliations. 0-1, 1 being a perfect match.')
args = parser.parse_args()
#-----------------------------

args.authorsource = os.path.abspath(os.path.expanduser(args.authorsource))

outputdir = args.outputdir
if args.outputdir == "autofiles": #default
    if not os.path.exists(args.outputdir):
        outputdir = os.path.join(os.path.split(args.authorsource)[0],"autofiles")

symbols = [
    "&Dagger;",
    "\*",
    "&dagger;",
    "&sect;",
    "&Vert;",
    "&para;",
    ] + list(string.ascii_lowercase) + list(string.ascii_uppercase)

nonmatchingfilename = "nonmatching_affiliations.csv"
alltextfilename = "alltext.md"

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

def getlabels(theselabels):
    ''' split a comma separated string and collapse into a single list'''
    return list(set([x.strip() for x in theselabels.split(',') if len(x)>0]))

# READ AUTHOR SOURCE FILE
if args.authorsource.endswith(".csv"):
    df = pd.read_csv(args.authorsource)
elif args.authorsource.endswith(".xlsx"):
    df = pd.read_excel(args.authorsource, dtype=str)
df = df.replace(np.nan, '', regex=True)
df = df.applymap(str) # make sure everything is a string.
for i in df.columns: # strip whitespace
    if df[i].dtype == 'object':
        df[i] = df[i].map(str.strip)
print (df)
if "Name" not in df.columns:
    df["Name"] = df["firstnames_or_initials"]+ " " + df["surname"]
affiliationcols = [x for x in df.columns if x.startswith("affiliation")]
nonaffiliationcols = [x for x in df.columns if not x.startswith("affiliation")]
for affcol in affiliationcols:
    df[affcol] = df[affcol].str.rstrip('.') # remove trailing '.'
df[args.labelcol] = df[args.labelcol].str.strip()
sourcelabels = [getlabels(x) for x in df[args.labelcol].dropna() if len(x)>0]
sourcelabels = pd.Series([x.strip() for sublist in sourcelabels for x in sublist]).unique()
if len(sourcelabels)>len(symbols):
    print ("Not enough symbols: {} needed, only {} specified".format(len(sourcelabels), len(symbols)))
labeldict = {x:symbols[i] for i,x in enumerate(sourcelabels)}
print (sourcelabels)
print (labeldict)

affiliations = []
for index, row in df.iterrows():
    for aff in row.drop(nonaffiliationcols):
        if aff not in affiliations:
            if check_affiliation(aff):
                affiliations.append(aff)

duplicates = {}
no_user_response_count = 0
for i,a in enumerate(affiliations):
    m = difflib.get_close_matches(a, affiliations, 6, args.similarity)[1:]
    if len(m)>0:
        if a not in duplicates.keys():
            print (a,duplicates.keys())
            m = [a]+m
            print ("\n*** Possible duplicate affiliation:\n{}".format("\n".join(["{}. {}".format(i,x) for i,x in enumerate(m)])))
            if args.fix_affiliations:
                i,o,e = select.select([sys.stdin],[],[],20) # 10 second timeout
                if no_user_response_count > 5: # user not interested. give up
                    continue
                if i==[] and o==[] and e==[]:
                    no_user_response_count += 1
                if i:
                    q = sys.stdin.readline().strip()
                    q = q.strip().upper()
                    try:
                        q=int(q)
                    except:
                        print ("Not an integer")
                        continue
                    try:
                        m[q]
                    except:
                        print ("{} not in list".format(q))
                        continue
                    print ("Choosing {}".format(m[q]))
                    new = m[q]
                    for x in m:
                        duplicates[x] = new
                    continue
            # if not replacing, keep everything the same
            for x in m:
                print (x)
                duplicates[x] = x

if args.fix_affiliations:
    df.replace(duplicates, inplace=True)
    if args.authorsource.endswith(".csv"):
        df.to_csv(args.authorsource, index=False)
    elif args.authorsource.endswith(".xlsx"):
        df.to_excel(args.authorsource, dtype=str)

dn = df[["Name"]+affiliationcols]
dn = dn.sort_values(by="Name")
dn = dn.drop_duplicates()
dn.to_csv(os.path.join(outputdir,nonmatchingfilename))

# MAKE A SINGLE LIST OF CONTRIBUTORS FOR EACH SECTION
authsecs = [x for x in df["author_section"].dropna().unique() if len(x)>0]
print (authsecs)
alltext = ""
for author_section in authsecs:
    these_authors = []
    these_affiliations = []
    these_labels = []
    outputfilename = ''.join(e for e in author_section if e.isalnum()).lower() + ".md"
    s = df.loc[df['author_section'] == author_section]
    subs = {x:[] for x in s["author_subsection"].dropna().unique()} #NB includes empty
    print (subs)
    for index, row in s.iterrows():
        thisname = row["Name"]
        superscript = []
        this_author_aff = []
        if row[args.labelcol] != "":
            lab = sorted(getlabels(row[args.labelcol]))
            superscript += [labeldict[x] for x in lab]
            these_labels += lab
        for aff in row.drop(nonaffiliationcols):
            if check_affiliation(aff):
                this_author_aff.append(affiliations.index(aff)+1+args.numstartaffiliation)
                these_affiliations.append(aff)
        thisname = "{}".format(thisname)
        superscript += [str(x) for x in sorted(this_author_aff)]
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
        subtext = ""
        if author_section != "main":
            subtext += "#### {}\n\n".format(author_section)
        for sub in subs:
            if sub != "main":
                if len(sub.strip())>0:
                    subtext += "##### {}\n\n".format(sub)
            subtext += ",\n".join(subs[sub])
            subtext += ".  \n\n"
        if args.small_affiliations:
            subtext += "\\scriptsize\n\n"
        for label in these_labels:
            subtext += "{} - {}  \n\n".format(labeldict[label], label)
        subtext += "  \n\n"
        for i,aff in enumerate(these_affiliations):
            subtext += "^{}^ {}  \n".format(affiliations.index(aff)+1+args.numstartaffiliation, aff)
        subtext += "\n\n"
        if args.small_affiliations:
            subtext += "\\normalsize\n\n"
        o.write(subtext)
        alltext += subtext

with open(os.path.join(outputdir, alltextfilename),"w") as o:
    o.write(alltext)

# Write author contributions statement
main = df.loc[df["author_section"] == "main"]
if len(main)>0:
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


print (len(df["Name"].unique()), "unique people named")




















