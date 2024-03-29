﻿#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

import io
import os
import re
import sys
import pandas as pd
import numpy as np
from tabulate import tabulate
from requests import get  # to make GET request
#-----------------------------
'''
RULES:

html:
INCLUDESECTION IN LINE
#INCLUDESECTION NOT IN LINE
WHOLE LINE IS REPLACED
FILENAME IS NEXT WORD AFTER INCLUDESECTION

md:
OR use {!blah.txt!}
OR use {!blah.xlsx!}

{{< include _content.qmd >}}

xlsx:
    - if £ in column name, format currency
    - drop columns beginning with #

'''
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#-----------------------------
comments = r'<!--[\s\S]+?-->'
tex_include_formats = [
            r'\\input{.+?}',
            r'\\import{.+?}',
            r'\\include{.+?}'
            ] 
md_include_format = r'{!.+?!}' # support for markdown-include format: pip install markdown-include
redactionpattern = r"\[STARTREDACT\][\s\S]+?\[ENDREDACT\]"
#-----------------------------
class cd:
    """Context manager for changing the current working directory"""
    """ by Brian M. Hunt https://stackoverflow.com/users/19212/brian-m-hunt"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def newext(filepath, thisext):
    return filepath[:filepath.rfind('.')] + thisext

def preext(filepath, thisext):
    lastdot = filepath.rfind('.')
    outlist = [filepath[:lastdot], thisext, filepath[lastdot+1:]]
    return ".".join(outlist)

def get_filename(thisinclude_instruction):
    x = thisinclude_instruction.split(' ')
    return x[x.index("INCLUDESECTION")+1].strip().replace("'",'').replace('"','').replace('-->','')

def get_includes(thistext):
    includedict = {}
    for x in re.findall(comments, thistext):
        if "INCLUDESECTION" in x and not "#INCLUDESECTION" in x:
            includedict[x] = get_filename(x)
    for tex_include_format in tex_include_formats:
        for x in re.findall(tex_include_format, thistext):
            includedict[x] = os.path.expanduser(x[7:-1].strip())+".tex"
    for x in re.findall(md_include_format, thistext):
        includedict[x] = os.path.expanduser(x[2:-2].strip())
    return includedict

def clear_nan(thisdf):
    thisdf = thisdf.replace('£nan', '', regex=True)
    thisdf = thisdf.replace('nan', '')
    thisdf = thisdf.replace('NaN', '', regex=True)
    thisdf = thisdf.replace(np.nan, '', regex=True)
    return thisdf

def allintegers(this_series):
    return np.array_equal(this_series.dropna(), this_series.dropna().astype(int))

'''
# THIS FUNCTION LEFT HERE FOR NOW BECAUSE IT IS USED FOR BUDGET ELSEWHERE
def include_df(thisfile, filetype="xlsx"):
    # BUDGET
    if filetype == "xlsx":
        df = pd.read_excel(thisfile, index_col=0, dtype=str)
    elif filetype == "csv":
        df = pd.read_csv(thisfile, index_col=0)
    elif filetype == "tsv":
        df = pd.read_csv(thisfile, index_col=0, sep="\t")
    else:
        print ("{} not read as {}".format(thisfile, filetype))
        return
    print (df)
    for colname in df.columns:
        if colname.startswith('#'):
            df = df.drop(colname, axis=1)
        if "£" in colname:
            df[colname] = pd.to_numeric(df[colname])
            df[colname] = df[colname].map("£{0:,.2f}".format)
    # format number
    for thiscol in []:
        df[thiscol] = df[thiscol].map("{}".format)
    df = clear_nan(df)
    return tabulate(df, tablefmt="grid", headers="keys")
'''
def include_df(thisfile, filetype="xlsx", tf="pipe"):
    '''
    parse and include table file

    table format options
    tablefmt="pipe"
    tablefmt="grid"
    tablefmt="github"
    tablefmt="simple"
    tablefmt="fancygrid"
    tablefmt="latex"
    plain, simple, github, grid, fancy_grid, pipe,
                          orgtbl, rst, mediawiki, html, latex, latex_raw,
                          latex_booktabs, tsv
    plain, simple, github, grid, fancy_grid, pipe, orgtbl, jira, presto, 
    psql, rst, mediawiki, moinmoin, youtrack, html, latex, latex_raw, 
    latex_booktabs, textile, 
    '''
    if filetype == "xlsx":
        df = pd.read_excel(thisfile, index_col=None, dtype=str)
    elif filetype == "csv" or filetype =="tsv":
        if filetype == "tsv":
            df = pd.read_csv(thisfile, index_col=None, sep="\t")
        else:
            df = pd.read_csv(thisfile, index_col=None)
    else:
        print ("{} not read as {}".format(thisfile, filetype))
        return
    df = clear_nan(df)
    # if any colname starts with "chr:pos" then sort by it
    chrlist = [x for x in list(df.columns) if x.startswith("chr:pos")]
    if len(chrlist)>0:
        sortcol = chrlist[0]
        tempchr = "chrtempforsorting123412341234"
        temppos = "postempforsorting123412341234"
        df[[tempchr,temppos]] = df[sortcol].str.split(":", expand=True)
        df[tempchr] = pd.to_numeric(df[tempchr])
        df[temppos] = pd.to_numeric(df[temppos])
        df.sort_values(by=[tempchr,temppos], inplace=True)
        df.drop([tempchr,temppos], axis=1, inplace=True)
    for i,colname in enumerate(df.columns):
        if colname.startswith('#'):
            df = df.drop(colname, axis=1)
        if colname == "[rowcolors]":
            # set colors for latex
            for rowIndex, row in df.iterrows(): #iterate over rows
                if row["[rowcolors]"] != "":
                    for columnIndex, value in row.items():
                        if df.loc[rowIndex, columnIndex] != '':
                            df.loc[rowIndex, columnIndex] = '\color{%s}%s'%(row["[rowcolors]"], df.loc[rowIndex, columnIndex])
            df = df.drop(colname, axis=1)
        if colname == "[rowbckgdcolors]":
            # set colors for latex NB need - \usepackage{colortbl}
            for rowIndex, row in df.iterrows(): #iterate over rows
                if row["[rowbckgdcolors]"] != "":
                    for columnIndex, value in row.items():
                        if df.loc[rowIndex, columnIndex] != '':
                            df.loc[rowIndex, columnIndex] = '\cellcolor{%s}%s'%(row["[rowbckgdcolors]"], df.loc[rowIndex, columnIndex])
            df = df.drop(colname, axis=1)
    newcolnames = {x:x for x in df.columns}
    alignmentinstructions = []
    for i,colname in enumerate(df.columns):
        if "[align:" in colname:
            newcolnames[colname] = colname.split("[align:")[0]
            alignmentinstructions.append(colname.split("[align:")[1].replace("]","").strip())
        else:
            alignmentinstructions.append("")
        if "£" in colname:
            df[colname] = pd.to_numeric(df[colname])
            df[colname] = df[colname].map("£{0:,.2f}".format)
        else:
            try:
                df[colname] = pd.to_numeric(df[colname])
            except:
                continue
            # round all numeric columns to 2sf
            #print (colname, allintegers(df[colname]))
            if allintegers(df[colname]):
                df[colname] = df[colname].map("{0:.0f}".format)
            else:
                df[colname] = df[colname].map("{0:.2g}".format)
    df.rename(columns=newcolnames, inplace=True)
    # format number
    for thiscol in []:
        df[thiscol] = df[thiscol].map("{}".format)
    # replace all spaces with newline
    #df.rename(columns={x:x.replace(" ","\n") for x in df.columns}, inplace=True)
    #print (df)
    df = clear_nan(df)
    out = tabulate(df, tablefmt=tf, headers="keys", showindex=False, colalign=alignmentinstructions)
    if tf in ["github","simple","grid","pipe","fancygrid"]:
        # markdown output. Need to manually change latex formatting for italics
        ilist = re.findall(r"\\textit\{.+?\}", out)
        for imatch in ilist:
            out = out.replace(imatch, "*{}*".format(imatch[8:-1]))
        ilist += re.findall(r"\\emph\{.+?\}", out)
        for imatch in ilist:
            out = out.replace(imatch, "*{}*".format(imatch[6:-1]))
    return out

def parse_includes(thisfile, verbose=False, tbf="pipe"):
    '''
        read file and return text with includes
    '''
    if os.path.exists(thisfile):
        try:
            with io.open(thisfile, "r", encoding="utf-8") as f:
                text = f.read()
        except Exception as e:
            print ("\n\n*** INITIAL I/O PROBLEM (COULD BE DROPBOX): {} ***\n\n".format(thisfile))
            print(e)
            sys.exit()
            return ""
    else:
        print ("\n\n*** INITIAL FILE NOT FOUND: {} ***\n\n".format(thisfile))
        sys.exit()
        return ""

    filedir = os.path.split(os.path.abspath(thisfile))[0]
    with cd(filedir):
        includes = get_includes(text)
        if verbose:
            print ("fileincludes:", includes)
        additionalfiles = list(set(includes.values()))
        if verbose:
            print ("[include.py] working in: {}".format(filedir))
        for filepath in additionalfiles:
            if os.path.exists(filepath):
                if not filepath.endswith('.xlsx'): # don't try to read excel directly
                    if verbose:
                        print ("[include.py] including:", filepath)
                    newtext = parse_includes(filepath, tbf=tbf)
            elif filepath.startswith("http"):
                pass
            else:
                print ("\n\n*** INCLUDE FILE NOT FOUND: ***\n{}".format(filepath))
                print ("*** path: {} ***".format(os.path.abspath(filepath)))
                print ("*** cwd: {} ***\n".format(os.getcwd()))
                continue
            for inc in includes:
                if includes[inc] == filepath:
                    if filepath.startswith('http'):
                        newtext = get(filepath).text
                    elif filepath.endswith('.xlsx'):
                        newtext = include_df(filepath, filetype="xlsx", tf=tbf)
                    elif filepath.endswith('.csv'):
                        newtext = include_df(filepath, filetype="csv", tf=tbf)
                    elif filepath.endswith('.tsv'):
                        newtext = include_df(filepath, filetype="tsv", tf=tbf)
                    text = text.replace(inc, newtext)
    return text

def stripcomments(thistext):
    return re.sub(comments, "", thistext)

def save_new(thisfile, outputfile="auto", stripc=False, verbose=False, tableformattype = "pipe"):
    if thisfile.endswith(".tex"):
        tableformattype = "tex"
    if verbose:
        print ("tableformattype:", tableformattype)
    text = parse_includes(thisfile, verbose=verbose, tbf=tableformattype)
    if outputfile == 'auto':
        outputfile = preext(thisfile, 'inc')
    if stripc:
        text = stripcomments(text)
    with open(outputfile,'w') as o:
        o.write(text)
    return outputfile

def redact(text):
    text, n = re.subn(redactionpattern, " (redacted section) ", text)
    print ("{} redacted sections".format(n))
    text = re.sub(r"<!--[\s\S]+?-->", "", text) # remove html comments too
    return text

#-----------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-o', '--outputfile', default="auto",   help='outputfile')
    parser.add_argument('-tf', '--tableformat', default="pipe",   help='fancy_grid, github, grid, html, jira, latex, latex_booktabs, latex_raw, mediawiki, moinmoin, orgtbl, pipe, plain, presto, psql, rst, simple, textile, tsv, youtrack')
    parser.add_argument('-s', '--stripcomments', action="store_true", default=False, help='stripcomments')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='verbose')
    args = parser.parse_args()
    save_new(args.filename, outputfile=args.outputfile, stripc = args.stripcomments, verbose=args.verbose, tableformattype=args.tableformat)


















