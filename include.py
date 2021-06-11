#!/opt/local/bin/python
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
INCLUDESECTION IN LINE
#INCLUDESECTION NOT IN LINE
WHOLE LINE IS REPLACED
FILENAME IS NEXT WORD AFTER INCLUDESECTION
OR use {!blah.txt!}
OR use {!blah.xlsx!}

RULES:
    - if £ in column name, format currency
    - drop columns beginning with #

'''
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#-----------------------------
comments = r'<!--[\s\S]+?-->'
md_include_format = r'{![\s\S]+?!}' # support for markdown-include format: pip install markdown-include
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
    for x in re.findall(md_include_format, thistext):
        includedict[x] = os.path.expanduser(x[2:-2].strip())
    return includedict

def clear_nan(thisdf):
    thisdf = thisdf.replace('£nan', '', regex=True)
    thisdf = thisdf.replace('nan', '')
    thisdf = thisdf.replace('NaN', '', regex=True)
    thisdf = thisdf.replace(np.nan, '', regex=True)
    return thisdf

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
def include_df(thisfile, filetype="xlsx"):
    '''
    parse and include excel file

    table format options
    tablefmt="github"
    tablefmt="simple"
    tablefmt="grid"
    tablefmt="pipe"
    tablefmt="fancygrid"
    tablefmt="latex"
    plain, simple, github, grid, fancy_grid, pipe,
                          orgtbl, rst, mediawiki, html, latex, latex_raw,
                          latex_booktabs, tsv
    '''

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
    for colname in df.columns:
        if colname.startswith('#'):
            df = df.drop(colname, axis=1)
        elif "£" in colname:
            df[colname] = pd.to_numeric(df[colname])
            df[colname] = df[colname].map("£{0:,.2f}".format)
        else:
            try:
                df[colname] = pd.to_numeric(df[colname])
            except:
                continue
            # round all numeric columns to 2sf
            df[colname] = df[colname].map("{0:.2g}".format)
    # format number
    for thiscol in []:
        df[thiscol] = df[thiscol].map("{}".format)
    df = clear_nan(df)
    # replace all spaces with newline
    #df.rename(columns={x:x.replace(" ","\n") for x in df.columns}, inplace=True)
    #print (df)
    return tabulate(df, tablefmt="simple", headers="keys")

def parse_includes(thisfile, verbose=False):
    '''
        read file and return text with includes
    '''
    if os.path.exists(thisfile):
        try:
            with io.open(thisfile, "r", encoding="utf-8") as f:
                text = f.read()
        except:
            print ("\n\n*** INITIAL I/O PROBLEM (COULD BE DROPBOX): {} ***\n\n".format(thisfile))
            sys.exit()
            return ""
    else:
        print ("\n\n*** INITIAL FILE NOT FOUND: {} ***\n\n".format(thisfile))
        sys.exit()
        return ""

    filedir = os.path.split(os.path.abspath(thisfile))[0]
    with cd(filedir):
        includes = get_includes(text)
        additionalfiles = list(set(includes.values()))
        if verbose:
            print ("[include.py] working in: {}".format(filedir))
        for filepath in additionalfiles:
            if os.path.exists(filepath):
                if verbose:
                    print ("[include.py] including:", filepath)
                if not filepath.endswith('.xlsx'): # don't try to read excel directly
                    newtext = parse_includes(filepath)
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
                        newtext = include_df(filepath, filetype="xlsx")
                    elif filepath.endswith('.csv'):
                        newtext = include_df(filepath, filetype="csv")
                    elif filepath.endswith('.tsv'):
                        newtext = include_df(filepath, filetype="tsv")
                    text = text.replace(inc, newtext)
    return text

def stripcomments(thistext):
    return re.sub(comments, "", thistext)

def save_new(thisfile, outputfile="auto", stripc=False, verbose=False):
    text = parse_includes(thisfile, verbose=verbose)
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
    parser.add_argument('-s', '--stripcomments', action="store_true", default=False, help='stripcomments')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='verbose')
    args = parser.parse_args()
    save_new(args.filename, outputfile="auto", stripc = args.stripcomments, verbose=args.verbose)


















