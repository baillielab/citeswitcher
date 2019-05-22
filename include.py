#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

import io
import os
import re
#-----------------------------
'''
RULES:
INCLUDESECTION IN LINE
#INCLUDESECTION NOT IN LINE
WHOLE LINE IS REPLACED
FILENAME IS NEXT WORD AFTER INCLUDESECTION
'''
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
#-----------------------------
comments = r'<!--[\s\S]+?-->'
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
    return x[x.index("INCLUDESECTION")+1].strip().replace("'",'').replace('"','')

def get_includes(thistext):
    includelist = [x for x in re.findall(comments, thistext) if "INCLUDESECTION" in x and not "#INCLUDESECTION" in x]
    return includelist

def parse_includes(thisfile, verbose=False):
    '''
        read file and return text with includes
    '''
    with io.open(thisfile, "r", encoding="utf-8") as f:
        text = f.read()

    filedir = os.path.split(os.path.abspath(thisfile))[0]
    with cd(filedir):
        includes = get_includes(text)
        additionalfiles = list(set([get_filename(x) for x in includes if x != None]))
        if verbose:
            print ("working in: {}".format(filedir))
            print ("Adding", additionalfiles)
        for filepath in additionalfiles:
            filepath = os.path.normpath(filepath)
            if os.path.exists(filepath):
                newtext = parse_includes(filepath)
            else:
                print ("\n\n*** INCLUDE FILE NOT FOUND: {} ***\n\n".format(filepath))
                continue
            for inc in includes:
                if get_filename(inc) == filepath:
                    text = text.replace(inc, newtext)
    return text

def stripcomments(thistext):
    return re.sub(comments, "", thistext)

def save_new(thisfile, outputfile="auto", stripc=False):
    text = parse_includes(thisfile)
    if outputfile == 'auto':
        outputfile = preext(thisfile, 'inc')
    if stripc:
        text = stripcomments(text)
    with open(outputfile,'w') as o:
        o.write(text)
    return outputfile

#-----------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-s', '--stripcomments', action="store_true", default=False, help='stripcomments')
    args = parser.parse_args()
    save_new(args.filename, outputfile="auto", stripc = args.stripcomments)


















