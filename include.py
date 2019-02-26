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

def get_include(thisline):
    if "INCLUDESECTION" in thisline and not "#INCLUDESECTION" in thisline:
        x = thisline.split(' ')
        return x[x.index("INCLUDESECTION")+1].strip().replace("'",'').replace('"','')

def parse_includes(thisfile, verbose=False):
    '''
        read file and return lines with includes
    '''
    with io.open(thisfile, "r", encoding="utf-8") as f:
        lines = f.readlines()
    filedir = os.path.split(os.path.abspath(thisfile))[0]
    with cd(filedir):
        additionalfiles = [get_include(x) for x in lines]
        additionalfiles = list(set([x for x in additionalfiles if x != None]))
        if verbose:
            print ("working in: {}".format(filedir))
            print ("Adding", additionalfiles)
        for filepath in additionalfiles:
            filepath = os.path.normpath(filepath)
            if os.path.exists(filepath):
                newlines = parse_includes(filepath)
                text = ''.join(newlines)
            else:
                print ("\n\n*** INCLUDE FILE NOT FOUND: {} ***\n\n".format(filepath))
                continue
            for i,line in enumerate(lines):
                if get_include(line) == filepath:
                    lines[i] = text
        return lines

def stripcomments(thistext):
    comments = r'<!--[\s\S]+?-->'
    return re.sub(comments, "", thistext)

def save_new(thisfile, outputfile="auto", stripc=False):
    lines = parse_includes(thisfile)
    if outputfile == 'auto':
        outputfile = preext(thisfile, 'inc')
    text = ''.join(lines)
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
    save_new(args.filename, stripc = args.stripcomments)


















