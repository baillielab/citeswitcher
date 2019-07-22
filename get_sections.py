#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
'''
save a new file containing:
header
+ each section, up until the next section with at the same level or less  (# ## ### #### #####)
(NB comments stripped)
'''

import io
import os
import include
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

def get_header(theselines):
    '''
        return **line numbers*** for markdown header
    '''
    firstline = theselines[0].strip()
    if firstline=="---" or firstline =="...":
        for i,line in enumerate(lines[1:]):
            if line.strip() == '...':
                #return lines[:i+2]
                return range(0,i+2)
    return []

def get_section(theselines, section_name):
    '''
        return **line numbers*** for sections to be included
    '''
    sectionlines = []
    if section_name.startswith('#'):
        section_name = ' '.join(section_name.split(' ')[1:])
    section_name  = section_name.strip()
    for i in range(len(lines)):
        line = lines[i]
        if line.startswith('#'):
            indexhashcount = line.split(' ')[0].count('#')
            hline = ' '.join(line.split(' ')[1:]).strip()
            if hline == section_name:  # section has been found. Find next section
                for j in range(i+1, len(lines)):
                    line = lines[j]
                    if line.startswith('#'):
                        hashcount = line.split(' ')[0].count('#')
                        if hashcount <= indexhashcount:
                            sectionlines += range(i,j)
                            break
                i=j # skip last section
    return sectionlines

#-----------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-o', '--outputfile', default='auto',   help='filename')
    parser.add_argument('-g', '--getfilename', action="store_true", default=False,    help='returns filename only')
    parser.add_argument('-s', '--sections', action='append', default=[], help='use this to append as many values as you want')
    args = parser.parse_args()
    if args.outputfile == 'auto':
        args.outputfile = 'auto_'+args.filename.split('.')[0]
        for s in args.sections:
            args.outputfile += '_'+s.replace(' ','')
        args.outputfile += '.txt'

    if args.getfilename:
        print (args.outputfile)
    else:
        linelist = []
        with cd(os.path.split(os.path.abspath(args.filename))[0]):
            text = include.parse_includes(args.filename)
            text = include.stripcomments(text)
            lines= text.split('\n')
            linelist += get_header(lines)
            for this_section in args.sections:
                linelist += get_section(lines, this_section)

            with open(args.outputfile,'w') as o:
                for i in sorted(list(set(linelist))):
                    o.write(lines[i]+'\n')

















