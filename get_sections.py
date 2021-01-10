#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
'''
save a new file containing:
header
+ each section, up until the next section with at the same level or less  (# ## ### #### #####)
(NB comments stripped)

if provided with a list of people or json file:
search for a list of people (lowercase search) and make a new file containing each line with that person's name in it
together with the relevant section hierarchy
'''

import io
import os
import re
import json
import collections
#-----------------------------
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
                section_open = True
                for j in range(i+1, len(lines)):
                    line = lines[j]
                    if line.startswith('#'):
                        hashcount = line.split(' ')[0].count('#')
                        if hashcount <= indexhashcount:
                            sectionlines += range(i,j)
                            section_open = False
                            break
                if section_open:
                    sectionlines += range(i,j)
                    section_open = False
                i=j # skip previous section

    return sectionlines

def get_sections_for_names(theselines, thesenames):
    '''
        theselines = list
        thesenames = list of lists or simple list
    '''
    namedic = {}
    for n in thesenames:
        if type(n) ==list:
            namedic[n[0]] = n
        else:
            namedic[n]=[n]
    outdic = collections.OrderedDict()
    for n in namedic:
        outdic[namedic[n][0]] = collections.OrderedDict()
    maxsections = 9
    running_section_list = ['' for i in range(maxsections)]
    for line in theselines:
        if line.startswith("#"):
            hashcount = line.split(' ')[0].count('#')
            running_section_list = running_section_list[:hashcount-1] + [line[hashcount+1:]] + ['' for i in range(maxsections-hashcount)]
        for n in namedic:
            for synonym in namedic[n]:
                if synonym.lower() in line.lower():
                    sec = ': '.join([x for x in running_section_list if x!=''])
                    try:
                        outdic[n][sec]
                    except:
                        outdic[n][sec]=[]
                    if line not in outdic[n][sec]:
                        outdic[n][sec].append(line)
                    break # only need to count a synonym once
    return outdic

#-----------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-o', '--outputfile', default='auto',   help='filename')
    parser.add_argument('-g', '--getfilename', action="store_true", default=False,    help='returns filename only')
    parser.add_argument('-s', '--sections', action='append', default=[], help='use this to append as many values as you want')
    parser.add_argument('-po', '--peopleoutputfile', default='auto',   help='filename')
    parser.add_argument('-sn', '--searchnames', action='append', default=[], help='use this to append as many values as you want')
    parser.add_argument('-sj', '--searchnames_jsonfile', default='auto', help="json file containing synonyms dict")
    args = parser.parse_args()

    args.filename = os.path.abspath(args.filename)
    if args.outputfile == 'auto':
        fn = os.path.split(args.filename)[1]
        ofn = 'auto_'+fn.split('.')[0]
        for s in args.sections:
            ofn += '_'+s.replace(' ','')
        ofn += '.txt'
        args.outputfile = os.path.join(os.path.split(os.path.abspath(args.filename))[0], ofn)
    else:
        args.outputfile = os.path.abspath(args.outputfile)


    if args.peopleoutputfile == 'auto':
        args.peopleoutputfile = os.path.join(os.path.split(os.path.abspath(args.searchnames_jsonfile))[0],'people.md')
    else:
        args.peopleoutputfile = os.path.abspath(args.peopleoutputfile)


    if len(args.searchnames)>0:
        sn = args.searchnames
    elif os.path.exists(args.searchnames_jsonfile):
        with open(args.searchnames_jsonfile) as f:
            j = json.load(f)
        sn = j['names']
    else:
        sn=[]

    if args.getfilename:
        print (args.outputfile)
    elif len(args.sections)>0 or len(sn)>0:
        print ("Sections:", args.sections)
        print ("Searchnames:", sn)
        linelist = []
        with cd(os.path.split(os.path.abspath(args.filename))[0]):
            text = include.parse_includes(args.filename)
            text = include.stripcomments(text)
            lines= text.split('\n')
            #linelist += get_header(lines)
            for this_section in args.sections:
                linelist += get_section(lines, this_section)
            linelist = list(set(linelist))

            if len(linelist)>0:
                with open(args.outputfile,'w') as o:
                    for i in sorted(linelist):
                        o.write(lines[i]+'\n')

            if len(sn)>0:
                people = get_sections_for_names(lines, sn)

                sqb_format = '\(.*?\)'
                sqb = re.findall(sqb_format, text)
                sqb = [x.replace("\)",")")[1:-1].split(",") for x in sqb]
                sqb = [item.strip() for sublist in sqb for item in sublist]
                allnames = [item.strip() for sublist in sn for item in sublist]
                new_sqb = list(set(sqb) - set(allnames))
                new_sqb = ",\n".join(["\t\t[\"{}\"]".format(x) for x in new_sqb])
                if len(new_sqb)>0:
                    print ("\n\n Possible new names to add to name list?\n\n{}\n".format(new_sqb))

                with open(args.peopleoutputfile,'w') as o:
                    for p in people:
                        if len(people[p])>0:
                            o.write("# {}\n\n".format(p))
                            for sec in people[p]:
                                o.write("## {}\n\n".format(sec))
                                for line in people[p][sec]:
                                    line = line.strip()
                                    if line.startswith("-"):
                                        line = line[1:]
                                        line = line.strip()
                                    o.write("- {}\n".format(line.strip()))
                                o.write("\n")



















