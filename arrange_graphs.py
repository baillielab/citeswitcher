#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

import os
import re
import sys
import math
import subprocess
#-----------------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory',    default="/Users/jkb/Dropbox/2_proposals_and_publications/2_SUBMITTED_PAPERS/genomicc-manuscript3/manuscript/supplementary-figures/forestplots-hospitalised")
parser.add_argument('-c', '--columns',  type=int,  default=2)
parser.add_argument('-r', '--rows',  type=int,  default=4, help="max rows per div")
parser.add_argument('-pad', '--padding',  type=int,  default=0, help="percentage to add as padding")
parser.add_argument('-cap', '--caption',   default="")
parser.add_argument('-e', '--extension',  default="png")
parser.add_argument('-o', '--outputfile',  default="autoarranged.md", help="must end with .md")
parser.add_argument('-p', '--pandoc_outputs',    action='append', default=["pdf"],       help='append as many pandoc formats as you want: pdf docx html txt md tex')
args = parser.parse_args()
#-----------------------------
def natural_key(string_):
    """See https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
#-----------------------------

outtext = ""
figlabel = os.path.split(args.directory)[1]
files = [x for x in os.listdir(args.directory) if x.endswith(args.extension)]
files.sort(reverse=False, key=natural_key)
for i,f in enumerate(files):
    print (i,f)
totrows = math.ceil(float(len(files))/args.columns)
totdivs = math.ceil(float(totrows)/args.rows)
width = math.floor(float(100-args.padding)/args.columns)
i=0
for d in range(totdivs):
    outtext += '<div id="fig:%s-%s-auto">\n\n'%(figlabel,d)
    for r in range(args.rows):
        for c in range(args.columns):
            if i>=len(files):
                break
            print (i,d,r,c)
            link = "![](%s){#fig:%s-%s-%s-auto width=%s%s}"%(os.path.relpath(os.path.join(args.directory, files[i])), figlabel, d, i, width, '%')
            if c>0:
                link = "\\ "+link
            outtext += link
            i+=1
        outtext += "\n\n"
    outtext += '\n{}\n\n</div>\n\n'.format(args.caption)
with open(args.outputfile,"w") as o:
    o.write(outtext)

if "tex" in args.pandoc_outputs:
    cmd = "pandoc {} -o {}".format(args.outputfile, args.outputfile.replace(".md",".tex"))
    print (cmd)
    subprocess.call(cmd, shell=True)

if "pdf" in args.pandoc_outputs:
    cmd = "pandoc {} -o {}".format(args.outputfile, args.outputfile.replace(".md",".pdf"))
    print (cmd)
    subprocess.call(cmd, shell=True)



















