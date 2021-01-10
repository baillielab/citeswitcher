#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

'''
    download a google doc from url
'''

import os
import sys
import platform
import subprocess
from datetime import datetime
import citefunctions
#-----------------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--direct',    help='enter direct hackmd link')
parser.add_argument('-l', '--localbibonly', action="store_true", default=False, help='use only local bib file')
parser.add_argument('-img', '--convertimages', action="store_true", default=False, help='automatically convert images')
parser.add_argument('-ptp', '--pathtopandoc', default='pandoc', help='specify a particular path to pandoc if desired')
args = parser.parse_args()
#-----------------------------
outputformats = [
        ".pdf",
        ".docx",
        ".html",
    ]
#-----------------------------
archive_dir = "archive"
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
config = citefunctions.getconfig()
#-----------------------------
class cd:
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)
    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)
    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

#-----------------------------

def make_output(thispath, thisfile, pathtopandoc=args.pathtopandoc):
    # run fix citations in messy mode
    extra_args = "-x -svg -pm -flc "
    # -x indicates xelatex mode. Handles special characters. Crashes sid.
    # -pm indicates that pandoc mermaid is used
    if args.localbibonly:
        extra_args += (" -l")
    cmd = '{} ~/Dropbox/3_scripts_and_programs/citeswitcher/fixcitations.py {} -f {} -m {} -ptp {} '.format(
        sys.executable,
        extra_args,
        os.path.join(thispath, thisfile),
        " ".join(["-p "+x.replace(".","") for x in outputformats]),
        pathtopandoc
        )
    print (cmd)
    subprocess.call(cmd, shell=True)

#-----------------------------
if args.direct:
    args.direct = args.direct.replace("/edit", "")
    dfile = 'hackmd_download.md'
    cmd = 'wget {}/download -O {}'.format(args.direct, dfile)
    print (cmd)
    subprocess.call(cmd, shell=True)
    make_output(dfile)
else:
    hackmd = os.path.expanduser(config['hackmddir'])
    with cd(hackmd):
        dirnames = [x for x in os.listdir() if os.path.isdir(x)]
        print ("running over:", dirnames)
        for d in dirnames:
            subdirs = [x for x in os.listdir(d) if os.path.isdir(os.path.join(d,x))]
            print ("chosendirs:", args.chosendirs)
            if len(args.chosendirs) > 0:
                subdirs = [x for x in subdirs if x in args.chosendirs]
            print ("subdirs chosen:", subdirs)
            for s in subdirs:
                thispath = os.path.join(d,s)
                archive_dir = os.path.abspath(os.path.join(thispath, 'archive'))
                if not os.path.exists(archive_dir):
                    os.mkdir(archive_dir)
                mdfile = s+".md"
                datestring = datetime.now().strftime("%Y%m%d")
                # archive today's file
                if os.path.exists(os.path.join(thispath,mdfile)):
                    cmd = 'cp {} {}'.format(os.path.join(thispath,mdfile), os.path.join(archive_dir, s+datestring+'.md'))
                    print (cmd)
                    subprocess.call(cmd, shell=True)

                # download md file
                hackmd_link = "https://hackmd.io/@{}/{}".format(d,s)
                cmd = 'wget {}/download -O {}'.format(hackmd_link, os.path.join(thispath,mdfile))
                print (cmd)
                subprocess.call(cmd, shell=True)

                make_output(thispath, mdfile)























