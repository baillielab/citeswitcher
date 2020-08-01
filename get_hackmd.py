#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

'''
    automatically read through the specified directory tree (either direct or in config.json)
    download all md files corresponding to directory names
    pdf, docx and html them
    NB:
    - read settings must be open to the world
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
parser.add_argument('-c', '--chosendirs',    action='append', default=[], help='chosen directory names')
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
config = citefunctions.getconfig(os.path.join(scriptpath, 'config.json'))
#-----------------------------
path_to_python = "python"
try:
    path_to_python = config['pathstopython'][platform.system()]
except:
    pass
print ("OS={} pathtopython={}".format(platform.system(), path_to_python))
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
def run_pandoc(thispath, thisfile, fixrefs=True):
    if fixrefs:
        # run fix citations in messy mode
        cmd = '{} ~/Dropbox/3_scripts_and_programs/citeswitcher/fixcitations.py -l -f {} -m {}'.format(
            path_to_python,
            os.path.join(thispath, thisfile),
            " ".join(["-p "+x.replace(".","") for x in outputformats])
            )
        print (cmd)
        subprocess.call(cmd, shell=True)
    else:
        # just run pandoc
        with cd(thispath):
            for ext in (outputformats):
                cmd = 'pandoc {} -o {}'.format(thisfile, thisfile.replace(".md",ext))
                print (cmd)
                subprocess.call(cmd, shell=True)
#-----------------------------
if args.direct:
    args.direct = args.direct.replace("/edit", "")
    dfile = 'hackmd_download.md'
    cmd = 'wget {}/download -O {}'.format(args.direct, dfile)
    print (cmd)
    subprocess.call(cmd, shell=True)
    run_pandoc(dfile)
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

                run_pandoc(thispath, mdfile)























