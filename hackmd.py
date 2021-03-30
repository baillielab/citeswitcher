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
if args.direct:
    args.direct = args.direct.replace("/edit", "")
    dfile = 'hackmd_download.md'
    cmd = 'wget {}/download -O {}'.format(args.direct, dfile)
    print (cmd)
    subprocess.call(cmd, shell=True)
    citefunctions.make_output(dfile, pathtopandoc=args.pathtopandoc, localbibonly=args.localbibonly, outputformats=outputformats)
else:
    hackmd = os.path.expanduser(config['hackmddir'])
    with citefunctions.cd(hackmd):
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

                citefunctions.make_output(os.path.join(thispath, mdfile), pathtopandoc=args.pathtopandoc, localbibonly=args.localbibonly, outputformats=outputformats)
                # TODO: push latex compilation output to the same folder so that users can see it























