import subprocess
import os

import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filepath',    help='filepath', default=config['testfile'])
args = parser.parse_args()
#-------------------
workingdir, filename = os.path.split(os.path.abspath(args.filepath))
#-------------
os.chdir(workingdir)
citefunctions.callpandoc(filename, '.docx')
citefunctions.callpandoc(filename, '.pdf')

