import subprocess
import os

import citefunctions
config = citefunctions.getconfig(os.path.join(os.path.dirname(__file__), 'config.json'))
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
citefunctions.callpandoc(filename, '.tex')
citefunctions.callpandoc(filename, '.txt')





