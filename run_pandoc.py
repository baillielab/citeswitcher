import subprocess
import os

import citefunctions
config = citefunctions.getconfig()
#-------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filepath',    help='filepath', default=config['testfile'])
parser.add_argument('-d', '--outputsubdir',    help='outputdir - always a subdir of the working directory', default=config['pandocsubdir'])
args = parser.parse_args()
#-------------------
workingdir, filename = os.path.split(os.path.abspath(args.filepath))
outputsubdir = os.path.join(workingdir, args.outputsubdir)
#-------------
def getext(filepath):
    return os.path.split(filepath)[-1].split('.')[-1]

def newext(filepath, thisext):
    return filepath[:filepath.rfind('.')] + thisext

def fix_permissions(this_path):
    os.system("/bin/chmod 755 %s"%(this_path))
    
def check_dir(this_dir):
    if not os.path.isdir(this_dir):
        os.mkdir(this_dir)
    fix_permissions(this_dir)

def callpandoc(f, out_ext, out_dir=outputsubdir):
    cmd = 'pandoc --filter pandoc-citeproc {} -o {}'.format(f, os.path.join(out_dir, newext(f, out_ext)))
    subprocess.call(cmd, shell=True)

#-------------

os.chdir(workingdir)
check_dir(outputsubdir)

callpandoc(filename, '.docx')
callpandoc(filename, '.pdf')

