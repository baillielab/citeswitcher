import os
import io
import subprocess
import citefunctions
config = citefunctions.getconfig(os.path.join(os.path.dirname(__file__), 'config.json'))
#-------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file1',    help='filepath')
parser.add_argument('-g', '--file2',    help='filepath')
parser.add_argument('-m', '--messy', action="store_true", default=False, help='disable clean up of intermediate files')
args = parser.parse_args()
#-------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
workingdir, filename = os.path.split(os.path.abspath(args.file2))
#-------------
diff_file = "diff.tex"
auxdir = "diffaux"
#-------------
f1 = os.path.abspath(args.file1)
f2 = os.path.abspath(args.file2)
t1 = citefunctions.newext(f1, ".citemd.tex")
t2 = citefunctions.newext(f2, ".citemd.tex")
with citefunctions.cd(scriptpath):
    cmd = "python fixcitations.py -p tex -f {}".format(f1)
    subprocess.call(cmd, shell=True)
    cmd = "python fixcitations.py -p tex -f {}".format(f2)
    subprocess.call(cmd, shell=True)
    cmd = "latexdiff --flatten --encoding=utf8 {} {} > {}".format(t1, t2, os.path.join(workingdir, "diff.tex"))
    print (cmd)
    subprocess.call(cmd, shell=True)

with citefunctions.cd(workingdir):
    cmd = "pdflatex -interaction=batchmode {}".format(diff_file)
    print (cmd)
    subprocess.call(cmd, shell=True)

    if not args.messy:
        for filepath in [t1,t2]:
            cmd = "rm {}".format(filepath)
            subprocess.call(cmd, shell=True)
