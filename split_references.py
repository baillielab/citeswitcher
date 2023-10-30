import os
from get_sections import get_section

'''
For Nature papers. 
Take reference list from end of a pandoc-generated markdown file
Split it at the first 50 references. 
Insert the first n references wherever this tag appears in the input file:
<!--INSERTREFS:<n>-->
Only one tag is allowed
'''

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

def preext(filepath, thisext):
    lastdot = filepath.rfind('.')
    outlist = [filepath[:lastdot], thisext, filepath[lastdot+1:]]
    return ".".join(outlist)
#-----------------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-o', '--outputfile', default="auto",   help='outputfile')
    parser.add_argument('-m', '--messy', action="store_true", default=False, help='keep intermediate markdown files')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='verbose')
    args = parser.parse_args()

    if args.outputfile == "auto":
        args.outputfile = preext(args.filename, "splitrefs")

    with cd(os.path.split(os.path.abspath(args.filename))[0]):
        with open(args.filename) as f:
            lines = f.readlines()

        numrefs = 0
        insertline = 0
        tagstart = "<!--INSERTREFS:"
        for i,line in enumerate(lines):
            if line.startswith(tagstart):
                numrefs = int(line.replace(tagstart,"").replace("-->",""))
                insertline = i
                print (numrefs, insertline)
                break # no reason why we'd split twice. 

        reflines = get_section(lines, "References")
        reftext = "".join(lines[x] for x in reflines)
        refs = reftext.split("\n\n")

        lines[insertline] = "\n\n".join(refs[:numrefs])
        del lines[reflines[0]:reflines[-1]]
        lines[reflines[0]] = "\n\n".join(refs[numrefs:])

        with open(args.outputfile, "w") as o:
            o.write("\n".join(lines))














