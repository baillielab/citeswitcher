import os
import re
import json
import slugify
import subprocess
from get_figure_numbers import get_id, mdregex
import include

'''
remove figures from a markdown file
put figures in separate files named by figure number from replace.json
but leave the captions in place
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

def fix_permissions(this_path):
    os.system("/bin/chmod 755 %s"%(this_path))

def check_dir(this_dir):
    if not os.path.isdir(this_dir):
        os.mkdir(this_dir)
    fix_permissions(this_dir)

def make_pdf(outlabel, this_markdown):
    filepath = os.path.join(args.outputdir, slugify.slugify(outlabel)+".md")
    outfile = filepath[:-3]+".pdf"
    # write markdownfile
    with open(filepath,'w') as o:
        o.write(
'''---
figureTitle: Item
geometry: top=2.3cm, bottom=2.3cm, left=2.3cm, right=2.3cm
header-includes:
  - \\sloppy
  - \\usepackage{cleveref}
---'''
)
        o.write('\n\pagenumbering{gobble}\n\n')
        o.write(this_markdown)
    # call pandoc
    cmd = "pandoc {} -o {} --filter pandoc-crossref".format(filepath, outfile)
    print (cmd)
    subprocess.call(cmd, shell=True)
    if args.croplegend:
        svgfile = filepath[:-3]+".svg"
        cmd = "inkscape {} -D --export-type=svg --export-overwrite --pdf-poppler --export-margin=10".format(outfile)
        print (cmd)
        subprocess.call(cmd, shell=True)
    if args.inkscapeshrink:
        cmd = "inkscape {} -D --export-type=pdf --export-overwrite --pdf-poppler --export-margin=10".format(outfile)
        print (cmd)
        subprocess.call(cmd, shell=True)
    if args.eps:
        cmd = "inkscape {} -D --export-type=eps --export-overwrite --pdf-poppler --export-margin=10".format(outfile)
        print (cmd)
        subprocess.call(cmd, shell=True)
    if args.tiff:
        tiffile = filepath[:-3]+".svg"
        cmd = "gm convert -density 1400x1400  {} {}.tif".format(outfile, tiffile)
        print (cmd)
        subprocess.call(cmd, shell=True)
    if not (args.messy):
        os.remove(filepath)


#-------------------

labeldic = {
    "fig":"Figure",
    "tbl":"Table",
    "eq":"Equation",
    "edt":"Extended Data Table",
    "edf":"Extended Data Figure",
}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-o', '--outputfile', default="auto",   help='outputfile')
    parser.add_argument('-d', '--outputdir', default="auto-generated",   help='outputfile')
    parser.add_argument('-m', '--messy', action="store_true", default=False, help='keep intermediate markdown files')
    parser.add_argument('-r', '--removecaptions', action="store_true", default=False, help='remove captions - not great usually')
    parser.add_argument('-c', '--croplegend', action="store_true", default=False, help='use inkscape to cut to image only')
    parser.add_argument('-i', '--inkscapeshrink', action="store_true", default=False, help='use inkscape to cut to image only')
    parser.add_argument('-e', '--eps', action="store_true", default=False, help='use inkscape to cut and generate eps')
    parser.add_argument('-t', '--tiff', action="store_true", default=False, help='use inkscape to cut and generate tiff')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='verbose')
    parser.add_argument('-q', '--quick', action="store_true", default=False, help='quick run with no image conversion')
    parser.add_argument('-cap', '--captions_only', action="store_true", default=False, help='report captions only')
    args = parser.parse_args()

    if args.outputfile == "auto":
        outputfile = preext(args.filename, "nofig")
    else:
        outputfile = args.outputfile

    check_dir(args.outputdir)

    with cd(os.path.split(os.path.abspath(args.filename))[0]):
        '''
        with open(args.filename) as f:
            text = f.read()
        '''

        text = include.parse_includes(args.filename)

        captionlist = []

        replacedict = {}
        replacefile = os.path.join(os.path.split(args.filename)[0],"replace.json")
        if os.path.exists(replacefile):
           with open(replacefile) as f:
                replacedict = json.load(f)

        figureformat = '<div[\s\S]+?<\/div>'
        figs = re.findall(figureformat, text)
        for divfig in figs:
            figlines = divfig.split("\n")
            outlabel = get_id(figlines[0], "div")
            if outlabel in replacedict:
                outlabel = replacedict[outlabel]
            captionlines = [x for x in figlines[1:-1] if not x.startswith("![]")]
            captionlines = [x for x in captionlines if len(x.strip())>0]
            caption = "".join(captionlines)
            newtext = "{}. {}\n\n".format(outlabel,caption.strip())
            text = text.replace(divfig, newtext)
            captionlist.append(newtext)
            divfig = []
            cap = True
            for x in figlines:
                if x not in captionlines:
                    divfig.append(x)
                elif cap:
                    divfig.append(outlabel)
                    cap = False
            divfig = "\n".join(divfig)
            if not args.quick:
                make_pdf(outlabel, divfig)

        figureformat = '!\[[\s\S]+?\]\([\s\S]+?\).*?\n'
        figs = re.findall(figureformat, text)
        for mdfig in figs:
            captions = re.findall('!\[[\s\S]+?\]\(', mdfig)
            caption = ""
            if len(captions)>0:
                caption = captions[0][2:-2]
            outlabel = get_id(re.findall(mdregex, mdfig)[0], "md")
            if outlabel in replacedict:
                outlabel = replacedict[outlabel]
            newtext = "{}. {}\n\n".format(outlabel,caption.strip())
            text = text.replace(mdfig, newtext)
            captionlist.append(newtext)
            if not args.quick:
                make_pdf(outlabel, mdfig.replace(caption,""))

        '''
        # unnecessary because we now parse includes
        tableformat = '\.csv\!\}[ ]*\n: [\s\S]*?\{#.*?:.*?}'
        tabs = re.findall(tableformat, text)
        for tab in tabs:
            outlabel = get_id(re.findall(mdregex, tab)[0], "md")
            if outlabel in replacedict:
                outlabel = replacedict[outlabel]
            mdlabelgex = r'\{#.*?:.*?}'
            newtab = re.sub(mdlabelgex, "", tab, count=1)
            text = text.replace(tab, ".csv!}\n\n%s. %s\n"%(outlabel,newtab[9:]))
        '''

        tableformat = '\n(?:\|.*?\|\n?)+\n: [\s\S]*?\{#.*?:.*?}'
        tabs = re.findall(tableformat, text)
        for tab in tabs:
            captions = re.findall('\n: [\s\S]*?\{', tab)
            caption = ""
            if len(captions)>0:
                caption = captions[0][2:-2]
            outlabel = get_id(re.findall(mdregex, tab)[0], "md")
            if outlabel in replacedict:
                outlabel = replacedict[outlabel]
            mdlabelgex = r'\{#.*?:.*?}'
            newtab = re.sub(mdlabelgex, "", tab, count=1)
            newtext = "\n%s. %s\n"%(outlabel,caption.strip())
            text = text.replace(tab, newtext)
            captionlist.append(newtext)

        with open(outputfile,"w") as o:
            if args.captions_only:
                o.write("\n\n".join(captionlist))
            else:
                o.write(text)















