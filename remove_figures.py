import os
import re
import json

'''
remove figures from a markdown file
but leave the captions in place
'''

#-------------------
def preext(filepath, thisext):
    lastdot = filepath.rfind('.')
    outlist = [filepath[:lastdot], thisext, filepath[lastdot+1:]]
    return ".".join(outlist)
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
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='verbose')
    args = parser.parse_args()

    if args.outputfile == "auto":
        outputfile = preext(args.filename, "nofig")
    else:
        outputfile = args.outputfile

    with open(args.filename) as f:
        text = f.read()

    replacefile = os.path.join(os.path.split(args.filename)[0],"replace.json")
    if os.path.exists(replacefile):
       with open(replacefile) as f:
            replacedict = json.load(f)

    figureformat = '!\[[\s\S]+?\]\([\s\S]+?\).*?\n'
    figs = re.findall(figureformat, text)

    for x in figs:
        captions = re.findall('!\[[\s\S]+?\]\(', x)
        caption = ""
        if len(captions)>0:
            caption = captions[0][2:-2]
        labels = re.findall('\)\{.*?\}', x)
        outlabel = ""
        if len(labels)>0:
            label = labels[0][2:-1]
            labeltype = label[1:].split(":")[0]
            if labeltype in labeldic.keys():
                outlabel += labeldic[labeltype]+" "
            labelref = "@{}".format(label[1:])
            if labelref in replacedict.keys():
                outlabel += replacedict[labelref]
        text = text.replace(x, "{}. {}\n\n".format(outlabel,caption.strip()))

    with open(outputfile,"w") as o:
        o.write(text)





