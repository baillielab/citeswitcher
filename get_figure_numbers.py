import os
import re
import json
import collections
import oyaml as yaml
#--------------
import include
from citefunctions import getyaml, mergeyaml

'''
read a tex or md file
split it into sections if <!--NEWYAML\n<any-number-of-yaml-lines-with-no-comments>\n--> appears
write a json file for replacing contents
'''

htmldivstart = r'<div.*?id.*?"fig:.*?>'
htmldiv = htmldivstart+r'[\S\s]*?<\/div>'
mdregex = r'#.*?:.*?[ }]'
latexregex = r'\\label\{.*?\}'

def get_divtype(thistext):
    return thistext[4:-1].replace("'",'"').split("id")[1].split('"')[1].strip().split(":")[0]

def get_mdtype(thistext):
    return thistext[1:-1].split(":")[0]

def get_textype(thistext):
    return thistext[7:].split(":")[0]

def get_id(thistext, idtype="md"):
    if idtype == "md":
        return "@"+thistext[1:-1]
    if idtype == "div":
        return "@"+thistext[4:-1].replace("'",'"').split("id")[1].split('"')[1].strip()
    if idtype == "tex":
        return "@"+thistext[7:-1]
    else:
        return "get_id() failed because idtype not recognised"


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-k', '--keep_original',   action="store_true", default=False,  help='keep contents of original file')
    parser.add_argument('-o', '--outputfile', default="replace.json",   help='outputfile')
    parser.add_argument('-v', '--verbose',    action="store_true", default=False,    help='increases verbosity')
    args = parser.parse_args()

    with open(args.filename) as f:
        text = f.read()
    text = include.parse_includes(args.filename)
    documentsections = text.split("<!--NEWYAML")

    results = collections.OrderedDict()
    for dsi, text in enumerate(documentsections): # dsi = documentsectionindex
        sortdict = {}
        stemlabels = {
            "fig": ["Figure", "figPrefix", 0], # label, prefix, counterstart
            "tbl": ["Table", "tblPrefix", 0],
            "eqn": ["Equation", "eqnPrefix", 0],
        }
        if dsi == 0: # this is the first section and has regular YAML
            thisyaml = getyaml(text)
        else:
            sectionyaml = getyaml("---\n{}".format(text.replace("-->","---",1)))
            thisyaml = mergeyaml(sectionyaml, thisyaml)
        for lab in stemlabels:
            prefix = stemlabels[lab][1]
            if prefix in thisyaml:
                if type(thisyaml[prefix]) is list or type(thisyaml[prefix]) is tuple:
                    stemlabels[lab][0] = min(thisyaml[prefix], key=len) # take the shortest item in list (the singular "Figure")
                else:
                    stemlabels[lab][0] = thisyaml[prefix] # take the string
        # read setcounter
        if "header-includes" in thisyaml:
            for y in thisyaml['header-includes']:
                if y.startswith("\\setcounter"):
                    y = y[12:-1].split("}{")
                    if y[0] == "table":
                        stemlabels["tbl"][2] = int(y[1])

        # remove all html comments
        text = re.sub(r"<!--[\s\S]+?-->", "", text)

        html_labels = re.findall(htmldivstart, text)
        stems = list(set([get_divtype(x) for x in html_labels]))
        workingtext = re.sub(htmldiv, "", text) # remove all htmldivs from text (so that subfig refs aren't counted)
        tex_labels = re.findall(latexregex, workingtext)
        stems = list(set([get_textype for x in tex_labels]+stems))
        md_labels = re.findall(mdregex, workingtext)
        stems = list(set([get_mdtype(x) for x in md_labels]+stems))

        labels = html_labels + tex_labels + md_labels
        lines = text.split("\n")
        
        for label in labels:
            sortdict[label] = [lines.index(x) for x in lines if label in x]

        for s in stems:
            if s not in stemlabels:
                stemlabels[s] = [s] # use stem if no label entered above
            # assign numbers to entities
            i=1
            for item in [(k, sortdict[k]) for k in sorted(sortdict, key=sortdict.get, reverse=False)]:
                label = item[0]
                if "{%s:"%(s) in label:
                    results[get_id(label, "tex")] = "{} {}".format(stemlabels[s][0], i + stemlabels[s][2])
                    i+=1
                elif label.startswith("#") and s in label:
                    results[get_id(label, "md")] = "{} {}".format(stemlabels[s][0], i + stemlabels[s][2])
                    i+=1
                elif label.startswith("<div") and s in label:
                    results [get_id(label, "div")] = "{} {}".format(stemlabels[s][0], i + stemlabels[s][2])
                    i+=1

    for i,item in enumerate(results.keys()):
        print (item, results[item])
        for j in range(i+1, len(results.keys())):
            other_item = list(results.keys())[j]
            if item.startswith(other_item) or other_item.startswith(item):
                print ("\n******\nWARNING: DUPLICATE ID STEM: {} {}\n******\n".format(item, other_item))

    alljson = {}
    if os.path.exists(args.outputfile):
        with open(args.outputfile) as f:
            alljson = json.load(f)
    for item in results:
        if args.keep_original:
            if item not in alljson:
                alljson[item]=results[item]
        else:
            alljson[item]=results[item]
    with open(args.outputfile,"w") as o:
        json.dump(alljson,o, indent=4)

    for i,item in enumerate(alljson.keys()):
        for j in range(i+1, len(alljson.keys())):
            other_item = list(alljson.keys())[j]
            if item.startswith(other_item) or other_item.startswith(item):
                print ("\n******\nWARNING: DUPLICATE ID STEM in replace.json: {} {}\n******\n".format(item, other_item))




