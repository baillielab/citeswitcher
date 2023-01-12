import os
import re
import json
import collections
#--------------
import include

'''
read a tex or md file
write a json file for replacing contents
'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', default=None,   help='filename')
parser.add_argument('-o', '--outputfile', default="replace.json",   help='outputfile')
args = parser.parse_args()

htmldivstart = r'<div.*?id.*?"fig:.*?>'
htmldiv = htmldivstart+r'[\S\s]*?<\/div>'

with open(args.filename) as f:
    text = f.read()
text = include.parse_includes(args.filename)

stemlabels = {
    "fig": "Supplementary Figure",
    "suppfig": "Supplementary Figure",
    "tbl": "Supplementary Table"
}

# remove all html comments
text = re.sub(r"<!--[\s\S]+?-->", "", text)
# remove all latex comments
text = re.sub(r"(?<!\\)%.*", "", text)

html_labels = re.findall(htmldivstart, text)
stems = list(set([x[4:-1].replace("'",'"').split("id")[1].split('"')[1].strip().split(":")[0] for x in html_labels]))
workingtext = re.sub(htmldiv, "", text) # remove all htmldivs from text (so that subfig refs aren't counted)
tex_labels = re.findall(r'\\label\{.*?\}', workingtext)
stems = list(set([x[7:].split(":")[0] for x in tex_labels]+stems))
md_labels = re.findall(r'#.*?:.*?[ }]', workingtext)
stems = list(set([x[1:-1].split(":")[0] for x in md_labels]+stems))

labels = html_labels + tex_labels + md_labels
lines = text.split("\n")
sortdict = {}
for label in labels:
    sortdict[label] = [lines.index(x) for x in lines if label in x]

results = collections.OrderedDict()
for s in stems:
    if s not in stemlabels:
        stemlabels[s] = s # use stem if no label entered above
    # assign numbers to entities
    i=1
    for item in [(k, sortdict[k]) for k in sorted(sortdict, key=sortdict.get, reverse=False)]:
        label = item[0]
        print (label)
        if "{%s:"%(s) in label:
            results["@"+label[7:-1]] = "{} {}".format(stemlabels[s], i)
            i+=1
        elif label.startswith("#") and s in label:
            results["@"+label[1:-1]] = "{} {}".format(stemlabels[s], i)
            i+=1
        elif label.startswith("<div") and s in label:
            results ["@"+label[4:-1].replace("'",'"').split("id")[1].split('"')[1].strip()] = "{} {}".format(stemlabels[s], i)
            i+=1

for item in results:
    print (item, results[item])

alljson = {}
if os.path.exists(args.outputfile):
    with open(args.outputfile) as f:
        alljson = json.load(f)
for item in results:
    alljson[item]=results[item]
with open(args.outputfile,"w") as o:
    json.dump(alljson,o, indent=4)






