import os
import re
import json
import collections

'''
read a tex or md file
write a json file for replacing contents
'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', default=None,   help='filename')
parser.add_argument('-o', '--outputfile', default="replace.json",   help='outputfile')
args = parser.parse_args()

with open(args.filename) as f:
    text = f.read()

stemlabels = {
    "suppfig": "Supplementary Figure",
    "tbl": "Supplementary Table"
}

# remove all latex and html comments
text = re.sub(r"(?<!\\)%.*", "", text)
text = re.sub(r"<!--[\s\S]+?-->", "", text)

labels = re.findall(r'\\label\{.*?\}', text)
stems = list(set([x[7:].split(":")[0] for x in labels]))

labels += re.findall(r'#.*?:.*?[ }]', text)
stems += list(set([x[1:-1].split(":")[0] for x in labels]))

results = collections.OrderedDict()
for s in stems:
    if s not in stemlabels:
        stemlabels[s] = s # use stem if no label entered above
    i=1
    for label in labels:
        if "{%s:"%(s) in label:
            results["@"+label[7:-1]] = "{} {}".format(stemlabels[s], i)
            i+=1
        elif label.startswith("#") and s in label:
            results["@"+label[1:-1]] = "{} {}".format(stemlabels[s], i)
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






