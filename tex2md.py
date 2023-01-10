import os
import re
import json
import subprocess
import collections

'''
convert a tex file to markdown
depends on all tables using latex csv input \\csvautotabular \\csvautotabularcenter
manually converts all TABLES and TABLE CROSSREFERENCES
then uses pandoc for the rest
'''

#-------------------
def getcurlies(thestring, leadintext=r"\\caption"):
    '''get matching text from a latex command within nested curly brackets'''
    res = []
    for stringmatch in re.finditer(leadintext, thestring):
        a = stringmatch.end(0)
        stack = 0
        for i in range(a, len(thestring)):
            if thestring[i] == "{":
                stack +=1
            elif thestring[i] == "}":
                if stack<=1:
                    res.append(thestring[a+1:i])
                    break
                stack -= 1
    return res

def refmdify(texref):
    if not texref.startswith("\\ref"):
        return texref
    if not texref.endswith("}"):
        return texref
    return texref.replace("\\ref{","@")[:-1]

def crossref_tex2md(thistext):
    ''' replace tex \\ref{} style references with @markdown-style ones'''
    rdic = collections.OrderedDict()
    s = r"\\ref\{[\S\s]+?\}"
    refs = re.findall(s, thistext)
    rtypes=set([])
    for r in refs:
        rdic[r]=""
        rtypes.add(r.split(":")[0].replace("\\ref{",""))
    for r in rdic:
        thistext = thistext.replace(
            r,
            refmdify(r),
            )
        if args.verbose:
            print ("replacing", r, "with", refmdify(r))
    #print ("\nReference types used:\n{}\n".format("\n".join(rtypes)))
    return thistext, rdic, rtypes

def escape_markdown_chars(thistext):
    thistext = thistext.translate(str.maketrans({
            "*":  r"\*",
            }))
    return thistext

runcounter = 0
def replace_panels(thistext, panellabel="table", paneltype="table", localrdic={}):
    global runcounter
    panelformat = r'\\begin\{%s\}[\s\S]+?\\end\{%s\}'%(panellabel,panellabel)
    panels = re.findall(panelformat, thistext)
    for t in panels:
        c = getcurlies(t, leadintext=r"\\caption")
        if len(c)>0:
            caption = c[0].strip().replace("\n"," ")
            caption = escape_markdown_chars(caption)
        else:
            caption = ""
        l = getcurlies(t, leadintext=r"\\label")
        if len(l)>0:
            label = l[0]
        else:
            label = "autolabel{}".format(runcounter)
            runcounter+=1
        if paneltype=="figure":
            p = getcurlies(t, leadintext=r"\\includegraphics")
            if len(p)>0:
                image = p[0]
            else:
                image = ""
            mdpanel = "\n![%s](%s){#%s}\n"%(caption, image, label) # if figure
        if paneltype=="table":
            p = getcurlies(t, leadintext=r"\\csvautotabular")
            if len(p)>0:
                csvfile = p[0]
            else:
                csvfile = ""
            mdpanel = "{!%s!}\n\n: %s {#%s}\n\n"%(csvfile, caption, label) # if table
        placeholder = "panelplaceholder-%s-panelplaceholder"%(label)
        localrdic[placeholder] = mdpanel
        thistext = thistext.replace(t, placeholder)
        if args.verbose:
            print ("Scheduled replacement:\n{}with:\n{}\n\n".format(t,mdpanel))
    return thistext, localrdic

#-------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-o', '--outputfile', default="auto",   help='outputfile')
    parser.add_argument('-m', '--messy', action="store_true", default=False,            help='disable clean up of intermediate files')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='verbose')
    parser.add_argument('-ptp', '--pathtopandoc', default='pandoc', help='specify a particular path to pandoc if desired')
    args = parser.parse_args()

    tempfile = args.filename.replace(".tex",".temp.tex")
    if args.outputfile == "auto":
        outputfile = args.filename.replace(".tex",".md")
    else:
        outputfile = args.outputfile
    replacefile = os.path.join(os.path.split(outputfile)[0],"replace.json")

    with open(args.filename) as f:
        text = f.read()

    # first manually convert tables and edfs etc because pandoc ignores them
    # NB use \csvautotabular{Tables/table.csv} or \csvautotabularcenter{Tables/table.csv} to include tables in latex original
    text = text.replace("csvautotabularcenter","csvautotabular")
    placeholder_replacedict={}
    for identifier,thistype in [
            ("table","table"),
            ("edt","figure"), # important: extended data tables are images of tables
            ("edf","figure"),
        ]:
        text, placeholder_replacedict = replace_panels(
                                            text,
                                            panellabel=identifier,
                                            paneltype=thistype,
                                            localrdic=placeholder_replacedict
                                            )


    text = text.replace("\\bibliography{cs}", "\nplaceholder-references-placeholder\n")
    text, refdic, reftypes = crossref_tex2md(text)
    placeholder_replacedict["placeholder-references-placeholder"] = "# References\n"
    with open(tempfile,"w") as o:
        o.write(text)

    cmd = "{} {} -o {}".format(
        args.pathtopandoc,
        tempfile,
        outputfile
        )
    print (cmd)
    subprocess.call(cmd, shell=True)

    with open(outputfile) as f:
        text = f.read()
    text = text.replace("\\@","@") # need to unescape references
    for p in placeholder_replacedict:
        panel = placeholder_replacedict[p]
        panel, prefdic, preftypes = crossref_tex2md(panel)
        refdic = refdic | prefdic # python 3.9 and above only
        reftypes = reftypes | preftypes
        text = text.replace(p, panel)
    with open(outputfile,"w") as o:
        o.write(text)

    if not(args.messy):
        cmd = "rm {}".format(tempfile)
        subprocess.call(cmd, shell=True)

    # fix replace.json file
    if len(reftypes)>0:
        print ("\nReference types used:\n{}\n\n".format("\n".join(reftypes)))
    replacements = collections.OrderedDict()
    if os.path.exists(replacefile):
        try:
            with open(replacefile) as f:
                replacements = json.load(f, object_pairs_hook=collections.OrderedDict) # load into ordered dict
        except:
            pass
    core_reftypes = ["tbl","fig","eq"] # put these first
    reftypes = core_reftypes + [r for r in reftypes if r not in core_reftypes]
    for rtype in reftypes:
        for item in refdic.keys():
            mditem = refmdify(item)
            try:
                replacements[mditem]
            except:
                if mditem.startswith("@{}".format(rtype)):
                    replacements[mditem]=mditem
    with open(replacefile,"w") as o:
        json.dump(replacements, o, indent=4)














