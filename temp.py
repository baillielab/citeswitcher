import os
import re

text = '''
here's some text

and then an image
![caption](img/yarr.eps){#fig:ooo}
![caption](img/yarr.svg){#fig:ooo}
![](img/yarr.svg){#fig:ooo width=45%}\ ![](img/yarr.svg){#fig:ooo width=45%}
![](img/yarr.svg){#fig:ooo width=45%}\ ![](img/yarr.pes){#fig:ooo width=45%}

'''


def replace_svgs(thistext, thispath):
    '''inkscape -D --export-filename=recruitmentbyweek.pdf recruitmentbyweek.svg '''
    imagelink = '!\[[\S]*?\]\([\S]*?.svg\)'
    svgcalls = re.findall(imagelink, thistext)
    for call in svgcalls:
        filename = re.findall('\([\S]*?.svg\)',call)[0][1:-1]
        pdf_filename = filename.replace(".svg",".pdf")
        callpdf = call.replace(filename, pdf_filename)
        cmd = "inkscape -D --export-filename={} {}".format(os.path.join(thispath,pdf_filename), os.path.join(thispath,filename))
        thistext = thistext.replace(call, callpdf)
    return thistext

print (replace_svgs(text, "/path/here/"))