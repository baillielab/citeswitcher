#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
#version 0.1

'''
Generate timeline from:
- a bib file (optional)
- a markdown file with the following features:
	# headers indicate new timelines and can contain settings in single line JSON format e.g. {color:red}
	- list items indicate timeline features and can contain PMID or DOI citations, or manual json-style entries:
		e.g. This text {img:img/path, date:mon-yy, link:<a href="https://doi.org/10.1038/nature10921">Nature (2012)</a>}
	One way or another, all list items must have date included
	Anything can be overwritten by putting a different value *after* it in a line.

	The following are acceptable json dict contents for an individual item:
	"name": "We discovered the first human gene associated with susceptibility to severe influenza",
	"labeltext": "<a href=\"https://doi.org/10.1038/nature10921\">Nature</a>",
	"img": "img/thisimage.png"
	"date": "2012-01-01"

'''

import os
import io
import re
import sys
import json
import datetime
import subprocess
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase
from PyPDF2 import PdfFileWriter,PdfFileReader
#-------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
import bib
bib.init()
import citefunctions
config = citefunctions.getconfig()
import include
sys.path.append(os.path.join(scriptpath,"Timeline-master"))
import make_timeline
#-------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filepath', default="timelines.md")
parser.add_argument('-d', '--basedir', default="auto")
parser.add_argument('-i', '--imagedir', default="img/")
parser.add_argument('-o', '--outputfile', default="auto-tl.html")
parser.add_argument('-b', '--bibfile', default="source/cs.bib")
parser.add_argument('-m', '--masterbibfile', default="/Users/jkb/Dropbox/2_proposals_and_publications/5_BIBLIOGRAPHY/lib.pmid.bib")
parser.add_argument('-s', '--defaultstart', default="2008-01-01")
parser.add_argument('-w', '--defaultwidth', default="750")
parser.add_argument('-l', '--localbibonly', action="store_true", default=False,     help='use only local bib file')
parser.add_argument('-v', '--verbose',    action="store_true", default=False,    help='increases verbosity')
args = parser.parse_args()
#-------------------
curlybrackets = '\{[\s\S]+?\}'
#-------------------
args.filepath = os.path.abspath(args.filepath)
#-------------------
def slugify(value, allow_unicode=False):
    value = str(value)
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')

def replace_timeline(thistext, use_whole=False):
	# pmid first as they are the most likely to have errors
	workingtext = thistext
	p, workingtext = citefunctions.get_mixed_citation_blocks(workingtext)
	print ("Number blocks using pmid or doi or md:{}".format(len(p)))
	replacedict = {}
	# there are slightly different procedures for each reference type:
	for b in p:
		citedhere = citefunctions.parse_citation_block(b)
		theseids, notfound = citefunctions.pmid2id(citedhere['pmids'], citedhere['notfound']) # ids added to bib
		theseids = list(set(theseids+citedhere['ids']))
		replacedict[b] = format_tl(theseids, notfound)
	for b in replacedict:
		if replacedict[b] == 'null':
			print ("{:>50} ... left alone".format(b))
			continue
		print ("{:>50} ==> {}".format(b, replacedict[b]))
		thistext = thistext.replace(b, replacedict[b])
	return thistext

def format_tl(theseids, thesemissing=[]):
	if len(theseids) == 0:
		return 'null'
	# add to the outputdatabase
	bib.cite(theseids)
	# make a blockstring
	blockstring = ''
	blockstring += '{'+"{}".format('; '.join([format_timeline(x) for x in theseids]))+'}'
	if len(thesemissing) > 0:
		blockstring += "[***{}]".format(', '.join(thesemissing))
	return blockstring

def format_timeline(thisid):
	try:
		bib.db.entries_dict[thisid]
	except:
		print (thisid, "not found in bib database")
		return ""
	# sort weird uppercase bib labels
	if "Year" in bib.db.entries_dict[thisid]:
		bib.db.entries_dict[thisid] =  {k.lower(): v for k, v in bib.db.entries_dict[thisid].items()}
	# DATE
	# first handle a bizarre string format introduced sometimes by bparser:
	if "month" in bib.db.entries_dict[thisid]:
		if isinstance(bib.db.entries_dict[thisid]['month'], bibtexparser.bibdatabase.BibDataStringExpression):
			bib.db.entries_dict[thisid]['month'] = bib.db.entries_dict[thisid]['month'].get_value()
		thisdate = '{} {}'.format(bib.db.entries_dict[thisid]['month'][:3], bib.db.entries_dict[thisid]['year'])
		date = datetime.datetime.strptime(thisdate, '%b %Y')
	else:
		date = datetime.datetime.strptime("{}".format(bib.db.entries_dict[thisid]['year']), '%Y')
	# IMAGE
	pdf = os.path.join(args.imagedir,"{}.pdf".format(thisid))
	pdf_page = os.path.join(args.imagedir,"{}.page.pdf".format(thisid))
	png = os.path.join(args.imagedir,"{}.png".format(thisid))
	if not os.path.exists(pdf):
		if "file" in bib.db.entries_dict[thisid]:
			files = bib.db.entries_dict[thisid]["file"].split(";")
			for f in files:
				if f.endswith(".pdf"):
					if os.path.exists(f):
						cmd = 'cp "{}" "{}"'.format(f, pdf)
						subprocess.call(cmd, shell=True)
	if os.path.exists(pdf):
		if not os.path.exists(png):
			inputpdf = PdfFileReader(open(pdf, "rb"))
			output = PdfFileWriter()
			output.addPage(inputpdf.getPage(0))
			with open(pdf_page, "wb") as outputStream:
				output.write(outputStream)
			cmd = 'gm convert "{}" "{}"'.format(pdf_page, png)
			print (cmd)
			subprocess.call(cmd, shell=True)
	# FORMAT NOW
	citation_contents = []
	citation_contents.append('date: {:%Y-%m-%d}'.format(date))
	if 'doi' in bib.db.entries_dict[thisid] and 'journal' in bib.db.entries_dict[thisid]:
		citation_contents.append('labeltext: <a href="https://doi.org/{}">{}</a>'.format(bib.db.entries_dict[thisid]['doi'], bib.db.entries_dict[thisid]['journal'].capitalize()))
	if os.path.exists(png):
		citation_contents.append("img: {}".format(os.path.relpath(png)))
	formatted_citation = ", ".join([x.replace(",","") for x in citation_contents]) # need to remove commas because these are the separator
	return formatted_citation

def getsettings(text):
	'''
		read a line of format:
		[# or -] some text {json:style, key:value_pairs}
	'''
	s = re.findall(curlybrackets, text)
	settingslist = []
	if len(s)>0:
		for thesesettings in s:
			text = text.replace(thesesettings,"").strip()
			settingslist += thesesettings[1:-1].split(",") # NB no commas allowed in settings items **including in citation**
	out = {"name":text}
	for entry in settingslist:
		e = entry.split(":",1) # maximum occurence =1 
		try:
			out[e[0].strip()] = e[1].strip()
		except Exception as error:
			print ("format_timeline error in getsettings() function. Malformed instruction: {}. Skipping".format(entry))
			print (entry, e)
			print (error)
			continue
	return out

def strip_hashes(thisstring):
	while thisstring.startswith("#"):
		thisstring = thisstring[1:]
	return thisstring

def apply_timeline_defaults(thesesettings):
	defaultsettings = {
		"start":args.defaultstart,
		"width":args.defaultwidth,
		"end":"{:%Y-%m-%d}".format(datetime.date.today()),
		"tick_format" : "%B-%Y",
		}
	for key in defaultsettings:
		if key not in thesesettings:
			thesesettings[key] = defaultsettings[key]
	thesesettings["width"]=int(thesesettings["width"])
	return thesesettings

def getid(thistitle):
	for illegalchar in [" ", ",", ".", ":", ";", "&", "|", "{", "}", "(", ")"]:
		thistitle = thistitle.replace(illegalchar,"")
	return thistitle.lower().strip()
#-------------------

if __name__ == "__main__":
	print (os.path.abspath(args.filepath))
	if args.basedir == "auto":
		args.basedir = os.path.join(os.path.split(args.filepath)[0],"auto-generated")
		if not os.path.exists(args.basedir):
			os.mkdir(args.basedir)
	with citefunctions.cd(args.basedir):
		if not os.path.exists(args.bibfile):
			args.bibfile = "cs.bib"
		if not os.path.exists(args.bibfile):
			print ("overwriting bib")
			with open(args.bibfile, "w") as f:
				f.write("")
		#------
		# read input files
		if args.localbibonly:
			bib.full_bibdat = citefunctions.read_bib_files([args.bibfile])
		else:
			bib.full_bibdat = citefunctions.read_bib_files([args.bibfile, args.masterbibfile])
		bib.make_alt_dicts()
		with io.open(args.filepath, "r", encoding="utf-8") as f:
			text = f.read()
		text = citefunctions.make_unicode(text)
		text = citefunctions.readheader(text)[1]
		#------
		text = include.stripcomments(text)
		text = replace_timeline(text, "md")
		# `text` is now modified with every reference (PMID or DOI) changed into a psuedo-json {key:value} string
		#------
		with open(args.bibfile, 'w') as bf:
			bibtexparser.dump(bib.db, bf)
		#------
		# now read through new text as if it were a list of instructions in markdown format `- text {key:value}`
		lines = text.split("\n")
		settings = {}
		eras = {}
		timelines = {}
		for line in lines:
			line = line.strip()
			if line.startswith("##"): # then this is an era
				line = strip_hashes(line).strip()
				this_era = getsettings(line)
				print ("ERA:", this_era)
				eras[this_tl["name"]].append(this_era)
			elif line.startswith("#"): # then this is a title
				line = strip_hashes(line).strip()
				this_tl = getsettings(line)
				this_tl = apply_timeline_defaults(this_tl)
				settings[this_tl["name"]] = this_tl
				eras[this_tl["name"]] = []
				timelines[this_tl["name"]] = [{"name":"start","date":"{}".format(this_tl["start"])}]
			if line.startswith("-"):
				line = line[1:].strip()
				this_element = getsettings(line)
				if "date" in this_element: # no date, no play
					timelines[this_tl["name"]].append(this_element)
				else:
					print ("ignoring {} as no date found".format(this_element))
		#------
		# check the contents are OK
		for t in timelines:
			# check that color is set
			try:
				settings[t]["background"]
			except:
				settings[t]["background"] = "black"
		timelines_out = {}
		for t in timelines:
			if len(timelines[t])>1:
				timelines_out[t] = timelines[t]
		#------
		# now write the output file
		htmlout = '''
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" /> 
<script src="src/d3.v2.min.js" type="application/javascript"></script>
<script src="src/timeknots-jkb.js" type="application/javascript"></script>
</head>

<body>
	'''
		for t in timelines_out:
			htmlout += '<h3 style="color:{}">{}</h3>\n<div id="{}"></div>'.format(
				settings[t]["background"],
				t,
				getid(t)
				)
		
		htmlout += '\n</body>\n\n<script type="text/javascript">\n'

		for t in timelines_out:
			htmlout += "var {} = ".format(getid(t))
			htmlout += json.dumps(timelines_out[t], indent=4)
			htmlout += ";\n\n"
			htmlout += 'TimeKnots.draw("#{}",\n'.format(getid(t))
			htmlout += '\t{},\n'.format(getid(t))
			htmlout += '\t\t{\n'
			htmlout += '\t\tdateFormat: "%B %Y",\n'
			htmlout += '\t\tlabelFormat: "%Y/%m/%d %H:%M:%S",\n'
			htmlout += '\t\tcolor: "#2971B0",\n'
			htmlout += '\t\tbackground: "{}",\n'.format(settings[t]["background"])
			htmlout += '\t\tshowLabels: true,\n'
			htmlout += '\t\tlabelFormat: "%Y", \n'
			htmlout += '\t\taddNow: true,\n'
			htmlout += '\t\twidth: 1000,\n'
			htmlout += '\t\theight: 100,\n'
			htmlout += '\t\tradius: 10,\n'
			htmlout += '\t\tlineWidth: 4,\n'
			htmlout += '\t\thorizontalLayout: true,\n'
			htmlout += '\t\t}\n'
			htmlout += '\t);\n\n'

		htmlout += '</script>'
		
		with open(args.outputfile,"w") as o:
			o.write(htmlout)

		# Make svg timeline
		print ("\n\n==SVG==\n\n")
		for t in timelines:
			svgdict = settings[t]
			svgdict["callouts"] = []
			svgdict["eras"] = []
			for event in timelines[t]:
				if "color" not in event:
					event["color"] = "black"
				svgdict["callouts"].append([
					event["name"],
					event["date"],
					event["color"]				
					])
			for era in eras[t]:
				if "color" not in era:
					era["color"] = "#C0C0FF"
				svgdict["eras"].append([
					era["name"],
					era["start"],
					era["end"],
					era["color"],
					])
			svgjson = json.dumps(svgdict, indent=4)
			svgfile = "{}.svg".format(slugify(t)[:20])
			print (svgfile)
			print (svgjson)
			timeline = make_timeline.Timeline(svgjson, svgfile)
			timeline.build()
			timeline.save(svgfile)
			# use inkscape to fix svg area
			cmd = 'inkscape --export-type="svg" --export-area-drawing -o {} {}'.format(svgfile,svgfile)
			print (cmd)
			subprocess.call(cmd, shell=True)


		











