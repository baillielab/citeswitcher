<!--
TODO

ADD stripreferences function to remove all confirmed references

-->


# Citeswitcher

**A new project, [manubot](https://github.com/manubot/manubot), performs some of the functions of citeswitcher and is actively developed by professionals, in contrast to my amateur hacking on this project. Worth a look before you get too deep into citeswitcher...**

For simple collaboration in academic writing. Cite papers as PMID, DOI, bib reference, or whole citation, and citeswitcher will find them and make a common bib file for you, and call pandoc to make a fully referenced pdf, docx or html file.

# TL;DR

`python fixcitations.py -f sup-files/test.md -p pdf`

You can cite papers in a plain text document like this [PMID: 24465370] or like this [DOI: 10.1371/journal.pone.0081229] or (if you share a .bib file) like this [@HallNetworkAnalysisReveals2014] or (if you really need to) like this [Hall, David P., Ian J. C. MacCormick, Alex T. Phythian-Adams, Nina M. Rzechorzek, David Hope-Jones, Sorrel Cosens, Stewart Jackson, et al. “Network Analysis Reveals Distinct Clinical Syndromes Underlying Acute Mountain Sickness.” PLoS ONE 9, no. 1 (January 22, 2014): e81229].

Citeswitcher will find these and format them properly as citations in the chosen output format.

# Requirements

## OS
- works on all linux and mac systems tested so far
- not tested on windows

## Essential

For basic function of fixcitations.py:

- python 3
- pyparsing 2.4.7
- requests-2.24.0
- biopython
- oyaml
- bibtexparser is also required but it is included in this git repository (because I was concerned that it may not be well maintained)

To install these:
pip install pyparsing
pip install requests
pip install biopython
pip install oyaml

## For pdf, docx and html outputs

You'll want to have pandoc properly installed with crossref so that you can automatically generate printable documents.

pandoc
pandoc-crossref

This can be done with conda, or pip.
`pip install pandoc-crossref`

## Optional additional requirements

Once you have the above installed, you're good to go. Additional dependencies for extra functions are listed at the end.

# Installation

Download or clone the git repository into a folder on your hard disk. Then run the scripts using python from the command line.
To make it easier to use, I use the unix `alias` command in my `~/.bash_profile` to automatically map common commands, like this:

`alias fix="/miniconda3/envs/thisenv/bin/python <full-path-to-git-repository>/citeswitcher/fixcitations.py"`

`alias hack="/miniconda3/envs/thisenv/bin/python <full-path-to-git-repository>/citeswitcher/hackmd.py"`

## Config

Edit the `config.json` file to include your own paths before running. Key elements:

- email: required for using the Pubmed API - enter your own email address here.
- filepaths: these need to be on you computer somewhere.

# Use

## fixcitations.py

Switches citation formats, finds missing citations online, and calls pandoc to create a fully referenced output file.

The default or specified .bib file will then be searched for the relevant citations.
Citations will be replaced with a citation in either .md or .tex format.

Any citations not present in the master .bib file will be downloaded from pubmed and added to a .bib file

Two new files will be written:

1. a new .md or .tex file with correct citations
2. a local .bib file containing all cited reference details.

Optionally, additional files can be written:

3. pandoc outputs in .docx, .pdf or .html formats

### Input files

Required:

- a .md or .tex document with references in square [] or curly {} or round brackets which will be searched for PMIDs, DOIs, bibtex refs, and whole citations

Optional:

- *bib file*.
	A bibtex file to use. If not specified, a new `cs.bib` file will be created in the same directory as your input file. You can also specify a global bib file for keeping a shared bank of references.
- *yaml file*.
	A yaml file called <filestem>.yaml will be automatically created in the same directory as your file. If it already exists, it will be augmented to add csl and bib files if needed.
	YAML is taken according to the following hierarchy: in-file YAML, <filestem>.yaml, other yaml files. All will be copied into the <filestem>.yaml file.
- *csl file*.
	You can specify a csl file using yaml.
- *replace.json file*.
	A file called `replace.json` will be read as a dictionary of text strings to find:replace. If you just create this file and leave it in the same directory as the input file, it will be read. An alternative is to specify text to be replaced in a YAML dict format.

The default (in `config.json`) or specified .bib file will be searched for the relevant citations.
Citations will be replaced with a citaion in either .md or .tex or PMID format.

Any citations not found in the .bib files will be downloaded from pubmed and added to the supplementary.bib file

### Output files

Two outputs are specified.

`-o` determines the CITATION FORMAT and is used for switching between different citation formats e.g. PMID==>DOI or PMID==>BIB

`-p` determines the PANDOC OUTPUT FORMAT and provides a quick workflow for generating fully formatted documents

If `-p` is specified, `-o` must be md. If both `-p` and `-o` are specified and `-o` is not md, this indicates that the user's expectation is different from what's going to happen so the script exits. Fix this by only specifying one of these options.

## Options

| arg | long arg| action|
|----|--------------|---------------------------|
| -b | --bibfile | bibfile |
| -y | --yaml | specify a yaml file to use as a source for yaml; use "normal" or "fancy" to copy from sup-files/yaml|
| -o | --outputstyle | choices=md,markdown,tex,latex,pubmed,pmid,inline; output references format |
| -p | --pandoc_outputs | append as many pandoc formats as you want: pdf docx html txt md tex |
| -l | --localbibonly | use only local bib file |
| -d | --outputsubdir | outputdir - always a subdir of the working directory || -img | --imagedir | imagedirectoryname |
| -i | --include | do NOT include files |
| -m | --messy | disable clean up of intermediate files |
| -mf | --move_figures | move all figures to the end and create captions section for submission to journal |
| -pm | --pandoc_mermaid | use pandoc-mermaid-filter |
| -lt | --latex_template | a latex template to be passed to pandoc |
| -wt | --word_template | a word template to be passed to pandoc |
| -ptp | --pathtopandoc | specify a particular path to pandoc if desired') # e.g. a letter forma |
| -redact | --redact | redact between <!-- STARTREDACT --> <!-- ENDREDACT --> tags |
| -s | --stripcomments | stripcomments in html format |
| -v | --verbose | verbose |
| -x | --xelatex | use xelatex in pandoc build |
| -w | --wholereference | try to match whole references. |
| -ch | --chaptermode | use pandoc --top-level-division=chapter |
| -svg | --convert_svg | convert svg images to pdf - replaces any pdf files with the same name |
| -ui | --uncomment_images | include images commented out using html syntax <!--![]()\{\} --> |
| -flc | --force_lowercase_citations | force all citation references into lowercase |

## hackmd.py

Uses `wget` to download your **published** files from hackmd.io and runs fixcitations.py. This makes working collaboratively with hackmd very easy, but takes a bit of manual setup:

1. Specify the filepath to a directory (`hackmddir`) in `config.json`. This should be an empty directory.
2. *Within* your new `hackmddir`, create a folder with the exact name of your [hackmd.io](hackmd.io) username or group name.
3. Create a hackmd file and publish it. It's easier if you choose the url name for it, such as "mytestfile"
4. Create yet another folder with the same name as your hackmd file. So your directory structure would be like this:
```
- my_hackmd_dir
|	+-- my_username
|		+-- mytestfile

```
5. Run `python hackmd.py`

If this works, your file should download into the mytestfile directory, and you should see the following files:
```
- my_hackmd_dir
|	+-- my_username
|		+-- mytestfile
|			+-- mytestfile.md
|			+-- mytestfie.citemd.md
|			+-- mytestfie.citemd.pdf
|			+-- mytestfie.citemd.docx
|			+-- mytestfie.citemd.html
```

hackmd.py - download hackmd files either by combining elements from a directory tree or a single file, then run fixcitations.py

# Additional utilities

The following additional scripts are included to provide some useful functions with varying levels of documentation:
- *include.py* - include files recursively according to a variety of syntaxes
- *put_pmid_in_bibtex.py* - search for PMIDs for every element in a .bib file. Designed to run as a cron job to keep your bib file up to date.
- *abbreviation_finder.py* - print abbreviations in a document (so that you can explain them)
- *wordcount.py* - function unknown at this time
- *svg2pdf.py* - batch convert a directory (-d) or file (-f) using cairosvg, which seems to perform better than imagemagick.
- *format_author_list.py* - this is a work in progress but it is intended to create author lists in a variety of formats from a spreadsheet in xlsx or csv format

## Additional dependencies for optional extra functions

### For including excel files in markdown
- pandas 1.1
- tabulate-0.8.7

`pip install pandas`
`pip install tabulate`

### For converting from svg to pdf

- inkscape

### To manage your own bibtex database

There are lots of bibtex database software tools available. The following work very well:
- zotero
- better bibtex (with automatic export feature)

# Intended use:
- Feed in a .md or .tex document
- references in square [] or curly {} brackets will be searched for:
-- PMID
-- doi
-- bibtext
-- whole ciation

### For making Gantt, flowcharts etc

- mermaids.cli (https://github.com/mermaid-js/mermaid-cli)
Mermaids is a bit of a hassle to install and needs to be in you PATH variable for this to work. But if it is in there, you can use the -pm option to call mermaids and automatically generate nice figures.

This installation procedure worked for me on a mac:

#### Install node:
https://www.npmjs.com/get-npm
#### Install pandiff:
https://github.com/davidar/pandiff

#### Install mermaid:
`npm install mermaid.cli`
Test it: 
`./node_modules/.bin/mmdc -h`
Add mermaid to ~/.bash_profile:
export MERMAID_BIN=/Users/<username>/node_modules/.bin/mmdc
export PATH=$PATH:/Users/<username>/node_modules/.bin
~




