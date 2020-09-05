# Citeswitcher

For simple collaboration in academic writing. Cite papers as PMID, DOI, bib reference, or whole citation, and citeswitcher will find them and make a common bib file for you.

# Example

`python fixcitations.py -f test/eme.md -p pdf`

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
- bibtexparser is also required but it is included in this git repository (because I was concerned that it may not be well maintained)

To install these:
`pip install pyparsing`
`pip install requests`
`pip install biopython`

## For pdf, docx and html outputs

You'll want to have pandoc properly installed with crossref and citeproc so that you can automatically generate printable documents.

pandoc
pandoc-crossref
pandoc-citeproc

This can be done with conda, or pip.
`pip install pandoc-crossref`
`pip install pandoc-citeproc`

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

The default or specified .bib file will then be searched for the relevant citations.
Citations will be replaced with a citaion in either .md or .tex format.


Any citations not present in the master .bib file will be downloaded from pubmed and added to a .bib file

Two new files will be written:
1. a new .md or .tex file with correct citations
2. a local .bib file containing all cited reference details.

### Input files

Required:

- a plain text file

Optional:

- bib file.
	A bibtex file to use. If not specified, a new `myfile.bib` file will be created in the same directory as your input file. You can also specify a global bib file for keeping a shared bank of references.
- csl file
	You can add a csl file to the same directory as the md file, and then direct pandoc to it using yaml. That simply means putting this exact text at the top of your md file:
	```
	---
	csl: myfile.csl
	---
	```
	or by specifying a separate YAML file
- yaml file
	You can specify an external YAML file
- replace.json file

Intended use:
Feed in a .md or .tex document with references in square [] or curly {} or round brackets will be searched for PMID

If there is a mixture of reference types in one set of parentheses, ids with the format PMID:NNNNN will be handled correctly but spaces etc will be ignored.

The default or specified .bib file will then be searched for the relevant citations.
Citations will be replaced with a citaion in either .md or .tex or PMID format.

Any citations not present in the master .bib file will be downloaded from pubmed and added to the supplementary.bib file

Two new files will be written:
1. a new .md or .tex file with correct citations
2. a local .bib file containing all cited reference details.

### OUTPUT FILE ARGUMENTS LOGIC

Two outputs are specified.
-o determines the CITATION FORMAT and is intended to be the primary use case - switching between different citation formats
-p determines the PANDOC OUTPUT FORMAT and provides a quick workflow for generating fully formatted documents

If -p is specified, -o must be md.
If both -p is specified an -o is not unset or md, this indicates that the user's expectation is different ==> so the script exits.




# Additional utilities

The following additional scripts are included to provide some useful functions:
hackmd.py - download hackmd files either by combining elements from a directory tree or a single file, then run fixcitations.py
include.py - include files recursively according to a variety of syntaxes
put_pmid_in_bibtex.py - search for PMIDs for every element in a .bib file. Designed to run as a cron job to keep your bib file up to date.
abbreviation_finder.py - print abbreviations in a document (so that you can explain them)
wordcount.py - function unknown at this time
svg2pdf.py - batch convert a directory (-d) or file (-f) using cairosvg, which seems to perform better than imagemagick.

## Additional dependencies for optional extra functions

### For including excel files in markdown
- pandas 1.1
- tabulate-0.8.7

`pip install pandas`
`pip install tabulate`

### For converting from svg to pdf

- cairosvg (`pip install cairosvg`)

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

- mermaids.cli (https://github.com/mermaidjs/mermaid.cli)
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




