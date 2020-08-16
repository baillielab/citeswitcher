# Citeswitcher

For simple collaboration in academic writing. Cite papers as PMID, DOI, bib reference, or whole citation, and citeswitcher will find them and make a common bib file for you.

# Example

`python fixcitations.py -f test/eme.md -p pdf`

You can cite papers in a plain text document like this [PMID: 24465370] or like this [DOI: 10.1371/journal.pone.0081229] or (if you share a .bib file) like this [@HallNetworkAnalysisReveals2014] or (if you really need to) like this [Hall, David P., Ian J. C. MacCormick, Alex T. Phythian-Adams, Nina M. Rzechorzek, David Hope-Jones, Sorrel Cosens, Stewart Jackson, et al. “Network Analysis Reveals Distinct Clinical Syndromes Underlying Acute Mountain Sickness.” PLoS ONE 9, no. 1 (January 22, 2014): e81229].

Citeswitcher will find these and format them properly as citations in the chosen output format.

# Requirements

## OS
- works on all linux and max systems tested so far
- not tested on windows

## For basic function of fixcitations.py

Python 3.7
- pyparsing 2.4.7 
- pandas 1.1 
- tabulate-0.8.7 
- requests-2.24.0 
- biopython 
- bibtexparser is also required but it is included in this git repository (because it may not be well maintained but it works)

To install these: 
pip install pyparsing
pip install pandas
pip install tabulate
pip install requests
pip install biopython

## For pdf, docx and html outputs

pandoc
pandoc-crossref
pandoc-citeproc

## For converting from svg to pdf
- cairosvg (`pip install cairosvg`)

## To manage your own bibtex database

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

The default or specified .bib file will then be searched for the relevant citations.
Citations will be replaced with a citaion in either .md or .tex format.

Any citations not present in the master .bib file will be downloaded from pubmed and added to a .bib file

Two new files will be written:
1. a new .md or .tex file with correct citations
2. a local .bib file containing all cited reference details.

# Input files

Required:

- a plain text file

Optional:

- bib file
- yaml file
- csl file
- [to do] find/replace json file

## BIB file

A bibtex file to use. If not specified, a new `myfile.bib` file will be created in the same directory as your input file.
You can also specify a global bib file for keeping a shared bank of references.

## CSL file

You can add a csl file to the same directory as the md file, and then direct pandoc to it using yaml. That simply means putting this exact text at the top of your md file:
```
---
csl: myfile.csl
---
```
or by specifying a separate YAML file


# Config

Edit the config.json file before running. Key elements:
email: required for using the Pubmed API - enter your own email address here.

# Shortcuts

- edit your `~/.bash_profile` file:

alias fix="python <path-to-fixcitations-script>"
alias hack="python <path-to-get_hackmd-script>"


# Additional utilities

The following additional scripts are included to provide some useful functions:
hackmd.py - download hackmd files either by combining elements from a directory tree or a single file, then run fixcitations.py
include.py - include files recursively according to a variety of syntaxes
put_pmid_in_bibtex.py - search for PMIDs for every element in a .bib file. Designed to run as a cron job to keep your bib file up to date.
abbreviation_finder.py - print abbreviations in a document (so that you can explain them)
wordcount.py - function unknown at this time
svg2pdf.py - batch convert a directory (-d) or file (-f) using cairosvg, which seems to perform better than imagemagick.

# To do:
- move "custom_find_replace" out of config file into an optional file in same dir as the input .md file, just like .bib & .yaml files
- [to do] automatic installation








