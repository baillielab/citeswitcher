---
title: Citeswitcher
author: Kenneth Baillie
date: \today
...
 

# Citeswitcher

Intended use:
- Feed in a .md or .tex document
- references in square [] or curly {} brackets will be searched for:
-- PMID
-- doi

The default or specified .bib file will then be searched for the relevant citations. 
Citations will be replaced with a citaion in either .md or .tex format.

Any citations not present in the master .bib file will be downloaded from pubmed and added to the supplementary.bib file

Two new files will be written:
1. a new .md or .tex file with correct citations
2. a local .bib file containing all cited reference details. 

# Examples

python fix_citations.py -f test/eme.md

# Requirements
Python 2.7 
bibtexparser

pandoc
pandoc-citeproc

latex
