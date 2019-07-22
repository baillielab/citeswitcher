---
csl: /Users/jkb/Dropbox/3_scripts_and_programs/citeswitcher/sup-files/csl/nature.csl
bibliography: /Users/jkb/Dropbox/3_scripts_and_programs/citeswitcher/README.bib
...
# Citeswitcher

For simple collaboration in academic writing. Cite papers as PMID, DOI, bib reference, or whole citation, and citeswitcher will find them and make a common bib file for you.

# Example

You can cite papers in a plain text document like this [@HallNetworkAnalysisReveals2014] or like this [@HallNetworkAnalysisReveals2014] or (if you share a .bib file) like this [@HallNetworkAnalysisReveals2014] or (if you really need to) like this [@HallNetworkanalysisreveals2014]. Citeswitcher will find these and format them properly as citations in any format you choose.

# Intended use:
- Feed in a .md or .tex document
- references in square [] or curly {} brackets will be searched for:
-- PMID
-- doi
-- bibtext
-- whole ciation

The default or specified .bib file will then be searched for the relevant citations.
Citations will be replaced with a citaion in either .md or .tex format.

Any citations not present in the master .bib file will be downloaded from pubmed and added to the supplementary.bib file

Two new files will be written:
1. a new .md or .tex file with correct citations
2. a local .bib file containing all cited reference details.

# Examples

python fix_citations.py -f test/eme.md

# Requirements
Python 3.7
- bibtexparser
- requests

pandoc
pandoc-citeproc

pdflatex

