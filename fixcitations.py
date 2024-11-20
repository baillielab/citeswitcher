#!/usr/bin/env python
# encoding: utf-8

import os
import io
import sys
import copy
import argparse
import re
import requests
import calendar
import xml.etree.ElementTree as ET

# Add necessary paths to sys.path
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath, 'python-bibtexparser-master/'))
sys.path.append(os.path.join(scriptpath, 'dependencies/'))

import bibtexparser
import oyaml as yaml

import bib
bib.init()
import citefunctions
config = citefunctions.getconfig()

# Patch for yaml error in Pythonista on iPad
import collections
if not hasattr(collections, 'Hashable'):
    import collections.abc
    collections.Hashable = collections.abc.Hashable

# ----------------- Utility Functions -----------------

def getconfig(cfgfile="null"):
    if cfgfile == "null":
        for cfgname in ['config_local.json', "config.json"]:
            cfgfile = os.path.join(scriptpath, cfgname)
            if os.path.exists(cfgfile):
                break
    with open(cfgfile) as json_data_file:
        data = json.load(json_data_file)
    for item in data:
        if isinstance(data[item], str) and data[item].startswith("{{scriptpath}}"):
            data[item] = os.path.join(scriptpath, data[item].replace("{{scriptpath}}", ""))
    if "@" not in data['email']:
        print("No email in config file: {}".format(cfgfile))
        sys.exit()
    return data

def make_unicode(inputstring):
    nbspace = re.compile(u"\N{NO-BREAK SPACE}", re.IGNORECASE)
    if not isinstance(inputstring, str):
        inputstring = inputstring.decode('utf-8')
    inputstring = nbspace.sub(" ", inputstring)
    inputstring = inputstring.replace(u"\u2003", " ")  # non-breaking space
    inputstring = inputstring.replace(u"\u2009", " ")  # thin space
    inputstring = inputstring.replace(u"\ufeff", "")   # BOM character
    return inputstring

def getyaml(text):
    h, r = readheader(text)
    if h.startswith("---"):
        h = h[4:-3]  # strip yaml identifiers
    yml = yaml.load(h, Loader=yaml.Loader)
    if yml:
        return yml
    else:
        return {}

def mergeyaml(priority_yaml, extra_yaml):
    # Keep the contents of the first yaml if there is a conflict
    if not priority_yaml:
        priority_yaml = {}
    if extra_yaml:
        for item in extra_yaml:
            if item not in priority_yaml or priority_yaml[item] is None:
                priority_yaml[item] = extra_yaml[item]
    return priority_yaml

def readheader(filecontents):
    '''
        Read a valid markdown yaml header
        Input is full text of file
        Returns list of header items + full text of remainder
    '''
    t = filecontents.strip().replace('\r', '\n')
    h = ""
    remainder = filecontents
    if t.startswith('---'):
        h_match = re.match(r'---[\s\S]+?---', t)
        if h_match:
            h = h_match.group(0)
            remainder = filecontents.replace(h, '', 1)
    return h, remainder

def read_bib_files(bibfiles):
    bfs = ""
    for bibfile in bibfiles:
        if os.path.exists(bibfile):
            try:
                with open(bibfile, encoding="utf-8") as bf:
                    bfs += bf.read()
            except:
                with open(bibfile, encoding="latin1") as bf:
                    bfs += bf.read()
        else:
            print("File does not exist: {}".format(bibfile))
    try:
        parser = bibtexparser.bparser.BibTexParser(common_strings=True, homogenize_fields=True, interpolate_strings=False)
        return parser.parse(bfs, partial=False)
    except:
        return bibtexparser.bibdatabase.BibDatabase()

def search_pubmed(search_string, restrictfields=""):
    """
    Return a list of PMIDs for articles found by a search string.

    Args:
        search_string (str): The search term for querying PubMed.
        restrictfields (str): Optional field restriction for the search (e.g., 'title').

    Returns:
        list: A list of PMIDs for matching articles.
    """
    if not search_string or not search_string.strip():
        print("Empty search string provided.")
        return []

    max_returns = 500
    print(f"Searching PubMed for: {search_string}")

    # Modify search string for DOI searches
    if restrictfields.lower() == 'doi':
        search_string = f"{search_string}[AID]"
        restrictfields = ""  # Clear restrictfields

    # Build the query parameters
    params = {
        'db': 'pubmed',
        'term': search_string,
        'retmax': max_returns,
        'retmode': 'xml',
        'email': config.get('email', 'your_email@example.com'),  # Use email from config
        # 'api_key': 'your_api_key',        # Uncomment and set if you have an API key
    }
    # Do not include 'field' parameter if empty
    if restrictfields:
        params['field'] = restrictfields

    # Make the GET request to the ESearch API
    try:
        response = requests.get(
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
            params=params
        )
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"PubMed search failure: {e}")
        return []

    # Parse the XML response
    try:
        root = ET.fromstring(response.text)
    except ET.ParseError as e:
        print(f"Error parsing PubMed response: {e}")
        return []

    # Extract PMIDs from the response
    idlist = root.find('IdList')
    pmid_list = [id_elem.text for id_elem in idlist.findall('Id')] if idlist is not None else []

    foundcount = len(pmid_list)
    print(f"{foundcount} records found.")

    return pmid_list

def p2b(pmidlist):
    ''' Convert PMID to BibTeX entry '''
    if not isinstance(pmidlist, list):
        pmidlist = [str(pmidlist)]

    # Fetch XML data from Entrez.
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    url = '{}?db=pubmed&id={}&rettype=abstract'.format(efetch, ','.join(pmidlist))
    try:
        r = requests.get(url)
        r.raise_for_status()
    except requests.exceptions.RequestException:
        return []

    # Parse the XML using xml.etree.ElementTree
    bibout = []
    try:
        root = ET.fromstring(r.text)
    except ET.ParseError:
        return []

    for PubmedArticle in root.findall('.//PubmedArticle'):
        PMID = PubmedArticle.find('./MedlineCitation/PMID')
        ISSN = PubmedArticle.find('./MedlineCitation/Article/Journal/ISSN')
        Volume = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/Volume')
        Issue = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/Issue')
        Year = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
        Month = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Month')
        Title = PubmedArticle.find('./MedlineCitation/Article/Journal/Title')
        ArticleTitle = PubmedArticle.find('./MedlineCitation/Article/ArticleTitle')
        MedlinePgn = PubmedArticle.find('./MedlineCitation/Article/Pagination/MedlinePgn')

        # Additional IDs
        PMCID = None
        DOI = None
        theseids = PubmedArticle.findall('./PubmedData/ArticleIdList/ArticleId')
        for thisid in theseids:
            if thisid.attrib.get('IdType') == 'pmc':
                PMCID = thisid
            elif thisid.attrib.get('IdType') == 'doi':
                DOI = thisid

        # Format author list
        authors = []
        for Author in PubmedArticle.findall('.//Author'):
            LastName = Author.find('LastName')
            ForeName = Author.find('ForeName')
            if LastName is not None and ForeName is not None:
                authors.append('{}, {}'.format(LastName.text, ForeName.text))

        # Use InvestigatorList if no authors found
        if not authors:
            for Investigator in PubmedArticle.findall('.//Investigator'):
                LastName = Investigator.find('LastName')
                ForeName = Investigator.find('ForeName')
                if LastName is not None and ForeName is not None:
                    authors.append('{}, {}'.format(LastName.text, ForeName.text))

        # Handle missing Year and Month
        if Year is None:
            MedlineDate = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate')
            if MedlineDate is not None and MedlineDate.text:
                Year = MedlineDate.text[:4]
                try:
                    month_abbr = MedlineDate.text[5:8]
                    Month = '{:02d}'.format(list(calendar.month_abbr).index(month_abbr))
                except ValueError:
                    Month = None
            else:
                Year = 'Unknown'
                Month = None
        else:
            Year = Year.text
            if Month is not None:
                Month = Month.text

        # Create a unique ID
        bib_entry = {}
        authorname = authors[0].split(',')[0] if authors else ''
        try:
            titlewords = [x for x in ArticleTitle.text.split(' ') if len(x) > 3]
        except AttributeError:
            print("PUBMED ERROR - no article title for PMID:{}.".format(PMID.text if PMID else 'Unknown'))
            continue

        if len(titlewords) > 2:
            titlestring = ''.join(titlewords[:3])
        elif titlewords:
            titlestring = ''.join(titlewords)
        else:
            titlestring = ''

        if not (authorname + titlestring):
            titlestring = "PMID{}_".format(PMID.text if PMID else 'Unknown')

        new_id = '{}{}{}'.format(authorname, titlestring, Year).lower()
        new_id = re.sub(r'\W+', '', new_id)
        bib_entry["ID"] = new_id

        # Populate bibtex fields
        bib_entry["Author"] = ' and '.join(authors)
        bib_entry["Title"] = ArticleTitle.text if ArticleTitle is not None else ''
        bib_entry["Journal"] = Title.text if Title is not None else ''
        bib_entry["Year"] = Year
        if Volume is not None:
            bib_entry["Volume"] = Volume.text
        if Issue is not None:
            bib_entry["Number"] = Issue.text
        if MedlinePgn is not None:
            bib_entry["Pages"] = MedlinePgn.text
        if Month is not None:
            bib_entry["Month"] = Month
        if PMCID is not None:
            bib_entry["pmcid"] = PMCID.text
        if DOI is not None:
            bib_entry["doi"] = DOI.text
        if ISSN is not None:
            bib_entry["ISSN"] = ISSN.text
        if PMID is not None:
            bib_entry["pmid"] = PMID.text

        # Append to output list
        bibout.append(bib_entry)

    return bibout

# ----------------- Main Function -----------------

def main(
    filepath,
    bibfile,
    yaml_file,
    wholereference,
    outputstyle,
    overwrite,
    localbibonly,
    force_lowercase_citations,
    outputsubdir,
    verbose
):
    # Determine source path and filenames
    sourcepath, filename = os.path.split(filepath)
    outpath = os.path.join(sourcepath, outputsubdir)
    if outpath != '':
        citefunctions.check_dir(outpath)
    filestem = '.'.join(filename.split('.')[:-1])

    # Read input file
    with io.open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    # Read YAML according to the hierarchy: infile, local.yaml, other
    if not yaml_file.endswith(".yaml"):
        yaml_file = yaml_file + ".yaml"
    yamlfile = os.path.join(sourcepath, filestem + ".yaml")
    infileyaml = getyaml(text)  # READ FROM INPUT FILE FIRST
    workingyaml = copy.copy(infileyaml)
    if os.path.exists(yamlfile):
        with open(yamlfile) as f:
            workingyaml = mergeyaml(workingyaml, getyaml(f.read()))
    if yaml_file in os.listdir(config['yamldir']):
        with open(os.path.join(config['yamldir'], yaml_file)) as f:
            workingyaml = mergeyaml(workingyaml, getyaml(f.read()))
    # CSL - hierarchy - yaml-specified, sup-files
    if 'csl' in workingyaml:
        if not workingyaml['csl'].endswith('.csl'):
            workingyaml['csl'] = workingyaml['csl'] + '.csl'
        if workingyaml['csl'].startswith("http"):
            pass
        elif not os.path.exists(workingyaml['csl']):
            # Try to find it in the sup-files folder
            newcsl = os.path.relpath(os.path.join(config['csldir'], os.path.split(workingyaml['csl'])[-1]))
            if os.path.exists(newcsl):
                workingyaml['csl'] = newcsl
    else:
        workingyaml['csl'] = os.path.relpath(os.path.join(config["csldir"], config["csldefault"]))  # HARD OVERWRITE
    # BIB - read them all and copy into one local version
    bibfile = os.path.abspath(os.path.expanduser(bibfile))
    if 'bibliography' in workingyaml.keys():
        print('Using YAML-specified bib:', workingyaml['bibliography'])
    else:
        workingyaml['bibliography'] = config['default_localbib']  # HARD OVERWRITE
    print("Using {} as bibout".format(workingyaml['bibliography']))

    if verbose:
        print("Filepath:", filepath)
    input_file_extension = filepath.split('.')[-1]
    if filepath.endswith((".md", ".txt", ".qmd")):
        if outputstyle == 'null':
            outputstyle = 'md'
    elif filepath.endswith(".tex"):
        if outputstyle == 'null':
            outputstyle = 'tex'

    # Name output file
    if overwrite:
        print("Overwriting original file")
        citelabel = "."
    elif outputstyle in ('pubmed', 'pmid'):
        print("Outputstyle = pmid")
        citelabel = ".citepmid."
    elif outputstyle in ('markdown', 'md'):
        print("Outputstyle = md")
        citelabel = ".citemd."
    elif outputstyle in ('latex', 'tex'):
        print("Outputstyle = tex")
        citelabel = ".citetex."
    elif outputstyle == 'inline':
        print("Outputstyle = inline")
        citelabel = "."
    else:
        citelabel = "."

    outputfile = os.path.join(outpath, filestem + citelabel + input_file_extension)

    text = make_unicode(text)
    text = readheader(text)[1]
    if verbose:
        print("Read input file:", filepath)

    # Prepare bibliography
    localbibpath = os.path.join(sourcepath, workingyaml['bibliography'])
    bib.db = read_bib_files([localbibpath])
    if localbibonly:
        print("\n*** Reading local bib file only: {} ***\n".format(localbibpath))
        bib.full_bibdat = bib.db
    else:
        if verbose:
            print("Reading bibfiles:", bibfile, localbibpath)
        bib.full_bibdat = read_bib_files([bibfile, localbibpath])
    if force_lowercase_citations:
        print("Forcing lowercase citations")
        bib.id_to_lower()
    bib.make_alt_dicts()

    # Replace the ids in the text with the outputstyle
    text = citefunctions.replace_blocks(text, outputstyle, use_whole=wholereference, flc=force_lowercase_citations)

    # Save bibliography
    print('\nSaving bibliography for this file here:', localbibpath)
    outbib = bibtexparser.dumps(bib.db)
    outbib = make_unicode(outbib)
    with open(localbibpath, "w", encoding="utf-8") as bf:
        bf.write(outbib)

    # Save new text file
    with io.open(outputfile, 'w', encoding='utf-8') as file:
        if outputstyle == "md":
            file.write('---\n{}\n---'.format(yaml.dump(workingyaml)).replace("\n\n", "\n"))
        file.write(text + "\n\n")
        print("Outputfile:", outputfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Essential arguments
    parser.add_argument("-f", '--filepath', help='Path to the input file', default="../genomicc-manuscript/manuscript.tex")
    # Additional files to specify
    parser.add_argument('-b', '--bibfile', default=config['default_bibfile'], help='BibTeX file')
    parser.add_argument('-y', '--yaml', default='auto', help='YAML file to use; use "normal" or "fancy" to use templates')
    # Other options
    parser.add_argument('-w', '--wholereference', action="store_true", default=False, help='Try to match whole references.')
    parser.add_argument('-o', '--outputstyle', type=str, choices=['md', 'markdown', 'tex', 'latex', 'pubmed', 'pmid', 'inline'], default='null', help='Output references format')
    parser.add_argument('-ow', '--overwrite', action="store_true", default=False, help='Overwrite input file with new version')
    parser.add_argument('-l', '--localbibonly', action="store_true", default=False, help='Use only local BibTeX file')
    parser.add_argument('-flc', '--force_lowercase_citations', action="store_true", default=False, help='Force all citation references into lowercase')
    parser.add_argument('-d', '--outputsubdir', default=config['outputsubdirname'], help='Output directory (subdirectory of the working directory)')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='Verbose output')
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(
        filepath=os.path.abspath(os.path.expanduser(args.filepath)),
        bibfile=args.bibfile,
        yaml_file=args.yaml,
        wholereference=args.wholereference,
        outputstyle=args.outputstyle,
        overwrite=args.overwrite,
        localbibonly=args.localbibonly,
        force_lowercase_citations=args.force_lowercase_citations,
        outputsubdir=args.outputsubdir,
        verbose=args.verbose
    )
