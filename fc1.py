def findcitation(info, infotype='pmid', additionalinfo='', force_search=False):
    """
    Searches for a citation in full_bibdat or online.
    """
    global no_user_response_count
    info = info.strip()
    if infotype == 'pmid':
        if not force_search:
            entry = full_bibdat.get_entry_by_pmid(info)
            if entry:
                return entry
            else:
                print(f"PMID {info} not found in full_bibdat.")
        # Search online
        pub_entries = p2b([info])
        if pub_entries:
            new_entry = pub_entries[0]
            if new_entry:
                full_bibdat.add_entry(new_entry)
                return new_entry
        print(f"PMID {info} not found online.")
        return None
    elif infotype == 'doi':
        if not force_search:
            entry = full_bibdat.get_entry_by_doi(info)
            if entry:
                return entry
            else:
                print(f"DOI {info} not found in full_bibdat.")
        # Search online
        pmids = search_pubmed(info, "doi")
        if pmids:
            pub_entries = p2b(pmids)
            if pub_entries:
                new_entry = pub_entries[0]
                if new_entry:
                    full_bibdat.add_entry(new_entry)
                    return new_entry
        print(f"DOI {info} not found online.")
        return None
    elif infotype == 'title':
        if no_user_response_count > 2:
            print("No user response to previous 3 queries. Running in silent mode.")
            return None
        print(f"Searching PubMed for title: {info}")
        pmids = search_pubmed(info, "title")
        if len(pmids) == 1:
            pmid = pmids[0]
            pub_entries = p2b([pmid])
            if pub_entries:
                new_entry = pub_entries[0]
                question = "--------------\n\
New citation (PMID:{}) found in Pubmed. Please check that input is the same as the found citation for this reference block: \n\
\n{}\n\n{}\n{}\n\n{}\n\n\
Enter y/n within 10 seconds".format(
                        pmid,
                        additionalinfo,
                        "{:>12}:    {}".format('Input Title', info),
                        "{:>12}:    {}".format('Found Title', pubent[0]['Title']),
                        '\n'.join( ["{:>12}:    {}".format(x,pubent[0][x]) for x in pubent[0] if x != "Title"])
                        )
                    #q = input (question)
                    print (question)
                    i,o,e = select.select([sys.stdin],[],[],10) # 10 second timeout
                    if i==[] and o==[] and e==[]:
                        no_user_response_count += 1
                    if i:
                        q = sys.stdin.readline().strip()
                        q = q.strip().upper()
                        if q == "Y":
                            print ('--confirmed--')
                            full_bibdat.add_entry(new_entry)
                            return new_entry
        print(f"Title {info} not found online.")



