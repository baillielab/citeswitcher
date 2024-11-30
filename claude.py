from dataclasses import dataclass
from typing import Dict, List, Optional, Set
import copy
from collections import OrderedDict

@dataclass
class Citation:
    id: str
    title: str
    author: str
    year: str
    journal: str
    volume: Optional[str] = None
    number: Optional[str] = None
    pages: Optional[str] = None
    month: Optional[str] = None
    pmid: Optional[str] = None
    doi: Optional[str] = None
    pmcid: Optional[str] = None
    issn: Optional[str] = None
    entry_type: str = "article"

class CitationDatabase:
    def __init__(self):
        self._entries: Dict[str, Citation] = OrderedDict()
        self._pmid_index: Dict[str, str] = {}  # PMID -> ID mapping
        self._doi_index: Dict[str, str] = {}   # DOI -> ID mapping
        self._id_changes: Dict[str, str] = {}  # old_id -> new_id mapping

    def add_entry(self, entry: Citation) -> None:
        """Add or update a citation entry."""
        existing_id = None
        
        # Check for existing entry by DOI or PMID
        if entry.doi and entry.doi in self._doi_index:
            existing_id = self._doi_index[entry.doi]
        elif entry.pmid and entry.pmid in self._pmid_index:
            existing_id = self._pmid_index[entry.pmid]
            
        if existing_id:
            # Merge with existing entry
            existing_entry = self._entries[existing_id]
            merged_entry = self._merge_entries(existing_entry, entry)
            if merged_entry.id != existing_id:
                self._id_changes[existing_id] = merged_entry.id
            self._update_entry(merged_entry)
        else:
            self._update_entry(entry)

    def _update_entry(self, entry: Citation) -> None:
        """Update internal indexes for an entry."""
        self._entries[entry.id] = entry
        if entry.pmid:
            self._pmid_index[entry.pmid] = entry.id
        if entry.doi:
            self._doi_index[entry.doi] = entry.id

    def _merge_entries(self, entry1: Citation, entry2: Citation) -> Citation:
        """Merge two citations, preferring data from entry1."""
        merged = copy.deepcopy(entry1)
        for field in entry2.__dict__:
            if getattr(merged, field) is None and getattr(entry2, field) is not None:
                setattr(merged, field, getattr(entry2, field))
        return merged

    def get_by_id(self, id: str) -> Optional[Citation]:
        """Get citation by ID."""
        return self._entries.get(id)

    def get_by_pmid(self, pmid: str) -> Optional[Citation]:
        """Get citation by PMID."""
        if pmid in self._pmid_index:
            return self._entries[self._pmid_index[pmid]]
        return None

    def get_by_doi(self, doi: str) -> Optional[Citation]:
        """Get citation by DOI."""
        if doi in self._doi_index:
            return self._entries[self._doi_index[doi]]
        return None

    def get_id_changes(self) -> Dict[str, str]:
        """Get mapping of changed IDs."""
        return self._id_changes.copy()

    def to_bibtex(self) -> str:
        """Convert database to BibTeX format."""
        output = []
        for entry in self._entries.values():
            fields = []
            for field, value in entry.__dict__.items():
                if field != 'id' and value is not None:
                    fields.append(f"  {field} = {{{value}}}")
            
            entry_str = f"@{entry.entry_type}{{{entry.id},\n"
            entry_str += ",\n".join(fields)
            entry_str += "\n}"
            output.append(entry_str)
        
        return "\n\n".join(output)

    @classmethod
    def from_bibtex(cls, bibtex_str: str) -> 'CitationDatabase':
        """Create database from BibTeX string."""
        # Implementation would parse BibTeX and create Citation objects
        # Using existing bibtexparser or custom parser
        pass

# Example usage:
def example_usage():
    db = CitationDatabase()
    
    # Add a citation
    citation = Citation(
        id="smith2020advances",
        title="Advances in Research",
        author="Smith, John",
        year="2020",
        journal="Science Journal",
        doi="10.1234/science.1234",
        pmid="12345678"
    )
    db.add_entry(citation)
    
    # Retrieve by different identifiers
    by_id = db.get_by_id("smith2020advances")
    by_pmid = db.get_by_pmid("12345678")
    by_doi = db.get_by_doi("10.1234/science.1234")
    
    print (by_id)
    print (by_pmid)
    print (by_doi)

    # Get citation in BibTeX format
    bibtex_output = db.to_bibtex()

example_usage()



