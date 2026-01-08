# YAML File Search and File Creation Logic

This document describes how the `fixcitations.py` script searches for YAML files, discovers bibliography files, creates missing files, and handles CSL (Citation Style Language) files.

## Overview

The script processes citation references in markdown/LaTeX files and requires bibliography files (`.bib`) to work. It searches for bibliography specifications in YAML files (both inline in document headers and sidecar files) and creates missing files as needed.

## YAML File Search Logic

### Search Order

The script searches for YAML files in the following order:

1. **Inline YAML Header** (from input file)
   - Extracted from the document header using `readheader()` function
   - Looks for YAML between `---` markers or `---`/`...` markers
   - Only the first YAML document is parsed (handles multi-document YAML safely)

2. **Sidecar YAML Files** (in order of priority):
   - `{filestem}.yaml` in source directory
   - `{filestem}.yml` in source directory
   - `_quarto.yml` in source directory
   - `_quarto.yml` in parent directory (`../`)
   - `_quarto.yml` in grandparent directory (`../../`)

3. **Command-line specified YAML** (`-y`/`--yaml` flag):
   - If absolute path: uses as-is
   - If relative path: tries both absolute resolution and relative to sourcepath
   - **ISSUE**: Both resolved absolute path AND relative path are added to candidates, potentially causing duplicate processing

### YAML Processing

- Only the `bibliography` key is extracted from YAML files (other keys are ignored)
- Bibliography paths are resolved relative to:
  - **Inline YAML**: relative to `sourcepath` (input file's directory)
  - **Sidecar YAML**: relative to the YAML file's directory
- All bibliography paths are normalized to absolute paths for consistency
- The script tracks which YAML file (or inline) provided each bibliography entry via `bib_to_yaml_file` mapping

### YAML Parsing Function

`load_first_yaml_doc()`:
- Strips YAML front-matter fences (`---` or `...`)
- Uses `yaml.safe_load_all()` to handle multi-document YAML
- Returns only the first document that is a dictionary
- Returns `None` on any error (silent failure)

## Bibliography File Discovery and Creation

### Discovery Priority

Bibliography files are discovered in this priority order:

1. **Command-line override** (`-b`/`--bibfile` flag)
   - Highest priority
   - If relative path, checks both absolute resolution and relative to sourcepath

2. **YAML 'bibliography' entries**
   - In the order they appeared in YAML files
   - **ISSUE**: All YAML-specified bibliography files are created even if command-line override exists (wasteful)

3. **Default `cs.bib`** in source folder
   - Only if file exists

4. **Other `.bib` files** in source folder
   - Alphabetically sorted
   - Only existing files are considered

5. **Fallback**: `_bib/cs.bib` in source folder
   - Created if no other bibliography found

### File Creation Logic

#### YAML-Specified Bibliography Files

```python
if yaml_bibs:
    for b in yaml_bibs:
        bib_path = Path(b)
        if not bib_path.exists():
            bib_path.parent.mkdir(parents=True, exist_ok=True)
            bib_path.write_text("", encoding="utf-8")
```

- Creates empty `.bib` files for any YAML-specified bibliography that doesn't exist
- Creates parent directories as needed
- **ISSUE**: This happens before bibliography selection, so files may be created unnecessarily

#### Fallback Bibliography File

```python
if chosen_bib is None:
    chosen_bib_path = Path(sourcepath) / "_bib" / "cs.bib"
    chosen_bib_path.parent.mkdir(parents=True, exist_ok=True)
    if not chosen_bib_path.exists():
        chosen_bib_path.write_text("", encoding="utf-8")
```

- Creates `_bib/cs.bib` if no bibliography file is found
- Only creates if it doesn't already exist

### Deduplication

- Bibliography candidates are deduplicated using resolved absolute paths
- Uses string comparison of `Path.resolve()` results
- **ISSUE**: Path resolution may not catch all duplicates if symlinks are involved differently

## CSL (Citation Style Language) File Handling

### CSL Discovery

- Checks all YAML files (inline and sidecar) for `csl` key
- Stops at first found `csl` entry

### CSL Creation

**Condition**: CSL is created only if:
1. No `csl` key exists in any YAML file
2. The chosen bibliography came from YAML (not command-line or default)

**Process**:
1. Copies `csl/minimal.csl` from script directory to same directory as chosen bibliography
2. Calculates relative path from YAML file's directory to CSL file
3. Updates the YAML file that specified the bibliography:
   - **Inline YAML**: Updates `original_header` variable (used when writing output)
   - **Sidecar YAML**: Reads, updates, and writes back to file

### CSL Path Calculation

```python
# For inline YAML
csl_rel_path = str(target_csl.relative_to(Path(sourcepath)))

# For sidecar YAML  
csl_rel_path = str(target_csl.relative_to(yaml_dir))
```

- Uses `Path.relative_to()` with fallback to absolute path if paths aren't related
- **ISSUE**: `Path(sourcepath)` is redundant - `sourcepath` is already a Path object

### YAML File Writing

When updating sidecar YAML files:
- Attempts to preserve original format (fences, closing markers)
- Logic is fragile:
  ```python
  has_closing_fence = yaml_content.strip().endswith('---') or yaml_content.strip().endswith('...')
  if yaml_content.strip().startswith('---'):
      yaml_file_path.write_text("---\n" + yaml_str + ("---\n" if has_closing_fence else ""), encoding="utf-8")
  ```
- **ISSUE**: Doesn't preserve exact original formatting (whitespace, comments, etc.)

## Output File Creation

### Bibliography Output

```python
localbibpath_obj = Path(localbibpath)
bibstem = localbibpath_obj.stem
localbibpath = str(localbibpath_obj.parent / (bibstem + citelabel + "bib"))
localbibpath_path.parent.mkdir(parents=True, exist_ok=True)
localbibpath_path.write_text(outbib, encoding="utf-8")
```

- Creates output bibliography with suffix based on `citelabel` (e.g., `.citemd.bib`)
- Creates parent directories as needed

### Global Bibliography

- Only writes if content has changed (hash comparison)
- Creates parent directories as needed
- **ISSUE**: Uses `args.safemode` but `args` is not passed to `main()` function - this will cause an error!

### Text Output File

```python
outputfile_path = Path(outpath) / (filestem + citelabel + input_file_extension)
outputfile_path.parent.mkdir(parents=True, exist_ok=True)
```

- Creates output file with appropriate suffix
- Includes original YAML header if output style is markdown

## Identified Issues and Recommendations

### Critical Issues

1. **`args.safemode` Reference Error** (Line 1610) ⚠️ **CRITICAL BUG**
   - `args` is not available in `main()` function scope
   - `main()` function receives `safemode` as a parameter (line 1281)
   - Code incorrectly uses `args.safemode` which will raise `NameError`
   - **Fix**: Change `if not args.safemode:` to `if not safemode:` on line 1610

2. **YAML Candidate Duplication**
   - When `yaml_file` is provided, both absolute and relative paths added
   - **Fix**: Add only the resolved absolute path, or check for duplicates before adding

### Efficiency Issues

3. **Unnecessary File Creation**
   - YAML-specified bibliography files created even when command-line override exists
   - **Fix**: Only create files after bibliography selection, or skip creation if higher-priority option exists

4. **Broad Exception Handling**
   - `except Exception:` hides all errors during YAML processing
   - **Fix**: Catch specific exceptions and log errors appropriately

### Code Quality Issues

5. **Path/String Mixing**
   - Inconsistent mixing of Path objects and strings
   - **Fix**: Use Path objects consistently, convert to strings only when needed for external APIs

6. **Redundant Path Conversion**
   - `Path(sourcepath)` when `sourcepath` is already a Path
   - **Fix**: Use `sourcepath` directly

7. **Fragile YAML Format Preservation**
   - YAML file writing doesn't preserve exact formatting
   - **Fix**: Use a YAML library that preserves formatting, or implement more sophisticated format preservation

8. **Deduplication Logic**
   - String-based deduplication may miss symlink-related duplicates
   - **Fix**: Use `Path.resolve()` consistently and compare Path objects directly

### Suggested Improvements

9. **YAML File Tracking**
   - Good: Tracks which YAML file provided each bibliography
   - Could be improved: Also track order of YAML files for better priority handling

10. **Error Messages**
    - Silent failures in YAML parsing
    - **Fix**: Add logging/error messages for debugging

11. **Path Resolution Consistency**
    - Some places use `resolve()`, others don't
    - **Fix**: Standardize on using `resolve()` for all absolute path comparisons

## Code Flow Summary

```
1. Read input file → Extract inline YAML header
2. Search for sidecar YAML files (in priority order)
3. Extract 'bibliography' entries from all YAML sources
4. Create any missing YAML-specified bibliography files
5. Build bibliography candidate list (command-line, YAML, defaults)
6. Select bibliography based on priority
7. Create fallback bibliography if none found
8. Check for CSL in YAML files
9. Create CSL file and update YAML if needed
10. Process citations and create output files
```

## Key Data Structures

- `yaml_bibs`: List of absolute paths to bibliography files from YAML
- `bib_to_yaml_file`: Dict mapping bibliography path → YAML file path (or None for inline)
- `yaml_files_checked`: List of tuples `(yaml_path, yaml_content)` for CSL checking
- `bib_candidates`: List of tuples `(path, reason)` for bibliography selection
- `existing_candidates`: Filtered list of candidates that actually exist

