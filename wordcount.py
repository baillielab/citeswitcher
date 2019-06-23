#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

import re
import include
import citefunctions
#-----------------------------
def wordcount(text):
    return len(re.findall(r'\w+', text))
#-----------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', default=None,   help='filename')
    parser.add_argument('-i', '--include', action="store_false", default=True, help='do NOT include files')
    parser.add_argument('-c', '--count_comments', action="store_true", default=False, help='also count the words in comments')
    args = parser.parse_args()
    if args.include:
        text = include.parse_includes(args.filename)
    else:
        with io.open(args.filename, "r", encoding="utf-8") as f:
            text = f.read()
    header, text = citefunctions.readheader(text)
    if not args.count_comments:
        text = include.stripcomments(text)
    print (wordcount(text))


















