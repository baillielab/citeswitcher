#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

'''
find all capitalised abbreviations in a text file
'''

#-------------------
import string,os,sys
import io
import re
#-------------------
import argparse
parser = argparse.ArgumentParser()
# - essential
parser.add_argument('-f', '--filepath',    help='filepath')
args = parser.parse_args()
#-------------------
print(args.filepath)
#-------------------
# read input file
with io.open(args.filepath, "r", encoding="utf-8") as my_file:
     text = my_file.read()
#-------------------
abbr = re.findall("\w*[A-Z]\w*[A-Z]\w*",text)
abbr = sorted(list(set(abbr)))
for a in abbr:
    print (a)
#-------------------
