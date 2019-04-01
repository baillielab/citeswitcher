#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: utf-8

'''

move all square brackets before fullstops and commas to after


'''

#-------------------
import io
import re
#-------------------
import argparse
parser = argparse.ArgumentParser()
# - essential
parser.add_argument('-f', '--filepath',    help='filepath')
args = parser.parse_args()
#-------------------
squarebrackets = '\[[\s\S]+?\]\.'
punctuationmarks = ['.',',']
#-------------------

def move_citations(thistext):
	for p in punctuationmarks:
		allfound = re.findall(squarebrackets,thistext)
	for x in allfound:
		print ("-->", x)



#-------------------
if __name__ == "__main__":
	with io.open(args.filepath, "r", encoding="utf-8") as f:
		text = f.read()
	text = move_citations(text)
	'''
	with io.open(args.filepath, "w", encoding="utf-8") as o:
		o.write(text)
	'''
