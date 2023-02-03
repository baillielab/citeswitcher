import os
import cairosvg
#-----------------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filepath',    help='single file')
parser.add_argument('-d', '--dirpath',  default=None,   help='whole directory')
parser.add_argument('-o', '--overwrite',    action="store_true", default=False,    help='over-write existing files')
args = parser.parse_args()
#-----------------------------
if args.dirpath:
	thisdir = args.dirpath
	files = [x for x in os.listdir(thisdir) if x.endswith('.svg')]
else:
	thisdir = os.path.split(args.filepath)[0]
	files = [args.filepath]
for x in files:
	filepath = os.path.join(thisdir, x)
	outfilepath = filepath.replace(".svg",".pdf")
	if os.path.exists(outfilepath) and not args.overwrite:
		print ("File already exists: {}. Use '-o' argument to overwrite output pdf files".format(outfilepath))
		continue
	cairosvg.svg2pdf(url=filepath, write_to=outfilepath)

'''
alternative: use inkscape
'''