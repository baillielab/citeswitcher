import os
import re

'''
Get the comment section from the top of each script, and insert it into the README.md file below the "EVERYTHING BELOW IS AUTOMATICALLY OVERWRITTEN BY make_readme.py" tag
'''

tag = '<!--EVERYTHING BELOW IS AUTOMATICALLY OVERWRITTEN BY make_readme.py-->'
readme = "README_tools.md"
thisdir ="./"
pycomment = '\'\'\'[\s\S]+?\'\'\''

out = {}
for script in [x for x in os.listdir(thisdir) if x.endswith(".py")]:
	with open(script) as f:
		text = f.read()
	comments = re.findall(pycomment, text)
	if len(comments)>0:
		out[script]=comments[0]

with open(readme) as f:
	rt = [line.strip() for line in f.readlines()]
if tag in rt:
	rt = rt[:rt.index(tag)]
rt+=["\n"+tag+"\n# Scriptlist\n\n"]
for s in sorted(out.keys()):
	rt+=["## {}\n\n{}\n".format(s, out[s][3:-3].strip())]
with open(readme,"w") as o:
	o.write("\n".join(rt))


