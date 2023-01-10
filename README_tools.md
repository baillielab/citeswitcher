# Workflow summary

These scripts do a number of useful tasks in academic writing. Each script does a single, relatively simple, task. This means they are easier to use and understand than a single catch-all automated wrapper, but it does mean that you have to use a series of them to make your manuscript. To get around this, I find it useful to keep a note of the build commands for each manuscript in a comment at the top of the file.






<!--EVERYTHING BELOW IS AUTOMATICALLY OVERWRITTEN BY make_readme.py-->
# Scriptlist


## abbreviation_finder.py

find all capitalised abbreviations in a text file

## format_author_list.py

split a comma separated string and collapse into a single list

## get_figure_numbers.py

read a tex or md file
write a json file for replacing contents

## make_readme.py

Get the comment section from the top of each script, and insert it into the README.md file below the "EVERYTHING BELOW IS AUTOMATICALLY OVERWRITTEN BY make_readme.py" tag

## remove_figures.py

remove figures from a markdown file
but leave the captions in place

## tex2md.py

convert a tex file to markdown
depends on all tables using latex csv input \\csvautotabular \\csvautotabularcenter
manually converts all TABLES and TABLE CROSSREFERENCES
then uses pandoc for the rest
