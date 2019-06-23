#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

'''
temporary installation - goes away after restart
'''


import os
import subprocess

scriptpath = os.path.dirname(os.path.realpath(__file__))
mainscript = "fixcitations.py"

cmd = "chmod +x {}".format(os.path.join(scriptpath,mainscript))
print (cmd)
subprocess.call(cmd, shell=True)

cmd = 'export PATH=$PATH":{}"'.format(scriptpath)
print (cmd)
subprocess.call(cmd, shell=True)

