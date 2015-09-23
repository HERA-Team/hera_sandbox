"""
rename_ext.py

Rename file extensions in bulk for a specified folder
Usage: python rename_ext.py target_folder old_ext new_ext

"""

import os,sys,glob

args = sys.argv[1:]
assert(len(args)==3) #Usage: python rename_ext.py target_folder old_ext new_ext
assert(args[1][0]=='.') #old extension must start with "."
assert(args[2][0]=='.') #new extension must start with "."

for f in glob.glob('%s/*%s'%(args[0],args[1])):
	FF = f.split('.')
	f_new=''
	for i,part in enumerate(FF):
		if part!=args[1].split('.')[1]: f_new+=part+'.'
	f_new+=args[2][1:]
	print 'mv %s %s'%(f,f_new)
	os.system('mv %s %s'%(f,f_new))
	
