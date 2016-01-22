import os, glob

TEST = True
EXTENSION = 'uvcRREc'
PATH2FOLDERS='../psa_live/psa*'
MOVE = True #as opposed to copy

print 'TEST=',TEST
print 'GLOB=',EXTENSION
print 'PATH2TARGET=',PATH2FOLDERS
print 'MOVE instead of COPY?',MOVE

for psafolder in glob.glob(PATH2FOLDERS):
	psaxxxx = psafolder[-7:]
	if os.path.exists(psaxxxx): 
		print psaxxxx,'already here'
		continue
	else: 
		#print psaxxxx, 'not in pol_live'
		content = []
		for f in glob.glob(psafolder+'/*.%s'%EXTENSION): 
			content.append(f) #get all the files from that night (in no particular order)
		content = sorted(content)
		
		try: print '1st file',content[0][-25:]
		except IndexError:
			print 'No %s files'%EXTENSION
			continue
		
		#if we get to this point, things are looking OK to create a folder and
		#put the observations in there
		if MOVE: cmd_str = 'mv'
		else: cmd_str = 'cp'
		
		if TEST: 
			print 'mkdir '+psaxxxx
			for f in content:
				print cmd_str+' -r '+f+' '+psaxxxx
		else:
			print 'mkdir '+psaxxxx
			os.system('mkdir '+psaxxxx)
			for f in content:
				print cmd_str+' -r '+f+' '+psaxxxx
				os.system(cmd_str+' -r '+f+' '+psaxxxx)
		
			
