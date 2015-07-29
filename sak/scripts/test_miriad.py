import aipy as a
import sys
VERB=True
for f in sys.argv[1:]:
	print f,
	sys.stdout.flush()
	try:
		uv = a.miriad.UV(f)
		for all in uv.all():
			continue
		print "[OK]"
		sys.stdout.flush()
	except:
		print "[ERROR]"
		if VERB:
			print sys.exc_info()
			sys.stdout.flush()
