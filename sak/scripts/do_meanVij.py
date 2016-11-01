import sys, os
for uv in sys.argv[1:]:
    os.system('~/ReposForCanopy/capo/sak/scripts/meanVij.py -C psa6622_v003 -p xx %s'%uv)
