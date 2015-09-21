import os
import glob

"""
This is the STANDARD prep procedure. To make <...> I'll need to change things up a little.
<...> = |V_Q|^2/|V_I|^2
<...> = P = sqrt(|V_Q|^2 + |V_U|^2)
"""

for file in glob.glob('*.uv'):
	os.system('mk_stokes.py '+file+' --stokes=I,Q,U,V') #uv->uvP
	
	for pol in ['I','Q','U','V']:
		os.system('python pull_pol_SK.py '+file+'P -p '+pol) #uvP->uvP{I,Q,U,V}
		os.system("python pspec_prep_jcp.py --horizon=50 --nophs --nogain --window='blackman-harris' --model -p "+pol+" -C psa819_v008 "+file+"P"+pol) #uvP{I,Q,U,V}->uvP{I,Q,U,V}{B,F}
		print 'NOTICE ME'
		print "NOTICE ME"
		print "uv_addsub.py "+file+"P"+pol+"F "+file+"P"+pol+"BR"
		os.system("xrfi_simple.py -n 3 "+file+"P"+pol+"B") #uvP{I,Q,U,V}B -> uvP{I,Q,U,V}BR
		
		os.system("uv_addsub.py "+file+"P"+pol+"F "+file+"P"+pol+"BR") #-> uvP{I,Q,U,V}Fa
	
	#os.system('python mapI2xx.py ')
