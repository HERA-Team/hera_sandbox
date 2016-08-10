POLNUM = {
    'x':0, # factor to multiply ant index for internal ordering
    'y':1,
    'l':2,
    'r':3,
    'a':4,
    'b':5,
}

NUMPOL = dict(zip(POLNUM.values(), POLNUM.keys()))


class Antpol(int):
	def __init__(self, ant, pol, nant):
		self.ant = ant
		self.pol = pol
		self.nant = nant
	
	def __new__(cls, ant, pol, nant):
		obj = int.__new__(cls, POLNUM[pol] * nant + ant)
		obj.pol = pol
		obj.nant = nant
		return obj
	
	#which way around do I want str and repr? str should be the physical number and pol?
	def __str__(self): return str(self.ant)+self.pol
	def __repr__(self): return str(POLNUM[self.pol] * self.nant + self.ant)#+self.pol
	

ap1x = Antpol(1,'x',3)
ap1y = Antpol(1,'y',3)
ap2x = Antpol(2,'x',3)

print ap1x
print ap2x
print ap1x.pol
print ap1y.pol
print str(ap1x)
print str(ap1y)

bl1 = (ap1x,ap2x)
bl2 = (ap1x,ap1y)

list = [bl1,bl2]
