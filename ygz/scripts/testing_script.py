import aipy as a, numpy as n, pylab as p, ephem as e

#aa=a.cal.get_aa('psa898_v003',n.array([.15]))
aa=a.cal.get_aa('paper128',n.array([.15]))
nants=32
rad2deg=180/n.pi
src= a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")

p.figure()

for time in n.arange(2456240,2456241,.1):
        aa.set_jultime(time)
       #     aa.set_jultime(2456240.5)
            #print aa.get_jultime()
        src.compute(aa)
        print time, src.get_crds('top',ncrd=2)
 
