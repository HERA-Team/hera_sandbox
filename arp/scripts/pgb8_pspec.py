#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, sys

x = n.array([1e2,2e2,4e2,1e3,3e3,1e4])
y = n.array([1.,3,10,30,60,50])
C_l_sim = n.polyfit(n.log10(x), n.log10(y), deg=3)

def C_l_sync(L, freq):
    return L*(L+1)/(2*n.pi) * 700. * (1e3/L)**2.4 * (130./freq)**(2*2.8)

f = open('pgb8_pspec_dat.txt')
#f = open(sys.argv[-1])

#ells = n.array(map(float, f.readline().split()))
ells = n.array(map(float, f.readline().split())) / 6.25
print ells

chunks = f.read().split('\n\n')
chs, freqs, inds, ind_errs = [], [], [], []
for cnt, chunk in enumerate(chunks):
    c1, c2, c3 = chunk.split('-----------------------------------')
    words = c1.split()
    ch = float(words[2][3:-1])
    freq = float(words[3][5:])
    c_l1 = n.array(map(float, c1[c1.find('[')+1:c1.find(']')].split()))
    c_l2 = n.array(map(float, c2[c2.find('[')+1:c2.find(']')].split()))
    errs = n.array(map(float, c3[c3.find('[')+1:c3.find(']')].split()))

    if ch not in [87.5, 117.5, 147.5, 177.5]: continue
    print freq

    def pwrlaw(x):
        A, ind = x
        wgt = 1 / errs
        wgt /= wgt.sum()
        #wgt = 1
        dif = (c_l1 - A*ells**ind) * wgt
        return n.sqrt((dif**2).sum())
    A,ind = a.optimize.fmin(pwrlaw, n.array([0., 0.]))
    c_l3 = A*ells**ind
    sig = pwrlaw([A, ind])
    wgt = 1 / errs
    ind_err = sig / (A*ind*ells**(ind-1) * wgt).sum()
    print sig, ind, ind_err
    chs.append(ch); freqs.append(freq)
    inds.append(ind)
    ind_errs.append(ind_err)

    # Also calculate C_l = (L*(L+1)/2pi) c_l
    C_l1 = ells*(ells+1)/(2*n.pi) * c_l1 * 1e6
    C_l2 = ells*(ells+1)/(2*n.pi) * c_l2 * 1e6
    C_l3 = ells*(ells+1)/(2*n.pi) * c_l3 * 1e6
    Errs = ells*(ells+1)/(2*n.pi) * errs * 1e6

    color = .8 * float(cnt) / len(chunks)
    color = (color, color, color)
    Errs_dn = n.where(Errs >= C_l1-1e-6, C_l1-1e-6, Errs)
    p.errorbar(ells, C_l1, [Errs_dn, Errs],
        fmt='.-', label='ch%fcl' % ch, color=color)
    p.plot(ells, C_l2, '-.', label='ch%ftherm' % ch, color=color, linewidth=2)
    #p.plot(ells, C_l3, ':', label='ch%ffita' % ch, color=color)
    if False:
        ax = p.gca(); ax.set_xscale('log'); ax.set_yscale('log')
        p.xlim(1e2,1e4)
        p.ylim(1e-5,1e1)
        p.xlabel(r'$\ell=2\pi u$', fontsize=20)
        p.ylabel(r'$\frac{\ell(\ell+1)}{2\pi}C_\ell~\left[mK^2\right]$', 
            fontsize=20)
        p.show()

p.plot(ells, 10**n.polyval(C_l_sim, n.log10(ells)), 'k-', 
    label='z=9.2_santos2005', linewidth=4)
p.plot(ells, C_l_sync(ells,147.), 'k:', 
    label='sync')
ax = p.gca(); ax.set_xscale('log'); ax.set_yscale('log')
#p.xlim(1e2,1e4)
p.xlim(60,1e3)
#p.ylim(1e0,1e12)
p.ylim(1e-1,1e11)
p.xlabel(r'$\ell=2\pi u$', fontsize=20)
p.ylabel(r'$\ell(\ell+1)C_\ell/2\pi\ \ (mK^2)$', fontsize=20)
p.show()

p.errorbar(freqs, inds, 2*n.array(ind_errs), fmt='k.')
p.xlabel(r'Frequency (MHz)')
p.ylabel(r'$C_\ell$ Power-Law Index')
p.show()
