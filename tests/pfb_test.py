import capo.pfb as P
import numpy as n

def get_response(window):
    L = window.size
    x = 2*n.pi*n.arange(L, dtype=n.float) / L
    fqs = n.arange(0,L/8-1) + 0.5
    resp = n.abs([n.fft.rfft(window*n.sin(f*x))[0]/L for f in fqs])**2
    return fqs, resp
def plot_pfb_windows(L=4096, taps=4, fwidth=1):
    import pylab
    for w in P.WINDOWS:
        print w
        P.__set_pm__(L, w, taps, fwidth)
        try:
            pylab.subplot(121)
            pylab.plot(P.pm['window'], label=w)
            pylab.subplot(122)
            fqs,resp = get_response(P.pm['window'].squeeze())
            pylab.semilogx(fqs, 10*n.log10(resp))
        except(AttributeError): pass
    pylab.subplot(121)
    pylab.grid()
    pylab.legend()
    pylab.subplot(122)
    pylab.grid()
    pylab.show()

if __name__ == '__main__':
    plot_pfb_windows()        
