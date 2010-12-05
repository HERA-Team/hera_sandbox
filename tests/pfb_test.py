import capo.pfb as P

def plot_pfb_windows(L=4096, taps=4, fwidth=1):
    import pylab
    for w in P.WINDOWS:
        P.__set_pm__(L, w, taps, fwidth)
        pylab.plot(P.pm['window'], label=w)
    pylab.legend()
    pylab.show()

if __name__ == '__main__':
    plot_pfb_windows()        
