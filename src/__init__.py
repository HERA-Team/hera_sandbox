import pfb, pspec, dspec, red, fringe, miriad, linsolve, xrfi
import fringe as frf_conv # for backward compatibility
import oqe, hex, metrics
import warnings
try: import omni, uCal#, wyl
except(ImportError,NameError):
    warnings.warn("Warning: omnical not installed, not importing capo.omni")
import arp, jcp, dfm, dcj, zsa, ctc, sak#, wyl
