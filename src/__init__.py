import gb_weather, pfb, pspec, dspec, red, fringe, frf_conv, miriad
import oqe
import warnings
try: import omni, uCal
except(ImportError,NameError):
    warnings.warn("Warning: omnical not installed, not importing capo.omni")
import arp, jcp, dfm, dcj, zsa, ctc#, wyl
