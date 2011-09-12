"""
A set of functions for handling weather data.
DCJ 23 May 2009
"""

def get_gb_weather_channel(t,thing):
   import glob, numpy as n, pfits, aipy as a,ephem
   weather_dir = '/data1/paper/gb/weather/'
   tstring = str(ephem.Date(a.phs.juldate2ephem(t)).datetime().date()).replace('-','_')
   tstring = tstring +'_'+str(ephem.Date(a.phs.juldate2ephem(t)).datetime().time())[:3]   
   tfile = glob.glob(weather_dir+tstring+'*')[0]
   fitsfile = pfits.FITS(tfile)
   hdus = fitsfile.get_hdus()
   ttimes = hdus[1].get_data()['DMJD']+2400000.5
   things_data = hdus[1].get_data()[thing]
   thing = things_data[n.argwhere(n.abs(ttimes-t) == n.min(n.abs(ttimes-t)))[0]]
   return thing
def get_gb_temp(t):
   import glob, numpy as n, pfits, aipy as a,ephem
   weather_dir = '/data1/paper/gb/weather/'
   tstring = str(ephem.Date(a.phs.juldate2ephem(t)).datetime().date()).replace('-','_')
   tstring = tstring +'_'+str(ephem.Date(a.phs.juldate2ephem(t)).datetime().time())[:3]
   tfile = glob.glob(weather_dir+tstring+'*')[0]
   fitsfile = pfits.FITS(tfile)
   hdus = fitsfile.get_hdus()
   ttimes = hdus[1].get_data()['DMJD']+2400000.5
   temp_data = hdus[1].get_data()['TEMPERATURE']
   temp = temp_data[n.argwhere(n.abs(ttimes-t) == n.min(n.abs(ttimes-t)))[0]]
   return temp
