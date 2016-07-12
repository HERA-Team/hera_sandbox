#! /usr/bin/env python

"""
Pipe in traceroute output, return a google maps url

90% of the code borrowed from Kyle Flanagan 
http://kyleflan.wordpress.com/2012/12/13/a-visual-traceroute-program/

"""
from pyipinfodb import *
import urllib, urllib2, sys, os


# api key for ipinfodb.com
API_KEY = "687f07cd5c6a20f0d7a6890751f049a3745c95f98c7706500bfd1fce73f0a1d0"
# URL for Google Static Map API
google_map_base_url = "http://maps.googleapis.com/maps/api/staticmap?"

def getURLDict(ips):
    """
    Returns a dictionary consisting of URL paramenters to pass to the Google
    static map api. The IP addresses in ips are looked up using the IPInfo
    class from pyipinfodb and constructed into the url_dict.
    """
    url_dict = {'zoom' : '', 
            'maptype':'roadmap', 
            #'sensor':'false', 
            'size':'1200x600', 
            'markers':[], 
            'path':'color:0x0000ff|weight:3',
            'center':''}
    for x, ip in enumerate(ips):
        ipi = IPInfo(API_KEY, ip, city=True) # get city, state
        location = ipi.getCity() + ',' + ipi.getRegion()
        print "IP: " + ip + " Location: " + location 

        # the IPInfo class returns '-' if it cannot find a city or
        # state, so ignore these values
        if '-' in location:
            continue

        # we could use the heuristic below to find an approximate center
        # for the map, but since we're passing markers and a path, Google
        # will find the center for us
        ##if len(ips) / 2 <= x and url_dict['center'] == None:
        ##    url_dict['center'] = ipi.getCity() + ',' + ipi.getRegion() 

        # append markers
        if x == len(ips)-1: # end with red marker
            url_dict['markers'].append('label:' + str(x) + '|' + location)
        else: # else use a green one
            url_dict['markers'].append('label:' + str(x) + '|' + location)

        # append to the path route
        url_dict['path'] = url_dict['path'] + "|" + location
    return url_dict

lines = sys.stdin.readlines()
if len(lines)<2:
    print "input error."
    sys.exit(1)


lines = [l.strip() for l in lines]
print "found %d hops"%len(lines)


print "Looking up IP addresses. Please wait..."
# first line is traceroute info output, remove it
lines.pop(0)
# now get hostname, ip from each line
ips = []
for line in lines:
    if line != "":
        # if we didn't get a reply from a gateway, ignore that line in the
        # traceroute output
        if '*' in line:
            continue
        
        # Split the line and extract the hostname and IP address
        split_line = line.split('  ')
        ips.append(split_line[1].split(' ')[1][1:-1])

print "IP addresses to be looked up:"
for ip in ips:
    print ip

url_dict = getURLDict(ips)

urldata = urllib.urlencode(url_dict, True)
print "Google Map URL (copy and paste into your web browser):"
url = google_map_base_url + urldata
print url
print "length of api call: ",len(urldata)


