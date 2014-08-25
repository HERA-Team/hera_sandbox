#! /usr/bin/env python
import sys, os
import smtplib
import time
import datetime
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email import Encoders

def cataddr(addrs):
    address_string = ''
    for a in addrs:
        address_string+=a+','
    address_string.strip(',')
    return address_string

##----Set stuff----##
emails = ['ddeboer@berkeley.edu','william@ska.ac.za','daniel.c.jacobs@asu.edu','davidm@astro.berkeley.edu']
tooLongPing_minutes = 12.0
tooShortAlert_minutes = 60.0
watchFile = 'watch.dat'
alertFile = 'alerts.dat'

##----Read watchdog file----##
watchFileExists = True
try:
    fp = open(watchFile,'r')
    times = []
    ping = []
    elineno = 0
    inEntries = False
    for line in fp:
        line = line.strip()
        if line == 'Watching':
            elineno = 0
            inEntries = True
        elif elineno == 1:
            etime = time.strptime(line,'%a %b %d %H:%M:%S UTC %Y')
            etime = datetime.datetime.fromtimestamp(time.mktime(etime))
            times.append(etime)
        elif elineno == 3:
            ds = line.split()
            try:
                pingTime = float([x for x in ds if x[0:2]=='ti'][0].split('=')[1])
            except IndexError:
                pingTime = 0.0
            ping.append(pingTime)
        if inEntries:
            elineno+=1
    fp.close()
    ##----Check time of last entry----##
    unow = datetime.datetime.utcnow()
    diffLast = unow - times[-1]
    diffLast = diffLast.total_seconds()/60.0
except IOError:
    watchFileExists = False
    diffLast = tooLongPing_minutes*10.0

##----Take action if needed----##
if diffLast > tooLongPing_minutes:
    # Check when last alert sent
    sendAlert = False
    try:
        fp = open(alertFile,'r')
        alerts = []
        for line in fp:
            atime = time.strptime(line,'%Y-%m-%d %H:%M:%S.')
            atime = datetime.datetime.fromtimestamp(time.mktime(atime))
            alerts.append(atime)
        fp.close()
        diffAlerts = unow - alerts[-1]
        diffAlerts = diffAlerts.total_seconds()/60.0
        if diffAlerts > tooShortAlerts_minutes:
            sendAlert = True
    except IOError:
        sendAlert = True

    if sendAlert:
        # Log alert to file
        fp = open(alertFile,'a')
        fp.write(str(unow))
        fp.close()
        if watchFileExists:
            emailmsg = 'Too long between obs --> morph pings (%.1f min)\n\n' % (diffLast)
        else:
            emailmsg = 'Watchdog file %s does not exist' % (watchFile)
        from_addr = 'david.r.deboer@gmail.com'
        server = smtplib.SMTP('smtp.gmail.com',587)
        server.ehlo()
        server.starttls()
        server.ehlo()
        server.login('david.r.deboer@gmail.com','2001SpaceOdD')

        msg = MIMEMultipart()
        msg['Subject'] = 'Watchdog bites'
        msg['From'] = 'david.r.deboer@gmail.com'
        msg['To'] = cataddr(emails)
        msg.attach(MIMEText(emailmsg))
        part = MIMEBase('application',"octet-string")
        part.set_payload(open('watch.dat',"rb").read())
        Encoders.encode_base64(part)
        part.add_header('Content-disposition','attachment; filename="%s"' % os.path.basename('watch.dat'))
        msg.attach(part)
        server.sendmail(from_addr,emails,msg.as_string())
        server.quit()
