#! /usr/bin/env python
import sys, os
import smtplib
import time
import datetime
import WatchAlert
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
emails = []
emailmsg = '' 
try:
    fp = open('watchers.txt','r')
    for line in fp:
        if line[0] == '#':
            continue
        line = line.strip()
        emails.append(line)
    fp.close()
except:
    emailmsg+= 'No e-mails read\n\n'
    emails = ['ddeboer@berkeley.edu','david.r.deboer@gmail.com']
tooLongPing_minutes = 35.0
tooShortAlert_minutes = 60.0
watchFile = 'watch.txt'
alertFile = 'alerts.txt'

##----Read watchdog file----##
w = WatchAlert.Watch(watchFile)
w.readFile()
diffLast = w.elapsed()

##----Take action if needed----##
skip = True
if not skip:
#if diffLast > tooLongPing_minutes:
    unow = datetime.datetime.utcnow()
    # Check when last alert sent
    sendAlert = False
    a = WatchAlert.Alert(alertFile)
    a.readFile()
    diffAlerts = a.elapsed()
    if diffAlerts > tooShortAlert_minutes:
        sendAlert = True

    if sendAlert:
        # Log alert to file
        a.addEntry(unow)
        if w.watchFileExists:
            emailmsg+= 'Too long between obs --> morph pings (%.1f min)\n\n' % (diffLast)
        else:
            emailmsg+= 'Watchdog file %s does not exist' % (watchFile)
        from_addr = 'david.r.deboer@gmail.com'
        server = smtplib.SMTP('smtp.gmail.com',587)
        server.ehlo()
        server.starttls()
        server.ehlo()
        server.login('david.r.deboer@gmail.com','')

        msg = MIMEMultipart()
        msg['Subject'] = 'Watchdog bites'
        msg['From'] = 'david.r.deboer@gmail.com'
        msg['To'] = cataddr(emails)
        msg.attach(MIMEText(emailmsg))
        part = MIMEBase('application',"octet-string")
        part.set_payload(open(watchFile,"rb").read())
        Encoders.encode_base64(part)
        part.add_header('Content-disposition','attachment; filename="%s"' % os.path.basename(watchFile))
        msg.attach(part)
        server.sendmail(msg['From'],emails,msg.as_string())
        server.quit()
else:
    ##----Get bandwidth info----##
    b = WatchAlert.Bandwidth('watchBand.dat',bandTransmitSize=150000000)
    b.getbw()
    b.addEntry()

