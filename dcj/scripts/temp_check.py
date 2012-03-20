#!/usr/bin/env python

#
#if the log file is staler than 4 hours, 
#Check the input temp log for over limit internal temp and send Danny an email if it is.
#
#Usage check_temps.py <tempfile> <limit in [K]>
#
#Assumes we care about the labjack internal temp at column 2 (first actual temp)
#
#
#

import numpy as n
import sys,os,statvfs,aipy as a
from glob import glob
import smtplib
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from time import time,asctime


logfile = '/home/obs/logs/overtempwarnings.txt'
age = 4. #check interval
if os.path.exists(logfile):
    if (time() - os.path.getmtime(logfile))/3600<age:
        print "warning log too young, not checking again at this time"
        sys.exit()

cutoff = float(sys.argv[2])
file = sys.argv[1]

lines = open(file).readlines()
temps = n.array([line.split()[2] for line in lines]).astype(float)
overtemps = n.argwhere(temps>cutoff)

if len(overtemps)==0:
    print "peak temp of %7.3fK is below %7.3fK cutoff"%(temps.max(),cutoff)
    sys.exit(0)
else:
    warning = "-----%s-----\n%i records out of %d are above %7.3fK\nMost recent: %7.3fK at %s\n"%(
    asctime(),
    len(overtemps),
    len(lines),
    cutoff,
    temps[overtemps[-1]],
    lines[overtemps[-1]])

open(logfile,'a').write(warning)


msg = MIMEMultipart()


To = "wheresmytab@gmail.com"
From = "teampaper@gmail.com"



msg['Subject'] = "WARNING: PSA Overtemp"
msg['From'] = From
msg['To'] = To
msg.attach(MIMEText(warning))

os.system('python ~/bin/plot_all_temps.py `ls -tr /home/obs/Temperature_Logs/2*txt | tail -n 48`')
msg.attach(MIMEImage(open('/home/obs/Temperature_Logs/all_temps.png').read()))
print "sending email"
username = 'teampaper'
password = 'b00lardy'

# Send email
server = smtplib.SMTP('smtp.gmail.com:587')
server.starttls()
server.login(username,password)
server.sendmail(From, [To], msg.as_string())
server.quit()

print "email sent"
