#!/usr/bin/env python

#
#  Checks the correlator log to see if it is staler than 1 hour or the errorlog is fresher than 1 hour
#  If it is email the start_correlator log to Danny
#  usage: cn_rx_log_check.py  <log> <errorlog>
#

import numpy as n
import sys,os,statvfs,aipy as a
from glob import glob
import smtplib
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from time import time,asctime


#logfile = '/home/obs/logs/start_correlator.log'
logfile = sys.argv[1]
errlog = sys.argv[2]
agelim = 1. #check interval
if os.path.exists(logfile):
    age = (time() - os.path.getmtime(logfile))/3600
#    if age<agelim:
#        print logfile,' age =',age
#        sys.exit()
else:
    age = -1

if os.path.exists(errlog):
    errage = (time() - os.path.getmtime(errlog))/3600
else:
    errage = 24*360#no file, no error!
if age<agelim and errage>agelim:
    print "all is well: log age = ",age,"last error ",errage,"hours ago"
    sys.exit(0)


warning = "--------%s-------\ncn_rx_log_check.py: correlator log %s is %4.2f hours old\nlast error in %s was %4.1f hours ago"%(asctime(),
logfile,
age,
errlog,
errage
)
open(logfile,'a').write(warning)


msg = MIMEMultipart()


To = "wheresmytab@gmail.com"
From = "teampaper@gmail.com"



msg['Subject'] = "WARNING: Stale Correlator log"
msg['From'] = From
msg['To'] = To
msg.attach(MIMEText(warning))

#os.system('python ~/bin/plot_all_temps.py `ls -tr /home/obs/Temperature_Logs/2*txt | tail -n 48`')
#msg.attach(MIMEImage(open('/home/obs/Temperature_Logs/all_temps.png').read()))
msg.attach(MIMEText("/var/log/messages follows\n"+'-'*10))
msg.attach(MIMEText(open('/var/log/messages').read()))
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
