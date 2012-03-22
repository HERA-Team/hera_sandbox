#!/usr/bin/env python
"""
Send an email to danny when the correlator dies
"""

import sys,os,statvfs,aipy as a
from glob import glob
import smtplib
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from time import time

msg = MIMEMultipart()


To = "wheresmytab@gmail.com"
From = "teampaper@gmail.com"



msg['Subject'] = "The Correlator is Dead. Long live the Correlator"
msg['From'] = From
msg['To'] = To

message = '\n'.join(open('/home/obs/logs/roachlog.txt').readlines()[-10:])
msg.attach(MIMEText(message))
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
