#!/usr/bin/python
'''
Module to email users from airglow.monitor@gmail.com

History: 24 Jul 2013 - initial script written

Written by Daniel J. Fisher (dfisher2@illinois.edu)
'''

# Import required modules
import smtplib


def email(recipient,subject,msg):
    '''
    Summary
    -------
        email(recipient,subject,msg)
        emails someone a something

    Inputs
    ------
        recipient = subject of email
        subject = subject of email
        msg = body text of email

    History
    -------
        7/17/13 -- Written by DJF (dfisher2@illionis.edu)
    '''
    # Email Sender
    sender = 'airglow.monitor@gmail.com'
    sender_pass = 'sB6vhcDH'
    # Send the email
    headers = ["From: " + sender,  "Subject: " + subject, "To: " + recipient, "MIME-Version: 1.0", "Content-Type: text/plain"]  # NO HTML so /n works
    headers = "\r\n".join(headers)
    email = smtplib.SMTP('smtp.gmail.com', 587)
    email.ehlo()
    email.starttls()
    email.ehlo
    email.login(sender,sender_pass)
    email.sendmail(sender,recipient, headers + "\r\n\r\n" + msg)
    email.quit()
    print("!!! Email Sent")

def emailerror(recipients,subject,msg):
    '''
    Summary:
        emailerror(subject,recipients,msg)
        emails recipients a warning

    Inputs:
        recipients = list of email addresses
        subject    = subject of email
        msg        = body text of email

    History:
        7/17/13 -- Written by DJF (dfisher2@illionis.edu)
        3/31/16 -- Added email list as input
    '''
    # Append more addresses here for warnings
    for to in recipients:
        email(to,subject,msg)

#def emailerror(subject,msg):
#    '''
#    Summary
#    -------
#        emailerror(subject,msg)
#        emails Haardvark & Danimal a warning
#
#    Inputs
#    ------
#        subject = subject of email
#        msg = body text of email
#
#    History
#    -------
#        7/17/13 -- Written by DJF (dfisher2@illionis.edu)
#    '''
#    # Append more addresses here for warnings
#    recipient = ['bharding512@gmail.com', 'fishnchips1624+status@gmail.com']
#    for to in recipient:
#        email(to,subject,msg)

