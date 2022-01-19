"""
A scraper to pull emails en-masse from my Gmail account. Designed to grab 
GeoNet Moment Tensor emails sent by John Ristau from GNS
"""
import os
import email
import getpass
import imaplib

# user = input("Gmail Username: ")
user = "bryant.h.c@gmail.com"
pwd = getpass.getpass("Gmail Password: ")

m = imaplib.IMAP4_SSL("imap.gmail.com")
m.login(user, pwd)

# Selecting the specific label in my inbox
resp, count = m.select("POSTHOC")

# Simply grab all emails in this inbox which has been pre-sorted
resp, uids = m.uid("search", None, "ALL")

event_ids = []
for uid in uids[0].split()[::-1]:
    resp, data = m.fetch(uid, "(UID BODY[TEXT])")
    lines = str(data[0][1])[:3000].split("\\r\\n")  # trial and errored this 
    for line in lines:
        if "EVENT ID" in line:
            event_ids.append(line.split(" ")[-1])
            break
    else:
        print(f"No ID found {str(uid)}")


print(event_ids)

