#!/usr/bin/python
# Filename: identica.py
# loads defined notes from identi.ca and saves then to "timeline.txt" for further perusal

import urllib
import string

for counter in range(1,4000000,1000) : 
	url = "http://identi.ca/notice/" + str(counter)
	print "loading " + str(url)
	html = urllib.urlopen(url)
	while 1:
		htmlline = html.readline()
		if htmlline == "":
   			break
		if (htmlline.find("<title>") > -1) :
			htmlline = htmlline.replace("<title>"," ")
			htmlline = htmlline.replace("</title>"," ")
			htmlline = htmlline.replace("UTC - Identi.ca"," ")
			htmlline = htmlline.replace("'s status on ","")
			htmlline = htmlline.replace("Monday,",",")
			htmlline = htmlline.replace("Tuesday,",",")
			htmlline = htmlline.replace("Wednesday,",",")
			htmlline = htmlline.replace("Thursday,",",")
			htmlline = htmlline.replace("Friday,",",")
			htmlline = htmlline.replace("Saturday,",",")
			htmlline = htmlline.replace("Sunday,",",")
			htmlline = htmlline.replace("-08 ","-08, ")
			htmlline = htmlline.replace("-09 ","-09, ")
			htmlline = htmlline.strip(" ")
			print "saving title of " + str(url) + " to timeline.csv"
			timeline = open('timeline.csv', 'a')
			timeline.write(str(counter) + ", " + htmlline)
timeline.close()