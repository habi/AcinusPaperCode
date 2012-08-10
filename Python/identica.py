# Loads notices from identi.ca via the API, saves the dates of the notices and plots them

import urllib2
from xml.dom.minidom import parseString # from http://is.gd/IQUOew
import datetime
from pylab import *

howmany = 50

# Load freshest notice ID
try:
	file = urllib2.urlopen('http://identi.ca/api/statuses/public_timeline.xml')
	data = file.read()
	file.close()
	dom = parseString(data)
	Tag = dom.getElementsByTagName('id')[0].toxml()
	ID=int(Tag.replace('<id>','').replace('</id>',''))
except:
	print 'Bonkers, I cannot load the public timeline. Giving up!'
	exit()

# Load 'howmany' notices up to the freshest one and plot their date
print 'I will load',howmany,'notices from the first to ID',ID
print '---'
plt.figure()
counter = 1
for Notice in range(1,ID,int(ID/howmany)):
	try:
		file = urllib2.urlopen('http://identi.ca/api/statuses/show/'+str(Notice)+'.xml')
		data = file.read()
		file.close()
		dom = parseString(data)
		Tag = dom.getElementsByTagName('created_at')[0].toxml()
		#strip off the tag (<tag>data</tag>  --->   data):
		Date=str(Tag.replace('<created_at>','').replace('</created_at>',''))
		Date = str(Date[:-11])+str(Date[-5:])
		Date = datetime.datetime.strptime(Date,'%a %b %d %H:%M:%S %Y')
		print '[' + str(counter) + '/' + str(howmany) + '] ID',Notice,'was posted on',Date
		plt.plot_date(Date,Notice)
		plt.title('Identica Notice Number vs. Date')
	except:
		# Try the next one if this one doesnt exist
		try:
			Notice += 1
			file = urllib2.urlopen('http://identi.ca/api/statuses/show/'+str(Notice)+'.xml')
			data = file.read()
			file.close()
			dom = parseString(data)
			Tag = dom.getElementsByTagName('created_at')[0].toxml()
			#strip off the tag (<tag>data</tag>  --->   data):
			Date=str(Tag.replace('<created_at>','').replace('</created_at>',''))
			Date = str(Date[:-11])+str(Date[-5:])
			Date = datetime.datetime.strptime(Date,'%a %b %d %H:%M:%S %Y')
			print '[' + str(counter) + '/' + str(howmany) + '] ID',Notice,'was posted on',Date
			plt.plot_date(Date,Notice)
			plt.title('Identica Notice Number vs. Date')
		except:
			print '[' + str(counter) + '/' + str(howmany) + '] ID',ID,'does not exist, sorry!'
	counter += 1

plt.show()
