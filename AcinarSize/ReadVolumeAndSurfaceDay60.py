from pylab import *
import os
import glob
import csv
import numpy as np
import xlrd # from http://www.lexicon.net/sjmachin/xlrd.htm, to read XLS-Files (like the ones containing the lung volumes)

Drive = 'R:\SLS'

VoxelSize = 1.48
SliceDistance = 10

WhichRat = 'R108C60'
Rat = ['B_B1_mrg','Dt-mrg','Et-mrg']
Beamtime = ['2010c','2009f\mrg','2009f\mrg']
SliceDistance = 10

#~ WhichRat = 'R108C36'
#~ Rat = ['Dt-mrg']
#~ Beamtime = ['2010a\mrg']
#~ VoxelSize = 2.96
#~ SliceDistance = 4

# According to http://stackoverflow.com/a/2397192 instead of 
	#~ bar = []
	#~ for item in some_iterable:
		#~ bar.append(SOME EXPRESSION)
#~ # one should (and as seen below can) use
	#~ bar = [SOME EXPRESSION for item in some_iterable]

# Reading Volume Data of RUL from Stefans Data-File
XLSFile = xlrd.open_workbook('p:\doc\#R\AcinusPaper\Datenblattstefan.xls')
WorkSheet = XLSFile.sheet_by_index(5) # load worksheet 6 (remember, python starts at 0)
#~ RULVolume = [ WorkSheet.cell(24+int(Sample),9).value for Sample in range(5) ] # load all Samples
#~ Name = ['A','B','C','D','E'] # load all Samples
RULVolume = [ WorkSheet.cell(24+int(Sample),9).value for Sample in [1,3,4] ] # only load B,D and E
Name = ['B','D','E'] # only load B,D and E
for Sample in range(len(Name)):
	print 'The RUL volume of',WhichRat+Name[Sample],'is',RULVolume[Sample],'cm^3'
print '---'

# Reading Data from STEPanizer XLS-Files (which actually are just .csv)
STEPanizerDir = 'voxelsize'+str(VoxelSize)+'-every'+str(SliceDistance)+'slice'

# See how many .csv-Files we actually have
SamplePath = {}

for Sample in range(len(Rat)):
	SamplePath[Sample] = os.path.join(Drive,Beamtime[Sample],WhichRat + Rat[Sample])
	if os.path.exists(SamplePath[Sample]):
		print SamplePath[Sample],'does exist'
	else:
		print SamplePath[Sample],'does not exist'
	print 'and contains',len(glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir))),\
		'directories with acini'
print '---'

# Generating a list of the filenames
CSVFile = [glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir,'*.xls')) for Sample in range(len(Rat))]

# Find maximum number of counted acini
tmp=[]
for Sample in range(len(Rat)):
	for CurrentFile in CSVFile[Sample][:]:
		AcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
		tmp.append(AcinusNumber)
MaximumAcini = max(tmp) +1

# Read necessary data from each STEPanizer .csv-file, calculate volume of acinus
AcinarVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
SurfaceDensity = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
AbsoluteSurface = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	print 'reading and calculating values for R108C60' + Rat[Sample]
	for CurrentFile in CSVFile[Sample][:]:
		AcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
		FileData = csv.reader(open(CurrentFile,'rb'),dialect=csv.excel_tab)
		for line in FileData:
			if len(line) > 0:
				if line[1] == 'Interception with Tissue':
					Interceptions = int(line[3])
				if line[1] == 'Point inside Acinus':
					Points = int(line[3])
				if line[0] == 'Pixel size:':
					PixelSize = double(line[1])
				if line[0] == 'a(p):':
					Area = double(line[1])*PixelSize**2
				if line[0] == 'l(p):':
					Length = double(line[1])*PixelSize
		#~ print 'Acinus',str(AcinusNumber) + ':',Interceptions,'interceptions,',Points,'points inside & area a(p) of',int(np.round(Area)),'um^2'
		# Volume = Points * Area * Voxelsize * Slicedistance		
		AcinarVolume[Sample][AcinusNumber] = Points * Area * VoxelSize * SliceDistance
		#~ print 'The volume of acinus',AcinusNumber,'is',Points,'*',int(np.round(Area)),'*',VoxelSize,'*',SliceDistance,'i.e.',int(AcinarVolume[Sample][AcinusNumber]),'um^3'
			
		# Surface Density = 2 * Interceptions * Length
		SurfaceDensity[Sample][AcinusNumber] = 2 * Interceptions / Length
		#~ print 'The SurfaceDensity of acinus',AcinusNumber,'is 2 *',Interceptions,'*',int(Length)
		
		# Absolute Surface = absolute Volume * Volume Density * Surface Density
		AbsoluteSurface[Sample][AcinusNumber] = AcinarVolume[Sample][AcinusNumber] * SurfaceDensity[Sample][AcinusNumber]
		##############################
		# How and where can I get the volume density?
		##############################
print '---'

color=['r','b','y']

plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(3,3,Sample+1)#,axisbg=color[Sample])
	plt.boxplot([x for x in AcinarVolume[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.title(str(WhichRat) + str(Rat[Sample])+' (RUL-Volume='+str(RULVolume[Sample])+' cm^3)')
	plt.ylabel('Acinar Volume')
	plt.ylim([0,1e9])
	plt.subplot(3,3,Sample+1+3)#,axisbg=color[Sample])
	plt.boxplot([x for x in SurfaceDensity[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.ylabel('Surface Density')
	plt.ylim([0,40])
	plt.subplot(3,3,Sample+1+6)#,axisbg=color[Sample])
	plt.boxplot([x for x in AbsoluteSurface[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.ylabel('Absolute Surface')
	plt.ylim([0,3e10])
plt.savefig('boxplots.png',transparent=False)
plt.draw()

plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(311)
	plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c=color[Sample])
	plt.legend([WhichRat + Rat[0],WhichRat + Rat[1],WhichRat + Rat[2]])	
	plt.title('Volume (Points * Area * Voxelsize * Slicedistance)')
	plt.xlim([0,MaximumAcini])
	plt.subplot(312)
	plt.scatter(range(MaximumAcini),SurfaceDensity[Sample],c=color[Sample])
	plt.legend([WhichRat + Rat[0],WhichRat + Rat[1],WhichRat + Rat[2]])	
	plt.title('Surface Density (2 * Interceptions * Length)')
	plt.xlim([0,MaximumAcini])
	plt.subplot(313)
	plt.scatter(range(MaximumAcini),AbsoluteSurface[Sample],c=color[Sample])
	plt.legend([WhichRat + Rat[0],WhichRat + Rat[1],WhichRat + Rat[2]])	
	plt.title('Absolute Surface (absolute Volume * Volume Density * Surface Density)')
	plt.xlim([0,MaximumAcini])
	plt.ylim([0,3e10])
plt.savefig('volume.png',transparent=False)
plt.draw()

plt.show()

# Save Data
AcinusWriter = csv.writer(open('acini.csv', 'wb'), delimiter=';')
tmp = range(MaximumAcini)
tmp.insert(0,'Acinus')
AcinusWriter.writerow(tmp)
for Sample in range(len(Rat)):
	AcinarVolume[Sample].insert(0,'AcinarVolume ' + WhichRat + Rat[Sample])
	AcinusWriter.writerow(AcinarVolume[Sample][:])
for Sample in range(len(Rat)):
	SurfaceDensity[Sample].insert(0,'SurfaceDensity ' + WhichRat + Rat[Sample])
	AcinusWriter.writerow(SurfaceDensity[Sample][:])
for Sample in range(len(Rat)):
	AbsoluteSurface[Sample].insert(0,'AbsoluteSurface ' + WhichRat + Rat[Sample])
	AcinusWriter.writerow(AbsoluteSurface[Sample][:])
