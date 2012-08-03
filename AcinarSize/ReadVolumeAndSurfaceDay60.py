from pylab import *
import os
import glob
import csv
import numpy as np

Drive = 'R:\SLS'

Rat = ['B_B1_mrg','Dt-mrg','Et-mrg']
Beamtime = ['2010c','2009f\mrg','2009f\mrg']
VoxelSize = 1.48
SliceDistance = 10
STEPanizerDir = 'voxelsize'+str(VoxelSize)+'-every'+str(SliceDistance)+'slice'

SamplePath = {}

for Sample in range(len(Rat)):
	SamplePath[Sample] = os.path.join(Drive,Beamtime[Sample],'R108C60' + Rat[Sample])
	if os.path.exists(SamplePath[Sample]):
		print SamplePath[Sample],'does exist'
	else:
		print SamplePath[Sample],'does not exist'
	print 'and contains',len(glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir))),\
		'directories with acini'
print '---'

# According to http://stackoverflow.com/a/2397192 instead of 
	#~ bar = []
	#~ for item in some_iterable:
		#~ bar.append(SOME EXPRESSION)
#~ # one should (and as seen below can) use
	#~ bar = [SOME EXPRESSION for item in some_iterable]

print 'Finding all *.xls-Files in'
for Sample in range(len(Rat)):
	print SamplePath[Sample],'R108C60' + Rat[Sample]
CSVFile = [glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir,'*.xls')) for Sample in range(len(Rat))]
print '---'

# Find maximum number of counted acini
an=[]
for Sample in range(len(Rat)):
	for CurrentFile in CSVFile[Sample][:]:
		AcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
		an.append(AcinusNumber)
MaximumAcini = max(an) +1

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
		AcinarVolume[Sample][AcinusNumber] = Points * Area * VoxelSize * SliceDistance
		print 'The volume of acinus',AcinusNumber,'is',Points,'*',int(np.round(Area)),'*',VoxelSize,'*',SliceDistance,'i.e.',int(AcinarVolume[Sample][AcinusNumber]),'um^3'
		SurfaceDensity[Sample][AcinusNumber] = 2 * Interceptions / Length
		#~ print 'The SurfaceDensity of acinus',AcinusNumber,'is 2 *',Interceptions,'*',int(Length)
		AbsoluteSurface[Sample][AcinusNumber] = AcinarVolume[Sample][AcinusNumber] * SurfaceDensity[Sample][AcinusNumber]
print '---'

color=['r','g','y']

plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(3,3,Sample+1)#,axisbg=color[Sample])
	plt.boxplot([x for x in AcinarVolume[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.title('R108C60' + Rat[Sample])
	plt.ylabel('Acinar Volume')
	plt.ylim([0,1e9])
	plt.subplot(3,3,Sample+1+3)#,axisbg=color[Sample])
	plt.boxplot([x for x in SurfaceDensity[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.ylabel('Surface Density')
	plt.ylim([0,8e5])
	plt.subplot(3,3,Sample+1+6)#,axisbg=color[Sample])
	plt.boxplot([x for x in AbsoluteSurface[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.ylabel('Absolute Surface')
	plt.ylim([0,7e14])
plt.savefig('boxplots.png',transparent=False)
plt.draw()

plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(311)
	plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c=color[Sample])
	plt.legend(['R108C60' + Rat[0],'R108C60' + Rat[1],'R108C60' + Rat[2]])	
	plt.title('Volume (Points * Area * Voxelsize * Slicedistance)')
	plt.xlim([0,MaximumAcini])
	plt.subplot(312)
	plt.scatter(range(MaximumAcini),SurfaceDensity[Sample],c=color[Sample])
	plt.legend(['R108C60' + Rat[0],'R108C60' + Rat[1],'R108C60' + Rat[2]])	
	plt.title('Surface Density (2 * Interceptions * Length)')
	plt.xlim([0,MaximumAcini])
	plt.subplot(313)
	plt.scatter(range(MaximumAcini),AbsoluteSurface[Sample],c=color[Sample])
	plt.legend(['R108C60' + Rat[0],'R108C60' + Rat[1],'R108C60' + Rat[2]])	
	plt.title('Absolute Surface (absolute Volume * Volume Density * Surface Density)')
	plt.xlim([0,MaximumAcini])
	plt.ylim([0,7e14])
plt.savefig('volume.png',transparent=False)
plt.draw()

plt.show()
	
# Save Data
AcinusWriter = csv.writer(open('acini.csv', 'wb'), delimiter=';')
AcinusWriter.writerow(range(MaximumAcini))
for Sample in range(len(Rat)):
	AcinusWriter.writerow(AcinarVolume[Sample][:])
for Sample in range(len(Rat)):
	AcinusWriter.writerow(SurfaceDensity[Sample][:])
for Sample in range(len(Rat)):
	AcinusWriter.writerow(AbsoluteSurface[Sample][:])
