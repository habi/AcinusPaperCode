from pylab import *
import os
import glob
import csv
import numpy as np
import xlrd # from http://www.lexicon.net/sjmachin/xlrd.htm, to read XLS-Files (like the ones containing the lung volumes)

Drive = 'R:\SLS'

if os.path.exists(Drive) == False:
	print 'Cannot read ' + str(Drive) + '. Exiting!'
	exit()

VoxelSize = 1.48
SliceNumber = 10

WhichRat = 'R108C60'
Rat = ['A','B_B1_mrg','C','Dt-mrg','Et-mrg']
Beamtime = ['','2010c','','2009f\mrg','2009f\mrg'] # obviously no Beamtime for A and C.
SliceNumber = 10

# According to http://stackoverflow.com/a/2397192 instead of 
	#~ bar = []
	#~ for item in some_iterable:
		#~ bar.append(SOME EXPRESSION)
#~ # one should (and as seen below can) use
	#~ bar = [SOME EXPRESSION for item in some_iterable]

# Reading Volume Data of RUL from Stefans Data-File
print 'Reading values from p:\doc\#Dev\AcinarSize\Datenblattstefan.xls'
ShrinkageFactor = 1
print 'Calculating everything with a Shrinkagefactor of',str(ShrinkageFactor) +'x'
XLSFile = xlrd.open_workbook('p:\doc\#Dev\AcinarSize\Datenblattstefan.xls')
WorkSheet = XLSFile.sheet_by_index(0)
RULVolume = [ WorkSheet.cell(70+int(Sample),9).value*ShrinkageFactor for Sample in range(len(Rat)) ]
AbsoluteParenchymalVolume = [ WorkSheet.cell(70+int(Sample),15).value*ShrinkageFactor for Sample in range(len(Rat)) ]
AbsoluteAirspaceSurface = [ WorkSheet.cell(70+int(Sample),94).value*ShrinkageFactor for Sample in range(len(Rat)) ]

for Sample in range(len(Rat)):
	print WhichRat+Rat[Sample]+':'
	print 'RUL volume of',RULVolume[Sample],'cm^3 (Shrinkagefactor:',str(ShrinkageFactor) +'x)'
	print 'Absolute parenchymal volume of',AbsoluteParenchymalVolume[Sample],'cm^2 (Shrinkagefactor:',str(ShrinkageFactor) +'x)'
	print 'Absolute airspace surface of',AbsoluteAirspaceSurface[Sample],'cm^2 (Shrinkagefactor:',str(ShrinkageFactor) +'x)'
if ShrinkageFactor != 1:
	print '-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'
	print ' Check the calculation of the values above, since ShrinkageFactor != 1! '
	print '     Stefan thinks he did not use the shrinkage for the surface...      '
	print '-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'
print '---'
print 'The mean RUL volume is',np.mean(RULVolume),'cm^3.'
print 'The mean absolute parenchymal volume is',np.mean(AbsoluteParenchymalVolume),'cm^2.'
print 'The mean absolute airspace surface is',np.mean(AbsoluteAirspaceSurface),'cm^2.'
print '---'

# Reading Data from STEPanizer XLS-Files (which actually are just .csv)
STEPanizerDir = 'voxelsize'+str(VoxelSize)+'-every'+str(SliceNumber)+'slice'

# See how many .csv-Files we actually have
SamplePath = {}
for Sample in range(len(Rat)):
	SamplePath[Sample] = os.path.join(Drive,Beamtime[Sample],WhichRat + Rat[Sample])
	if os.path.exists(SamplePath[Sample]):
		print SamplePath[Sample],'exists and contains',len(glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir))),\
		'acinus directories'
	else:
		print SamplePath[Sample],'does not exist, we probably did not scan this rat.'
print '---'

# Generating a list of the filenames
CSVFile = [glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir,'*.xls')) for Sample in range(len(Rat))]

# See for which Rat we counted the most acini. We need this for scaling the plot afterwards
tmp=[]
for Sample in range(len(Rat)):
	for CurrentFile in CSVFile[Sample][:]:
		AcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
		tmp.append(AcinusNumber)
MaximumAcini = max(tmp) + 1

# Read Volumes from .dcm-Files generated with MeVisLab
MeVisLabVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print 'R108C60' + Rat[Sample] + ' was not measured'
	else:
		print Beamtime[Sample] + '\R108C60' + Rat[Sample] + ': Reading volumes from .dcm-Files'
		for CurrentFile in sorted(glob.glob(os.path.join(SamplePath[Sample],'*.dcm'))):
			AcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('.volume')])
			MeVisLabVolume[Sample][AcinusNumber] = float(CurrentFile[CurrentFile.find('volume')+len('volume'):CurrentFile.find('.pixelsize')])
			print 'MeVisLab-Volume of acinus',AcinusNumber,'is',MeVisLabVolume[Sample][AcinusNumber]

# Read necessary data from each STEPanizer .csv-file, calculate volume of acinus
AcinarVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
SurfaceDensity = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
AbsoluteSurface = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
TotalTestPoints = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
NonParenchymalPoints= [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
ParenchymalPoints = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
ParenchymalVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print 'R108C60' + Rat[Sample] + ' was not measured'
	else:
		print Beamtime[Sample] + '\R108C60' + Rat[Sample] + ': Reading and calculating',len(glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir))),'sets of values'
		for CurrentFile in sorted(CSVFile[Sample][:]):
			AcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
			TotalSlices = len(glob.glob(os.path.join((os.path.join(SamplePath[Sample],'acinus'+ str("%02d" % (AcinusNumber)),STEPanizerDir)),'*.jpg')))
			FileData = csv.reader(open(CurrentFile,'rb'),dialect=csv.excel_tab)
			for line in FileData:
				if len(line) > 0:
					if line[1] == 'Interception with Tissue':
						Interceptions = int(line[3])
					if line[1] == 'Point inside Acinus':
						AcinusTestPoints = int(line[3])
					if line[1] == 'Non-Parenchymal Points':
						NParenchymalPoints = int(line[3])
					if line[0] == 'Pixel size:':
						PixelSize = double(line[1])
					if line[0] == 'a(p):':
						Area = double(line[1])*PixelSize**2
					if line[0] == 'l(p):':
						LinePointLength = double(line[1])*PixelSize
					if line[0] == 'Number of test points:':
						AllAreaTestPoints = int(line[1])
			# Give out counted/assessed data
			print 'Acinus', "%02d" % (AcinusNumber,),\
				'|',"%03d" % (TotalSlices),'Files',\
				'|',AllAreaTestPoints,'test points per image',\
				'|',"%04d" % (Interceptions),'Interceptions',\
				'|',"%03d" % (AcinusTestPoints),'Points in Acinus',\
				'|',"%03d" % (NParenchymalPoints),'Nonparenchymal Points'

			#~ print 'Acinus',str(AcinusNumber) + ':',Interceptions,'interceptions,',TestPoints,'points inside & area a(p) of',int(np.round(Area)),'um^2'
			# Volume = AcinusTestPoints * Area * PixelSize * SliceNumber * VoxelSize		
			AcinarVolume[Sample][AcinusNumber] = AcinusTestPoints * Area * PixelSize * SliceNumber * VoxelSize # PixelSize == STEPanizer, VoxelSize == TOMCAT
			#~ print 'The volume of acinus',AcinusNumber,'is',AcinusTestPoints,'*',int(np.round(Area)),'*',PixelSize,'*',SliceNumber,'i.e.',int(AcinarVolume[Sample][AcinusNumber]),'um^3'
								
			# Total of all points = AllAreaTestPoints (in file) * Total of slices
			TotalTestPoints[Sample][AcinusNumber] = AllAreaTestPoints * TotalSlices
			
			# Parenchym-Volume = Parenchymal points * Area * PixelSize * SliceNumber * VoxelSize. Parenchymal points are total points minus NParenchymalpoints (which are all points outside the sample and in non-parenchyma)
			ParenchymalPoints[Sample][AcinusNumber] = TotalTestPoints[Sample][AcinusNumber] - NParenchymalPoints
			NonParenchymalPoints[Sample][AcinusNumber] = NParenchymalPoints
			ParenchymalVolume[Sample][AcinusNumber] = ParenchymalPoints[Sample][AcinusNumber] * Area * PixelSize * SliceNumber * VoxelSize
			
			# Surface Density = 2 * Interceptions / Length
			Length = LinePointLength * AcinusTestPoints # Linepointlength has to be calculated for the acinar reference space.
			SurfaceDensity[Sample][AcinusNumber] = 2 * Interceptions / Length
			
			# Absolute Surface = absolute Volume * Volume Density * Surface Density
			AbsoluteSurface[Sample][AcinusNumber] = SurfaceDensity[Sample][AcinusNumber] * AcinarVolume[Sample][AcinusNumber]
print '---'

color=['c','r','m','b','y']

# Plot MeVisLab- against STEPanizer-volumes
plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(len(Rat),1,Sample+1)
	plt.scatter(range(MaximumAcini),MeVisLabVolume[Sample],c='r')
	plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c='b')
	plt.legend(['MeVisLab','STEPanizer'])
	plt.title(Beamtime[Sample]+WhichRat+Rat[Sample])
	subplots_adjust(hspace=0.5)
	

plt.show()
exit()

# Boxplot of Volumes
plt.figure(figsize=(16,9))
subplots_adjust(hspace=0.01)
for Sample in range(len(Rat)):
	plot1 = plt.subplot(3,len(Rat),Sample+1)#,axisbg=color[Sample])
	if Beamtime[Sample] != '':
		plt.boxplot([x for x in AcinarVolume[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.title(str(WhichRat) + str(Rat[Sample])+'\nRUL-Volume='+str(RULVolume[Sample])+' cm^3')
	plt.ylabel('Acinar Volume')
	plt.ylim([0,6e9])
	plot2 = plt.subplot(3,len(Rat),Sample+1+len(Rat))#,axisbg=color[Sample])
	if Beamtime[Sample] != '':
		plt.boxplot([x for x in SurfaceDensity[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.ylabel('Surface Density')
	plt.ylim([0,0.09])
	plot3 = plt.subplot(3,len(Rat),Sample+1+2*len(Rat))#,axisbg=color[Sample])
	if Beamtime[Sample] != '':
		plt.boxplot([x for x in AbsoluteSurface[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.ylabel('Absolute Surface')
	plt.ylim([0,2.5e8])
	xticklabels = plot1.get_xticklabels() + plot2.get_xticklabels() + plot3.get_xticklabels()
	setp(xticklabels, visible=False)
plt.savefig('boxplot.png',transparent=False)
plt.draw()

# Plotting Volumes, Surface Density and Absolute Surfaces
plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plot1 = plt.subplot(311)
	if Beamtime[Sample] != '':
		plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c=color[Sample])
	plt.title('Volume (AcinusTestPoints * Area * PixelSize * SliceNumber * VoxelSize)')
	plt.xlim([0,MaximumAcini])
	plot1 = plt.subplot(312)
	if Beamtime[Sample] != '':
		plt.scatter(range(MaximumAcini),SurfaceDensity[Sample],c=color[Sample])
	plt.title('Surface Density (2 * Int / Length)')
	plt.xlim([0,MaximumAcini])
	plot1 = plt.subplot(313)
	if Beamtime[Sample] != '':
		plt.scatter(range(MaximumAcini),AbsoluteSurface[Sample],c=color[Sample])
	plt.title('Absolute Surface (SurfaceDensity * AcinarVolume)')
	plt.xlim([0,MaximumAcini])
	plt.legend([Beamtime[1] + "\\" + WhichRat + Rat[1],Beamtime[3] + '\\' + WhichRat + Rat[3],Beamtime[4] + '\\' + WhichRat + Rat[4]],\
		loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3)
subplots_adjust(hspace=0.25)
plt.savefig('volume.png',transparent=False)
plt.draw()

print 'Mean acinar volume for single samples'
MeanAcinarVolume = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		MeanAcinarVolume[Sample] = np.mean(np.ma.masked_array(AcinarVolume[Sample],np.isnan(AcinarVolume[Sample])))
		print WhichRat + Rat[Sample] + ':', MeanAcinarVolume[Sample]/1e12, 'cm^3'
print 'Mean acinar volume for all samples:', np.mean(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume)))/1e12, 'cm^3'
print '---'

print 'Mean acinar surface for single samples'
MeanAcinarSurface = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		MeanAcinarSurface[Sample] = np.mean(np.ma.masked_array(AbsoluteSurface[Sample],np.isnan(AbsoluteSurface[Sample])))
		print WhichRat + Rat[Sample] + ':', MeanAcinarSurface[Sample]/1e8, 'cm^2'
print 'Mean acinar surface for all samples:', np.mean(np.ma.masked_array(MeanAcinarSurface,np.isnan(MeanAcinarSurface)))/1e8, 'cm^2'
print '---'

# NonParenchymal absolute Volume / mean acinar volume = Number of Acini
print 'Number of acini (= nonparenchymal absolute volume from stefan / mean acinar volume)'
NumberOfAcini = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ' was not measured'
	else:
		NumberOfAcini[Sample] = int(np.round(AbsoluteParenchymalVolume[Sample] / (MeanAcinarVolume[Sample]/1e12)))
		print WhichRat + Rat[Sample],'contains',NumberOfAcini[Sample],'acini'
print ''		
print 'Rodriguez1987 (page 146) states a total of 4023 acini for the whole rat lungs.'
print 'We have a *mean* of',int(np.round(np.mean(np.ma.masked_array(NumberOfAcini,np.isnan(NumberOfAcini))))),'acini calculated for the whole lung.'
print '---'

# Number of Acini * Mean acinar Surface = Diffusionsurface
print 'Diffusion surface (=number of acini * mean acinar surface)'
DiffusionSurface = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ' was not measured'
	else:
		DiffusionSurface[Sample] = np.round(NumberOfAcini[Sample] * MeanAcinarSurface[Sample]/1e8,decimals=3)
		print WhichRat + Rat[Sample],'contains',DiffusionSurface[Sample],'cm^2 of diffusion surface'
print ''
print 'Stefan measured the absolute airspace surface with EM and got'
for Sample in range(len(Rat)):
	print WhichRat + Rat[Sample] +':',np.round(AbsoluteAirspaceSurface[Sample],decimals=3),'cm^2'
print ''
print 'Stefans mean absolute airspace surface is',np.round(np.mean(AbsoluteAirspaceSurface),decimals=3),'cm^2.'
print 'Our mean airspace surface is',np.round(np.mean(np.ma.masked_array(DiffusionSurface,np.isnan(DiffusionSurface))),decimals=3),'cm^2.'

plt.show()
