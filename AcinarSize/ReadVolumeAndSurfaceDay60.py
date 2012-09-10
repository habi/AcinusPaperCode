# -*- coding: utf8 -*-

from pylab import *
import os
import glob
import csv
import numpy as np
import xlrd # from http://www.lexicon.net/sjmachin/xlrd.htm, to read XLS-Files (like the ones containing the lung volumes)
from matplotlib2tikz import save as tikz_save

print "Hey ho, let's go"

Drive = 'R:\SLS'

if os.path.exists(Drive) == False:
	print 'Cannot read ' + str(Drive) + '. Exiting!'
	exit()

chatty = False # set to 'False' to suppress some output.
PlotTheData = True # True/False switches between showing and not showing the plots at the end
SaveFigures = True # save the plot to .png-Files
TikZTheData = True # save the data to .tikz-Files to import into LaTeX
TOMCATVoxelSize = 1.48
SliceNumber = 10
ShrinkageFactor = 0.6
color=['c','r','m','b','y'] # colors for Volumeplots
STEPanizerDir = 'voxelsize'+str(TOMCATVoxelSize)+'-every'+str(SliceNumber)+'slice'
print 'Calculating with a Shrinkagefactor of',str(ShrinkageFactor) +'x'
print '---'

WhichRat = 'R108C60'
Rat = ['A','B_B1_mrg','C','Dt-mrg','Et-mrg']
Beamtime = ['','2010c','','2009f\mrg','2009f\mrg'] # obviously no Beamtime for A and C.

# According to http://stackoverflow.com/a/2397192 instead of 
	#~ bar = []
	#~ for item in some_iterable:
		#~ bar.append(SOME EXPRESSION)
#~ # one should (and as seen below can) use
	#~ bar = [SOME EXPRESSION for item in some_iterable]

# Reading Volume Data of RUL from Stefans Data-File
print 'Reading values from p:\doc\#Dev\AcinarSize\Datenblattstefan.xls'
Datenblattstefan = xlrd.open_workbook('p:\doc\#Dev\AcinarSize\Datenblattstefan.xls').sheet_by_index(0)
RULVolume = [ Datenblattstefan.cell(70+int(Sample),9).value for Sample in range(len(Rat)) ]
AbsoluteParenchymalVolume = [ Datenblattstefan.cell(70+int(Sample),15).value for Sample in range(len(Rat)) ]
AbsoluteAirspaceVolume = [ Datenblattstefan.cell(70+int(Sample),19).value for Sample in range(len(Rat)) ]
AbsoluteAirspaceSurface = [ Datenblattstefan.cell(70+int(Sample),94).value for Sample in range(len(Rat)) ]

for Sample in range(len(Rat)):
	print WhichRat+Rat[Sample]+': RUL volume =',RULVolume[Sample],'cm^3'
	print WhichRat+Rat[Sample]+': Absolute parenchymal volume =',AbsoluteParenchymalVolume[Sample],'cm^2'
	print WhichRat+Rat[Sample]+': Absolute airspace volume =',AbsoluteAirspaceVolume[Sample],'cm^2'
	print WhichRat+Rat[Sample]+': Absolute airspace surface =',AbsoluteAirspaceSurface[Sample],'cm^2'

print ''
print 'DatenblattStefan.xls:'
print 'Mean RUL volume is',np.mean(RULVolume),'cm^3.'
print 'Mean absolute parenchymal volume is',np.mean(AbsoluteParenchymalVolume),'cm^2.'
print 'Mean absolute airspace volume is',np.mean(AbsoluteAirspaceVolume),'cm^2.'
print 'Mean absolute airspace surface is',np.mean(AbsoluteAirspaceSurface),'cm^2.'
print '---'

# Reading Data from MeVisLab and STEPanizer files
print 'Reading data from the MeVisLab .dcm and STEPanizer .xls files on',Drive
print 'See which samples exist and how many acini we assessed'

# See how many .csv-Files we actually have in the STEPanizer directories
SamplePath = {}
for Sample in range(len(Rat)):
	SamplePath[Sample] = os.path.join(Drive,Beamtime[Sample],WhichRat + Rat[Sample])
	if os.path.exists(SamplePath[Sample]):
		print SamplePath[Sample],'exists and contains',len(glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir))),\
		'directories with images generated for STEPanizering'
	else:
		print SamplePath[Sample],'does not exist, we probably did not scan this rat.'

# Generating a list of the .csv filenames
CSVFile = [glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir,'*.xls')) for Sample in range(len(Rat))]	
print '---'

# See for which Rat we counted the most acini. We need this for scaling the plots and pre-allocating empty variables
ListOfAcini=[]
for Sample in range(len(Rat)):
	for CurrentFile in CSVFile[Sample][:]:
		CurrentAcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
		ListOfAcini.append(CurrentAcinusNumber)
MaximumAcini = max(ListOfAcini) + 1
TotalAssessedAcini = len(ListOfAcini)

# MeVisLab ->
# Read Volumes from .dcm-Files generated with MeVisLab
print 'MeVisLab'
print 'Extracting acinar volumes from .dcm-Filenames written with MeVisLab'
MeVisLabVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]

for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print 'Sample',WhichRat + Rat[Sample] + ': not measured'
	else:
		print 'Sample',Beamtime[Sample] + '\R108C60' + Rat[Sample] + ': Reading volumes from',len(glob.glob(os.path.join(SamplePath[Sample],'*.dcm'))),'.dcm-Files'
		for CurrentFile in sorted(glob.glob(os.path.join(SamplePath[Sample],'*.dcm'))):
			Acinus = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('.volume')])
			MeVisLabVolume[Sample][Acinus] = float(CurrentFile[CurrentFile.find('volume')+len('volume'):CurrentFile.find('.pixelsize')])/1000 # normalize ul to cm^3: http://is.gd/XxU3Ei
			if chatty:
				print 'MeVisLab-Volume of acinus',Acinus,'is',MeVisLabVolume[Sample][Acinus],'cm^3'
print ''

print 'Mean acinar volume'
MeanMeVisLabVolume = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		MeanMeVisLabVolume[Sample] = np.mean(np.ma.masked_array(MeVisLabVolume[Sample],np.isnan(MeVisLabVolume[Sample])))
		print WhichRat + Rat[Sample] + ':', MeanMeVisLabVolume[Sample], 'cm^3'
print 'Mean acinar volume for all samples:', np.mean(np.ma.masked_array(MeanMeVisLabVolume,np.isnan(MeanMeVisLabVolume))), 'cm^3'
print 'Standard deviation of the mean acinar volume for all samples:', np.std(np.ma.masked_array(MeanMeVisLabVolume,np.isnan(MeanMeVisLabVolume)))
print '---'
# <- MeVisLab

# STEPanizer ->
print 'STEPanizer'
# Read data from each STEPanizer .csv-file and calculate the desired values
AcinarVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
SurfaceDensity = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
AbsoluteSurface = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
TotalTestPoints = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
NonParenchymalPoints= [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
ParenchymalPoints = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
ParenchymalVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print 'Sample',WhichRat + Rat[Sample] + ': not measured'
	else:
		print 'Sample',Beamtime[Sample] + '\R108C60' + Rat[Sample] + ': Reading and calculating',len(glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerDir))),'sets of values'
		for CurrentFile in sorted(CSVFile[Sample][:]):
			Acinus = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
			TotalSlices = len(glob.glob(os.path.join((os.path.join(SamplePath[Sample],'acinus'+ str("%02d" % (Acinus)),STEPanizerDir)),'*.jpg')))
			FileData = csv.reader(open(CurrentFile,'rb'),dialect=csv.excel_tab)
			for line in FileData:
				if len(line) > 0:
					if line[1] == 'Interception with Tissue':
						Interceptions = int(line[3])
					if line[1] == 'Point inside Acinus':
						AcinusTestPoints = int(line[3])
					if line[1] == 'Non-Parenchymal Points':
						NonParenchymalPoints[Sample][Acinus] = int(line[3])
					if line[0] == 'Pixel size:':
						STEPanizerPixelSize = double(line[1])
					if line[0] == 'a(p):':
						Area = double(line[1])*STEPanizerPixelSize**2
					if line[0] == 'l(p):':
						LinePointLength = double(line[1])*STEPanizerPixelSize
					if line[0] == 'Number of test points:':
						AllAreaTestPoints = int(line[1])
			# Give out counted/assessed data if desired
			if chatty:
				print 'Acinus', "%02d" % (Acinus,),\
					'|',"%03d" % (TotalSlices),'Files',\
					'|',AllAreaTestPoints,'test points per image',\
					'|',"%04d" % (Interceptions),'Interceptions',\
					'|',"%03d" % (AcinusTestPoints),'Points in Acinus',\
					'|',"%03d" % (NonParenchymalPoints[Sample][Acinus]),'Nonparenchymal Points'

			# Volume = AcinusTestPoints * Area * STEPanizerPixelSize * SliceNumber * TOMCATVoxelSize		
			AcinarVolume[Sample][Acinus] = ((( AcinusTestPoints * Area * STEPanizerPixelSize * SliceNumber * TOMCATVoxelSize ) / ShrinkageFactor) / 1e12) # scaling volume to cm^3: http://is.gd/wbZ81O
								
			# Total of all points = AllAreaTestPoints (in file) * Total of slices
			TotalTestPoints[Sample][Acinus] = AllAreaTestPoints * TotalSlices
			
			# Parenchym-Volume = Parenchymal points * Area * STEPanizerPixelSize * SliceNumber * TOMCATVoxelSize. Parenchymal points are total points minus NParenchymalpoints (which are all points outside the sample and in non-parenchyma)
			ParenchymalPoints[Sample][Acinus] = TotalTestPoints[Sample][Acinus] - NonParenchymalPoints[Sample][Acinus]
			ParenchymalVolume[Sample][Acinus] = ParenchymalPoints[Sample][Acinus] * Area * STEPanizerPixelSize * SliceNumber * TOMCATVoxelSize
			
			# Surface Density = 2 * Interceptions / Length
			Length = ( LinePointLength * AcinusTestPoints ) / 1e4 # Linepointlength has to be calculated for the acinar reference space, /1e4 to scale from um to cm, http://is.gd/7wu1UD
			SurfaceDensity[Sample][Acinus] = 2 * Interceptions / Length
			
			# Absolute Surface = Surface density * acinar volume
			AbsoluteSurface[Sample][Acinus] = SurfaceDensity[Sample][Acinus] * AcinarVolume[Sample][Acinus]
if ShrinkageFactor != 1:
	print 'Calculated with a Shrinkagefactor of',ShrinkageFactor,'x'			
print ''

# Give out interesting values
print 'Mean acinar volume for single samples'
MeanAcinarVolume = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		MeanAcinarVolume[Sample] = np.mean(np.ma.masked_array(AcinarVolume[Sample],np.isnan(AcinarVolume[Sample])))
		print WhichRat + Rat[Sample] + ':', np.round(MeanAcinarVolume[Sample],decimals=3), 'cm^3'
print 'Mean acinar volume for all samples:', np.round(np.mean(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume))),decimals=3), 'cm^3'
print 'Standard deviation of the mean acinar volume for all samples:', np.std(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume)))
print ''

print 'Mean acinar surface for single samples'
MeanAcinarSurface = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		MeanAcinarSurface[Sample] = np.mean(np.ma.masked_array(AbsoluteSurface[Sample],np.isnan(AbsoluteSurface[Sample])))
		print WhichRat + Rat[Sample] + ':', np.round(MeanAcinarSurface[Sample],decimals=3), 'cm^2'		
print 'Mean acinar surface for all samples:',np.round(np.mean(np.ma.masked_array(MeanAcinarSurface,np.isnan(MeanAcinarSurface))),decimals=3), 'cm^2'
print 'Standard deviation of the mean acinar surface for all samples:', np.std(np.ma.masked_array(MeanAcinarSurface,np.isnan(MeanAcinarSurface)))
print ''

# Absolute parenchymal Volume / mean acinar volume = Number of Acini
print 'Number of acini (= absolute airspace volume from stefan / mean acinar volume)'
NumberOfAcini = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ' was not measured'
	else:
		NumberOfAcini[Sample] = AbsoluteAirspaceVolume[Sample] / (MeanAcinarVolume[Sample])
		print WhichRat + Rat[Sample],'contains',int(np.round(NumberOfAcini[Sample])),'acini'
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
		DiffusionSurface[Sample] = NumberOfAcini[Sample] * MeanAcinarSurface[Sample]
		print WhichRat + Rat[Sample],'contains',np.round(DiffusionSurface[Sample],decimals=3),'cm^2 of diffusion surface'
print ''

print 'Stefan measured the absolute airspace surface with EM and got'
for Sample in range(len(Rat)):
	print WhichRat + Rat[Sample] +':',np.round(AbsoluteAirspaceSurface[Sample],decimals=3),'cm^2'
print ''
print 'Stefans mean absolute airspace surface is',np.round(np.mean(AbsoluteAirspaceSurface),decimals=3),'cm^2.'
print 'Our mean airspace surface is',np.round(np.mean(np.ma.masked_array(DiffusionSurface,np.isnan(DiffusionSurface))),decimals=3),'cm^2.'
#~ # <- STEPanizer

print ''
print 'MeVisLab volume compared to STEPanizer volume (STEPanizer/MeVisLab)'
STEPanizerMeVisLabDifference = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] != '':
		for Acinus in range(MaximumAcini):
			if isnan(MeVisLabVolume[Sample][Acinus]) == False:
				STEPanizerMeVisLabDifference[Sample][Acinus] = np.round((AcinarVolume[Sample][Acinus] / MeVisLabVolume[Sample][Acinus]),decimals=3)
				if chatty:
					print WhichRat + Rat[Sample],'Acinus',str(Acinus) + ':',STEPanizerMeVisLabDifference[Sample][Acinus]
		if chatty:
			print '(a mean of',np.round(np.mean(np.ma.masked_array(STEPanizerMeVisLabDifference[Sample],np.isnan(STEPanizerMeVisLabDifference[Sample]))),decimals=3),')'
	else:
		if chatty:
			print WhichRat + Rat[Sample],'was not measured'
	
print 'In mean, (STEPanizer-/MeVisLab-volume) is',np.round(np.mean(np.ma.masked_array(STEPanizerMeVisLabDifference,np.isnan(STEPanizerMeVisLabDifference))),decimals=3)

print ''
print 'Values for acinus.tex:'
print 'add to preamble around line 76.'
print '\\newcommand{\\numberofacini}{' + str(TotalAssessedAcini) + '}'
print '\\newcommand{\\volume}{' + str(np.mean(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume)))) + '} % cm^3, (mean acinar volume)'
print '\\newcommand{\std}{' + str(np.std(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume)))) + '} % (Standard deviation of acinar volumes)'
print '\\newcommand{\difference}{' + str(+np.mean(np.ma.masked_array(STEPanizerMeVisLabDifference,np.isnan(STEPanizerMeVisLabDifference)))) + '} % X times bigger (acinar volumes STEPanizer/MeVisLab-volumes)'

if PlotTheData == False:
	exit()

if TikZTheData:
	print "Using 'plt.plot' instead of 'plt.scatter', since 'tikz_save('file.tikz')' doesn't work otherwise"
	print "just add 'only marks' to the TikZ-code (in the \\addplot-command)"

# plot MeVisLabVolumes in one plot
plt.figure()
legend=[]
for Sample in range(len(Rat)):
	legend.append(WhichRat+Rat[Sample])
	plt.plot(range(MaximumAcini*Sample,MaximumAcini*(Sample+1)),MeVisLabVolume[Sample],c=color[Sample])
plt.title('MeVisLabVolume')
plt.legend(legend)
if SaveFigures:
	plt.savefig('plot_mevisvolumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_mevisvolumes.tikz')

# plot STEPanizerVolumes in one plot
plt.figure()
for Sample in range(len(Rat)):
	plt.plot(range(MaximumAcini*Sample,MaximumAcini*(Sample+1)),AcinarVolume[Sample],c=color[Sample])
plt.title('STEPanizerVolume')
plt.legend(legend)
if SaveFigures:
	plt.savefig('plot_stepanizervolumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_stepanizervolumes.tikz')

		
# Plot the interesting stuff
# Plot MeVisLab- against STEPanizer-volumes
plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(len(Rat),1,Sample+1)
	if Beamtime[Sample] != '':
		if TikZTheData:
			plt.plot(range(MaximumAcini),MeVisLabVolume[Sample],c='g')
			plt.plot(range(MaximumAcini),AcinarVolume[Sample],c='b')
		else:
			plt.scatter(range(MaximumAcini),MeVisLabVolume[Sample],c='g')
			plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c='b')
			plt.legend(['MeVisLab','STEPanizer'])
	if Beamtime[Sample] == '':
		plt.title('Volumes of acini of ' + WhichRat + Rat[Sample] + ' were not assessed')
	else:
		plt.title('Volumes of acini of ' + Beamtime[Sample] + '\\' + WhichRat+Rat[Sample])
	plt.xlim([0,MaximumAcini])
plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_mevis_vs_stepanizer_volumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_mevis_vs_stepanizer_volumes.tikz')


plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(len(Rat),1,Sample+1)
	if TikZTheData:
		plt.plot(range(MaximumAcini),STEPanizerMeVisLabDifference[Sample],c=color[Sample])
	else:
		plt.scatter(range(MaximumAcini),STEPanizerMeVisLabDifference[Sample],c=color[Sample])		
	plt.xlim([0,MaximumAcini])
	plt.ylim([0,13])
	locs, labels = plt.xticks()
	if Beamtime[Sample] != '':
		plt.boxplot([x for x in STEPanizerMeVisLabDifference[Sample] if not math.isnan(x)],1,positions=[MaximumAcini/2],widths=5) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	if Beamtime[Sample] == '':
		plt.title(WhichRat + Rat[Sample])
	else:
		plt.title(Beamtime[Sample] + '\\' + WhichRat + Rat[Sample] + '\nSTEPanizer/MeVisLab-Volume')
	plt.xticks(locs)
plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_mevis_vs_stepanizer_volumeratio.png',transparent=False)
if TikZTheData:
	tikz_save('plot_mevis_vs_stepanizer_volumeratio.tikz')

# Boxplot of Volumes
plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	#~ plt.xlim([0,2])
	plt.subplot(3,len(Rat),Sample+1)#,axisbg=color[Sample])
	if Beamtime[Sample] != '':
		plt.boxplot([x for x in AcinarVolume[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.title(str(WhichRat) + str(Rat[Sample])+'\nRUL-Volume='+str(RULVolume[Sample])+' cm^3')
	plt.ylabel('Acinar Volume')
	plt.subplot(3,len(Rat),Sample+1+len(Rat))#,axisbg=color[Sample])
	if Beamtime[Sample] != '':
		plt.boxplot([x for x in SurfaceDensity[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.ylabel('Surface Density')
	plt.subplot(3,len(Rat),Sample+1+2*len(Rat))#,axisbg=color[Sample])
	if Beamtime[Sample] != '':
		plt.boxplot([x for x in AbsoluteSurface[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	plt.ylabel('Absolute Surface')
plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_boxplot_of_volumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_boxplot_of_volumes.tikz')

# Plotting Volumes, Surface Density and Absolute Surfaces
plt.figure(figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(311)
	plt.xlim([0,MaximumAcini])
	if Beamtime[Sample] != '':
		if TikZTheData:
			plt.plot(range(MaximumAcini),AcinarVolume[Sample],c=color[Sample])
		else:
			plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c=color[Sample])
	plt.title('Volume (AcinusTestPoints * Area * STEPanizerPixelSize * SliceNumber * TOMCATVoxelSize)')
	plt.subplot(312)
	plt.xlim([0,MaximumAcini])
	if Beamtime[Sample] != '':
		if TikZTheData:
			plt.plot(range(MaximumAcini),SurfaceDensity[Sample],c=color[Sample])
		else:
			plt.scatter(range(MaximumAcini),SurfaceDensity[Sample],c=color[Sample])
	plt.title('Surface Density (2 * Int / Length)')
	plt.subplot(313)
	plt.xlim([0,MaximumAcini])
	plt.ylim([0,0.0006])
	if Beamtime[Sample] != '':
		if TikZTheData:
			plt.plot(range(MaximumAcini),AbsoluteSurface[Sample],c=color[Sample])
		else:
			plt.scatter(range(MaximumAcini),AbsoluteSurface[Sample],c=color[Sample])
	plt.title('Absolute Surface (SurfaceDensity * AcinarVolume)')
	plt.legend([Beamtime[1] + "\\" + WhichRat + Rat[1],Beamtime[3] + '\\' + WhichRat + Rat[3],Beamtime[4] + '\\' + WhichRat + Rat[4]],\
		loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3)
if SaveFigures:		
	plt.savefig('plot_acinarvolume_surfacedensity_absolutesurface.png',transparent=False)
if TikZTheData:
	tikz_save('plot_acinarvolume_surfacedensity_absolutesurface.tikz')

plt.show()
