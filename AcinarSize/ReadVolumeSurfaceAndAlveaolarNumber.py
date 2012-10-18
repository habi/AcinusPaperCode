# -*- coding: utf8 -*-

from pylab import *
import os
import glob
import csv
import numpy as np
import xlrd # from http://www.lexicon.net/sjmachin/xlrd.htm, to read XLS-Files (like the ones containing the lung volumes)
from matplotlib2tikz import save as tikz_save
import time

print "Hey ho, let's go"

Drive = 'R:\SLS'

if os.path.exists(Drive) == False:
	print 'Cannot read ' + str(Drive) + '. Exiting!'
	exit()

chatty = False # set to 'False' to suppress some output.
PlotTheData = True # True/False switches between showing and not showing the plots at the end
SaveFigures = False # save the plot to .png-Files
TikZTheData = False # save the data to .tikz-Files to import into LaTeX
TOMCATVoxelSize = 1.48
SliceNumber = 10
DisectorThickness = 5 # sklices
ShrinkageFactor = 0.6
color=['c','r','m','b','y'] # colors for Volumeplots
STEPanizerVolumeDir = 'voxelsize'+str(TOMCATVoxelSize)+'-every'+str(SliceNumber)+'slice'
STEPanizerAlveoliDir = 'voxelsize'+str(TOMCATVoxelSize)+'-every'+str(SliceNumber)+'slice-DisectorThickness-' + str('%.2f' % round(DisectorThickness*TOMCATVoxelSize,2)) + 'um-or' + str(DisectorThickness) + 'slices'

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

print
print 'DatenblattStefan.xls:'
print 'Mean RUL volume is',np.mean(RULVolume),'cm^3.'
print 'Mean absolute parenchymal volume is',np.mean(AbsoluteParenchymalVolume),'cm^2.'
print 'Mean absolute airspace volume is',np.mean(AbsoluteAirspaceVolume),'cm^2.'
print 'Mean absolute airspace surface is',np.mean(AbsoluteAirspaceSurface),'cm^2.'
print '________________________________________________________________________________'

# Reading Data from MeVisLab and STEPanizer files
print 'Reading data from the MeVisLab .dcm and STEPanizer .xls files on',Drive
print
print 'See which samples exist and how many acini we assessed:'

# See how many .csv-Files we actually have in the STEPanizer directories
SamplePath = {}
AssessedAcini = [0 for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	SamplePath[Sample] = os.path.join(Drive,Beamtime[Sample],WhichRat + Rat[Sample])
	if os.path.exists(SamplePath[Sample]):
		AssessedAcini[Sample] = len(glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerVolumeDir)))
		print SamplePath[Sample],'exists and contains',AssessedAcini[Sample],'directories with images for STEPanizering.'
	else:
		print SamplePath[Sample],'does not exist.'

# Generating a list of the .csv filenames
CSVFileVolume = [glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerVolumeDir,'*.xls')) for Sample in range(len(Rat))]	
CSVFileAlveoli = [glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerAlveoliDir,'*.xls')) for Sample in range(len(Rat))]	
print '________________________________________________________________________________'

# See for which Rat we counted the most acini. We need this for scaling the plots and pre-allocating empty variables
ListOfAciniVolume=[]
ListOfAciniAlveoli=[]
MaximumAcini=[]
TotalAssessedAcini=[]
for Sample in range(len(Rat)):
	for CurrentFile in CSVFileVolume[Sample][:]:
		CurrentAcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
		ListOfAciniVolume.append(CurrentAcinusNumber)
	for CurrentFile in CSVFileVolume[Sample][:]:
		CurrentAcinusNumber = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
		ListOfAciniAlveoli.append(CurrentAcinusNumber)
MaximumAcini.append(max(ListOfAciniVolume) + 1)
MaximumAcini.append(max(ListOfAciniAlveoli) + 1)
TotalAssessedAcini.append(len(ListOfAciniVolume))
TotalAssessedAcini.append(len(ListOfAciniAlveoli))

if TotalAssessedAcini[0] != TotalAssessedAcini[1]:
	print 'Total number of assessed acini for volume (',TotalAssessedAcini[0],') and alveoli (',TotalAssessedAcini[1],') does not match!'
	print 'Print see what is wrong and restart'
	sys.exit()
else:
	tmp = TotalAssessedAcini[0]
	TotalAssessedAcini = None
	TotalAssessedAcini = tmp
		
if MaximumAcini[0] != MaximumAcini[1]:
	print 'Maximum number of assessed acini for volume (',MaximumAcini[0],') and alveoli (',MaximumAcini[1],') does not match!'
	print 'Print see what is wrong and restart'
	sys.exit()
else:
	tmp = MaximumAcini[0]
	MaximumAcini = None
	MaximumAcini = tmp

# MeVisLab ->
# Read Volumes from .dcm-Files generated with MeVisLab
print 'MeVisLab'
print '  Extracting acinar volumes from .dcm-Filenames written with MeVisLab:'
MeVisLabVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		print Beamtime[Sample] + '\R108C60' + Rat[Sample] + ': Reading volumes from',AssessedAcini[Sample],'.dcm files'
		for CurrentFile in sorted(glob.glob(os.path.join(SamplePath[Sample],'*.dcm'))):
			Acinus = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('.volume')])
			MeVisLabVolume[Sample][Acinus] = float(CurrentFile[CurrentFile.find('volume')+len('volume'):CurrentFile.find('.pixelsize')])/1000 # normalize ul to cm^3: http://is.gd/XxU3Ei
			if chatty:
				print 'MeVisLab-Volume of acinus',Acinus,'is',MeVisLabVolume[Sample][Acinus],'cm^3'
print

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

# STEPanizer ->
print 'STEPanizer: Volume (David)'
print '  Extracting acinar volumes from .xls files exported from STEPanizer:'
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
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		print Beamtime[Sample] + '\R108C60' + Rat[Sample] + ': Reading',AssessedAcini[Sample],'.xls files'
		for CurrentFile in sorted(CSVFileVolume[Sample][:]):
			Acinus = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
			TotalSlices = len(glob.glob(os.path.join((os.path.join(SamplePath[Sample],'acinus'+ str("%02d" % (Acinus)),STEPanizerVolumeDir)),'*.jpg')))
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
						STEPanizerPixelSize_Vol = double(line[1])
					if line[0] == 'a(p):':
						Area_Vol = double(line[1])*STEPanizerPixelSize_Vol**2
					if line[0] == 'l(p):':
						LinePointLength_Vol = double(line[1])*STEPanizerPixelSize_Vol
					if line[0] == 'Number of test points:':
						AllAreaTestPoints_Vol = int(line[1])
			# Give out counted/assessed data if desired
			if chatty:
				print 'Acinus', "%02d" % (Acinus),\
					'|',"%03d" % (TotalSlices),'Files',\
					'|',AllAreaTestPoints_Vol,'test point(s) per image',\
					'|',"%04d" % (Interceptions),'Interceptions',\
					'|',"%03d" % (AcinusTestPoints),'Points in Acinus',\
					'|',"%03d" % (NonParenchymalPoints[Sample][Acinus]),'Nonparenchymal Points'

			# Volume = AcinusTestPoints * Area_Vol * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize		
			AcinarVolume[Sample][Acinus] = ((( AcinusTestPoints * Area_Vol * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize ) / ShrinkageFactor) / 1e12) # scaling volume to cm^3: http://is.gd/wbZ81O
								
			# Total of all points = AllAreaTestPoints_Vol (in file) * Total of slices
			TotalTestPoints[Sample][Acinus] = AllAreaTestPoints_Vol * TotalSlices
			
			# Parenchym-Volume = Parenchymal points * Area_Vol * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize. Parenchymal points are total points minus NParenchymalpoints (which are all points outside the sample and in non-parenchyma)
			ParenchymalPoints[Sample][Acinus] = TotalTestPoints[Sample][Acinus] - NonParenchymalPoints[Sample][Acinus]
			ParenchymalVolume[Sample][Acinus] = ParenchymalPoints[Sample][Acinus] * Area_Vol * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize
			
			# Surface Density = 2 * Interceptions / Length
			Length = ( LinePointLength_Vol * AcinusTestPoints ) / 1e4 # Linepointlength has to be calculated for the acinar reference space, /1e4 to scale from um to cm, http://is.gd/7wu1UD
			SurfaceDensity[Sample][Acinus] = 2 * Interceptions / Length
			
			# Absolute Surface = Surface density * acinar volume
			AbsoluteSurface[Sample][Acinus] = SurfaceDensity[Sample][Acinus] * AcinarVolume[Sample][Acinus]
if ShrinkageFactor != 1:
	print 'Calculated with a Shrinkagefactor of',ShrinkageFactor,'x'			
print

print 'STEPanizer: Alveoli (Evelyne)'
# Read data from each STEPanizer .csv-file and calculate the desired values
Bridges = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print 'Sample',WhichRat + Rat[Sample] + ': not measured'
	else:
		print 'Sample',Beamtime[Sample] + '\R108C60' + Rat[Sample] + ': Reading',AssessedAcini[Sample],'.xls files.'
		for CurrentFile in sorted(CSVFileAlveoli[Sample][:]):
			Acinus = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
			TotalSlices = len(glob.glob(os.path.join((os.path.join(SamplePath[Sample],'acinus'+ str("%02d" % (Acinus)),STEPanizerAlveoliDir)),'*.jpg')))
			FileData = csv.reader(open(CurrentFile,'rb'),dialect=csv.excel_tab)
			for line in FileData:
				if len(line) > 0:
					if line[0] == 'Num 1':
						Counts = int(line[3])
					if line[0] == 'Pixel size:':
						STEPanizerPixelSize_Acini = double(line[1])
					if line[0] == 'a(p):':
						Area_Acini = double(line[1])*STEPanizerPixelSize_Vol**2
					if line[0] == 'l(p):':
						LinePointLength_Alveoli = double(line[1])*STEPanizerPixelSize_Vol
			# Give out counted/assessed data if desired
			if chatty:
				print 'Acinus ' + str("%02d" % (Acinus)) +\
					' | ' + str("%03d" % (TotalSlices)) + ' Files' +\
					' | ' + str(Counts) + ' counts'
			
			# Bridges = Something to do with the Counts
			Bridges[Sample][Acinus] = Counts

if ShrinkageFactor != 1:
	print 'Calculated with a Shrinkagefactor of',ShrinkageFactor,'x'			
print

print '________________________________________________________________________________'

# Give out interesting values
print 'Mean acinar volume for single samples'
MeanAcinarVolume = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		MeanAcinarVolume[Sample] = np.mean(np.ma.masked_array(AcinarVolume[Sample],np.isnan(AcinarVolume[Sample])))
		print WhichRat + Rat[Sample] + ':', MeanAcinarVolume[Sample], 'cm^3'
print 'Mean acinar volume for all samples:', np.round(np.mean(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume))),decimals=3), 'cm^3'
print 'Standard deviation of the mean acinar volume for all samples:', np.std(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume)))
print

print 'Mean acinar surface for single samples'
MeanAcinarSurface = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not measured'
	else:
		MeanAcinarSurface[Sample] = np.mean(np.ma.masked_array(AbsoluteSurface[Sample],np.isnan(AbsoluteSurface[Sample])))
		print WhichRat + Rat[Sample] + ':', MeanAcinarSurface[Sample], 'cm^2'
print 'Mean acinar surface for all samples:',np.round(np.mean(np.ma.masked_array(MeanAcinarSurface,np.isnan(MeanAcinarSurface))),decimals=3), 'cm^2'
print 'Standard deviation of the mean acinar surface for all samples:', np.std(np.ma.masked_array(MeanAcinarSurface,np.isnan(MeanAcinarSurface)))
print

# Absolute parenchymal Volume / mean acinar volume = Number of Acini
print 'Number of acini (= absolute airspace volume from stefan / mean acinar volume)'
NumberOfAcini = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ' was not measured'
	else:
		NumberOfAcini[Sample] = AbsoluteAirspaceVolume[Sample] / (MeanAcinarVolume[Sample])
		print WhichRat + Rat[Sample],'contains',int(np.round(NumberOfAcini[Sample])),'acini'
print		
print 'Rodriguez1987 (page 146) states a total of 4023 acini for the whole rat lungs.'
print 'We have a *mean* of',int(np.round(np.mean(np.ma.masked_array(NumberOfAcini,np.isnan(NumberOfAcini))))),'acini calculated for the whole lung.'

# Number of Acini * Mean acinar Surface = Diffusionsurface
print 'Diffusion surface (=number of acini * mean acinar surface)'
DiffusionSurface = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ' was not measured'
	else:
		DiffusionSurface[Sample] = NumberOfAcini[Sample] * MeanAcinarSurface[Sample]
		print WhichRat + Rat[Sample],'contains',np.round(DiffusionSurface[Sample],3),'cm^2 of diffusion surface'
print

print 'Stefan measured the absolute airspace surface with EM and got'
for Sample in range(len(Rat)):
	print WhichRat + Rat[Sample] +':',np.round(AbsoluteAirspaceSurface[Sample],decimals=3),'cm^2'
print
print 'Stefans mean absolute airspace surface is',np.round(np.mean(AbsoluteAirspaceSurface),decimals=3),'cm^2.'
print 'Our mean airspace surface is',np.round(np.mean(np.ma.masked_array(DiffusionSurface,np.isnan(DiffusionSurface))),decimals=3),'cm^2.'
#~ # <- STEPanizer

print
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

print '________________________________________________________________________________'
print
print 'Values for acinus.tex:'
print 'add to preamble around line 76.'
print '  \\newcommand{\\numberofacini}{' + str(TotalAssessedAcini) + '}'
print '  \\newcommand{\\volume}{' + str(np.mean(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume)))) + '} % cm^3, (mean acinar volume)'
print '  \\newcommand{\std}{' + str(np.std(np.ma.masked_array(MeanAcinarVolume,np.isnan(MeanAcinarVolume)))) + '} % (Standard deviation of acinar volumes)'
print '  \\newcommand{\difference}{' + str(+np.mean(np.ma.masked_array(STEPanizerMeVisLabDifference,np.isnan(STEPanizerMeVisLabDifference)))) + '} % X times bigger (acinar volumes STEPanizer/MeVisLab-volumes)'

if PlotTheData == False:
	exit()

if TikZTheData:
	print "Using 'plt.plot' instead of 'plt.scatter', since 'tikz_save('file.tikz')' doesn't work otherwise"
	print "just add 'only marks' to the TikZ-code (in the \\addplot-command)"

PlotMeVisLabVolume = []
PlotSTEPAnizerVolume = []
for Sample in range(len(Rat)):
	for Acinus in range(MaximumAcini):
		if isnan(MeVisLabVolume[Sample][Acinus]) == False:
			PlotMeVisLabVolume.append(MeVisLabVolume[Sample][Acinus])
			PlotSTEPAnizerVolume.append(AcinarVolume[Sample][Acinus])
			
PlotNormalizedMeVisLabVolume = [float(i)/max(PlotMeVisLabVolume) for i in PlotMeVisLabVolume]
PlotNormalizedSTEPanizerVolume = [float(i)/max(PlotSTEPAnizerVolume) for i in PlotSTEPAnizerVolume]

# Plot MeVisLab- and STEPanizer-Volumes (both original and normalized
plt.figure(num=None,figsize=(16,9))
# plot MeVisLabVolumes
plt.subplot(221)
# Set a counter, so we can iterate through the number of Assessedacini for plotting.
# For each plot, we would like to plot from AssessedAcini[Sample-1]:AssesseAcini[Sample],
# but this does not work, since the array has only a length of 5. We thus
# accumulate the value with "counter = counter + AssessedAcini[Sample]" and use this for plotting
counter = 0
for Sample in range(len(Rat)):
	plt.plot(range(counter,counter+AssessedAcini[Sample]),PlotMeVisLabVolume[counter:counter+AssessedAcini[Sample]],c=color[Sample])
	counter = counter + AssessedAcini[Sample]
plt.ylabel('Volume [cm^3]')
plt.xlabel('Acinus')
plt.legend([WhichRat+Rat[1],WhichRat+Rat[3],WhichRat+Rat[4]],loc='best')
plt.title('MeVisLabVolume')
plt.xticks(arange(TotalAssessedAcini))
# plot STEPanizerVolumes
plt.subplot(222)
counter = 0
for Sample in range(len(Rat)):
	plt.plot(range(counter,counter+AssessedAcini[Sample]),PlotSTEPAnizerVolume[counter:counter+AssessedAcini[Sample]],c=color[Sample])
	counter = counter + AssessedAcini[Sample]
plt.ylabel('Volume [cm^3]')
plt.xlabel('Acinus')
plt.legend([WhichRat+Rat[1],WhichRat+Rat[3],WhichRat+Rat[4]],loc='best')
plt.title('STEPanizerVolumes')
plt.xticks(arange(TotalAssessedAcini))
# plot NORMALIZED MeVisLabVolumes
plt.subplot(223)
counter = 0
for Sample in range(len(Rat)):
	plt.plot(range(counter,counter+AssessedAcini[Sample]),PlotNormalizedMeVisLabVolume[counter:counter+AssessedAcini[Sample]],c=color[Sample])
	counter = counter + AssessedAcini[Sample]
plt.ylabel('Normalized Volume')
plt.xlabel('Acinus')
plt.legend([WhichRat+Rat[1],WhichRat+Rat[3],WhichRat+Rat[4]],loc='best')
plt.title('NORMALIZED MeVisLabVolume')
plt.xticks(arange(TotalAssessedAcini))
# plot NORMALIZED STEPanizerVolumes without NaNs
plt.subplot(224)
counter = 0
for Sample in range(len(Rat)):
	plt.plot(range(counter,counter+AssessedAcini[Sample]),PlotNormalizedSTEPanizerVolume[counter:counter+AssessedAcini[Sample]],c=color[Sample])
	counter = counter + AssessedAcini[Sample]
plt.ylabel('Normalized Volume')
plt.xlabel('Acinus')
plt.legend([WhichRat+Rat[1],WhichRat+Rat[3],WhichRat+Rat[4]],loc='best')
plt.title('NORMALIZED STEPanizerVolumes')
plt.xticks(arange(TotalAssessedAcini))
plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_mevislab_vs_stepanizervolumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_mevislab_vs_stepanizervolumes.tikz')

# plot MeVisLabVolumes in one plot
plt.figure(num=None,figsize=(16,9))
legend=[]
for Sample in range(len(Rat)):
	legend.append(WhichRat+Rat[Sample])
	plt.plot(range(MaximumAcini*Sample,MaximumAcini*(Sample+1)),MeVisLabVolume[Sample],c=color[Sample])
plt.title('MeVisLabVolume')
plt.legend(legend,loc='best')
if SaveFigures:
	plt.savefig('plot_mevisvolumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_mevisvolumes.tikz')

# plot STEPanizerVolumes in one plot
plt.figure(num=None,figsize=(16,9))
for Sample in range(len(Rat)):
	plt.plot(range(MaximumAcini*Sample,MaximumAcini*(Sample+1)),AcinarVolume[Sample],c=color[Sample])
plt.title('STEPanizerVolume')
plt.legend(legend,loc='best')
if SaveFigures:
	plt.savefig('plot_stepanizervolumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_stepanizervolumes.tikz')
			
# Plot the interesting stuff
# Plot MeVisLab- against STEPanizer-volumes
plt.figure(num=None,figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(len(Rat),1,Sample+1)
	if Beamtime[Sample] != '':
		if TikZTheData:
			plt.plot(range(MaximumAcini),MeVisLabVolume[Sample],c='g')
			plt.plot(range(MaximumAcini),AcinarVolume[Sample],c='b')
		else:
			plt.scatter(range(MaximumAcini),MeVisLabVolume[Sample],c='g')
			plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c='b')
			plt.legend(['MeVisLab','STEPanizer'],loc='best')
	if Beamtime[Sample] == '':
		plt.title('Volumes of acini of ' + WhichRat + Rat[Sample] + ' were not assessed')
	else:
		plt.title('Volumes of acini of ' + Beamtime[Sample] + '\\' + WhichRat+Rat[Sample])
	plt.xlim([0,MaximumAcini])
	plt.ylim([0,None])
	plt.xticks(arange(MaximumAcini))
plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_mevis_vs_stepanizer_volumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_mevis_vs_stepanizer_volumes.tikz')

plt.figure(num=None,figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(len(Rat),1,Sample+1)
	if TikZTheData:
		plt.plot(range(MaximumAcini),STEPanizerMeVisLabDifference[Sample],c=color[Sample])
	else:
		plt.scatter(range(MaximumAcini),STEPanizerMeVisLabDifference[Sample],c=color[Sample])		
	plt.xlim([0,MaximumAcini])
	plt.ylim([0,None])
	if Beamtime[Sample] == '':
		plt.title(WhichRat + Rat[Sample])
	else:
		plt.title(Beamtime[Sample] + '\\' + WhichRat + Rat[Sample] + ': STEPanizer/MeVisLab-Difference')
	plt.xticks(arange(MaximumAcini))
plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_mevis_vs_stepanizer_volumeratio.png',transparent=False)
if TikZTheData:
	tikz_save('plot_mevis_vs_stepanizer_volumeratio.tikz')

# Boxplot of Volumes
plt.figure(num=None,figsize=(16,9))
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
plt.figure(num=None,figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(311)
	plt.xlim([0,MaximumAcini])
	if Beamtime[Sample] != '':
		if TikZTheData:
			plt.plot(range(MaximumAcini),AcinarVolume[Sample],c=color[Sample])
		else:
			plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c=color[Sample])
	plt.title('Volume (AcinusTestPoints * Area * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize)')
	plt.xlim([0,MaximumAcini])
	plt.ylim([0,None])
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
	plt.ylim([0,None])
	if Beamtime[Sample] != '':
		if TikZTheData:
			plt.plot(range(MaximumAcini),AbsoluteSurface[Sample],c=color[Sample])
		else:
			plt.scatter(range(MaximumAcini),AbsoluteSurface[Sample],c=color[Sample])
	plt.title('Absolute Surface (SurfaceDensity * AcinarVolume)')
	plt.legend([Beamtime[1] + "\\" + WhichRat + Rat[1],Beamtime[3] + '\\' + WhichRat + Rat[3],Beamtime[4] + '\\' + WhichRat + Rat[4]],loc='best')#loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3)
	plt.xlim([0,MaximumAcini])
	plt.ylim([0,None])
	plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_acinarvolume_surfacedensity_absolutesurface.png',transparent=False)
if TikZTheData:
	tikz_save('plot_acinarvolume_surfacedensity_absolutesurface.tikz')
	
	
plt.show()
exit()	





	

# Plot Bridges vs. Volume (Evelyne vs. David)
plt.figure(num=None,figsize=(16,9))
plt.subplots_adjust(hspace=1)
for Sample in range(len(Rat)):
	ax = plt.subplot(len(Rat),1,Sample+1)
	plt.scatter(range(MaximumAcini),Bridges[Sample],c='r')
	plt.scatter(range(MaximumAcini),AcinarVolume[Sample],c='b')
	plt.title('Bridges ' + WhichRat + Rat[Sample] +': ' +\
		str(np.sum(np.ma.masked_array(Bridges[Sample],np.isnan(Bridges[Sample])))) +\
		' (total)')
	# Shink plot to make space for the legend: http://stackoverflow.com/a/4701285/323100
	#~ box = ax.get_position()
	#~ ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])		
	plt.legend(['Bridges (Evelyne) ','Volume (David)'],loc='best')
	plt.xlim([0,MaximumAcini])
	plt.xticks(arange(MaximumAcini))
	plt.ylim([0,None])
	plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_bridges_vs_volume.png',transparent=False)
if TikZTheData:
	tikz_save('plot_bridges_vs_volume.tikz')
	
#~ # http://is.gd/JWhkjn
#~ t = np.arange(10)
#~ plt.figure(num=None,figsize=(16,9))
#~ plt.plot(t, np.sin(t))
#~ print "Please click"
#~ x = plt.ginput(3)
#~ print "clicked",x[0][0]
#~ plt.show()
#~ 
#~ 
#~ sys.exit()


plt.show()
