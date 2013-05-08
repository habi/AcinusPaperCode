# -*- coding: utf8 -*-

from pylab import *
import numpy as np
from matplotlib2tikz import save as tikz_save
import os
import glob
import csv
import xlrd # from http://www.lexicon.net/sjmachin/xlrd.htm, to read XLS-Files (like the ones containing the lung volumes)
import time
import fileinput

def tellme(s):
    print s
    plt.title(s,fontsize=16)
    plt.draw()

print "Hey ho, let's go: http://youtu.be/Q2Yb145eIIE !"
print '________________________________________________________________________________'
print


########################################################################
#                                                                      #
#	SETUP                                                              #
#                                                                      #
########################################################################

Drive = 'R:\SLS'

if os.path.exists(Drive) == False:
	print 'Cannot read ' + str(Drive) + '. Exiting!'
	exit()

verbose = False # set to 'False' to suppress some output.
RemoveOutliers = True
PlotTheData = True # True/False switches between showing and not showing the plots at the end
SaveFigures = True # save the plot to .png-Files
TikZTheData = True # save the data to .tikz-Files to import into LaTeX
TOMCATVoxelSize = 1.48
SliceNumber = 10
DisectorThickness = 5 # slices
ShrinkageFactor = 0.61 # Volume-Shrinkage-Factor = 61% with STD=5, calculated by Sébastien: Volume TOMCAT / Waterdisplacement
color=['c','r','m','b','y'] # colors to use, so that all the plots match
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
	
########################################################################
#                                                                      #
#	READING DATA                                                       #
#                                                                      #
########################################################################	

# Reading Volume Data of RUL from Stefans Data-File
print 'Reading values from p:\doc\#Dev\AcinarSize\DatenblattStefan.xls'
DatenblattStefan = xlrd.open_workbook('p:\doc\#Dev\AcinarSize\DatenblattStefan.xls').sheet_by_index(0)
RULVolume = [ DatenblattStefan.cell(70+int(Sample),9).value for Sample in range(len(Rat)) ]
BodyWeight = [ DatenblattStefan.cell(70+int(Sample),4).value for Sample in range(len(Rat)) ]
TotalLungVolume = [ DatenblattStefan.cell(70+int(Sample),5).value for Sample in range(len(Rat)) ]
AbsoluteParenchymalVolume = [ DatenblattStefan.cell(70+int(Sample),15).value for Sample in range(len(Rat)) ]
AbsoluteAirspaceVolume = [ DatenblattStefan.cell(70+int(Sample),19).value for Sample in range(len(Rat)) ]
AbsoluteAirspaceSurface = [ DatenblattStefan.cell(70+int(Sample),94).value for Sample in range(len(Rat)) ]

for Sample in range(len(Rat)):
	print WhichRat+Rat[Sample] + ': - RUL volume =',np.round(RULVolume[Sample],decimals=3),'cm^3'
	print (len(WhichRat+Rat[Sample])+1) * ' ','- Absolute parenchymal volume =',np.round(AbsoluteParenchymalVolume[Sample],decimals=3),'cm^2'
	print (len(WhichRat+Rat[Sample])+1) * ' ','- Absolute airspace volume =',np.round(AbsoluteAirspaceVolume[Sample],decimals=3),'cm^2'
	print (len(WhichRat+Rat[Sample])+1) * ' ','- Absolute airspace surface =',np.round(AbsoluteAirspaceSurface[Sample],decimals=3),'cm^2'

print
print 'DatenblattStefan.xls: - Mean RUL volume is',np.round(np.mean(RULVolume),decimals=3),'cm^3.'
print '                      - Mean absolute parenchymal volume is',np.round(np.mean(AbsoluteParenchymalVolume),decimals=3),'cm^2.'
print '                      - Mean absolute airspace volume is',np.round(np.mean(AbsoluteAirspaceVolume),decimals=3),'cm^2.'
print '                      - Mean absolute airspace surface is',np.round(np.mean(AbsoluteAirspaceSurface),decimals=3),'cm^2.'
print '________________________________________________________________________________'
print 
print "Setting Animals A and C to NaN, since we don't want to use their data for the biological table!"
print '________________________________________________________________________________'

BodyWeight[0] = np.nan
BodyWeight[2] = np.nan
TotalLungVolume[0] = np.nan
TotalLungVolume[2] = np.nan
AbsoluteParenchymalVolume[0] = np.nan
AbsoluteParenchymalVolume[2] = np.nan

# Reading Alveolar Number data from XLS-File from Lilian/Stefan
print 'Reading values from p:\doc\#Dev\AcinarSize\R108NumberofAlveoliDEF.xls'
DatenblattLilian = xlrd.open_workbook('p:\doc\#Dev\AcinarSize\R108NumberofAlveoliDEF.xls').sheet_by_index(1)
NumberOfAlveoliLilian = DatenblattLilian.cell(5,4).value
print 'Lilian/Stefan found',int(round(NumberOfAlveoliLilian)),'acini in a rat lung at day 60.'
NumberOfAciniRodriguez = 4023
print 'Rodriguez1987 (page 146) states a total of',NumberOfAciniRodriguez,'acini for the whole rat lungs.'
NumberOfAlveoliPerAcinusEstimated = NumberOfAlveoliLilian/NumberOfAciniRodriguez
print 'Grossly estimated, we thus have around (=Lilian/Rodriguez)',int(round(NumberOfAlveoliPerAcinusEstimated)),'alveoli per acinus.'
print '________________________________________________________________________________'
print

# See how many .csv-Files we actually have in the STEPanizer directories
SamplePath = {}
AssessedAcini = [0 for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	SamplePath[Sample] = os.path.join(Drive,Beamtime[Sample],WhichRat + Rat[Sample])
	if os.path.exists(SamplePath[Sample]):
		AssessedAcini[Sample] = len(glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerVolumeDir)))

# Generating a list of the .csv filenames
CSVFileVolume = [glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerVolumeDir,'*.xls')) for Sample in range(len(Rat))]	
CSVFileAlveoli = [glob.glob(os.path.join(SamplePath[Sample],'acin*',STEPanizerAlveoliDir,'*.xls')) for Sample in range(len(Rat))]	

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
	print 'check what is wrong and restart'
	sys.exit()
else:
	tmp = TotalAssessedAcini[0]
	TotalAssessedAcini = None
	TotalAssessedAcini = tmp
		
if MaximumAcini[0] != MaximumAcini[1]:
	print 'Maximum number of assessed acini for volume (David:',MaximumAcini[0],') and alveoli (Evelyne:',MaximumAcini[1],') does not match!'
	print 'check what is wrong and restart'
	sys.exit()
else:
	tmp = MaximumAcini[0]
	MaximumAcini = None
	MaximumAcini = tmp

# Read Volumes from .dcm-Files generated with MeVisLab
print 'Reading MeVisLab, volumes'
AcinarVolumeMeVisLab = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	for CurrentFile in sorted(glob.glob(os.path.join(SamplePath[Sample],'*.dcm'))):
		Acinus = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('.volume')])
		AcinarVolumeMeVisLab[Sample][Acinus] = float(CurrentFile[CurrentFile.find('volume')+len('volume'):CurrentFile.find('.pixelsize')])/1000 # normalize ul to cm^3: http://is.gd/XxU3Ei
		if verbose:
			print WhichRat + Rat[Sample] + ': Acinus',str(Acinus) + ', Volume =',AcinarVolumeMeVisLab[Sample][Acinus],'cm^3'

print 'Reading STEPanizer, volumes (David)'
# Read data from each STEPanizer .csv-file and calculate the desired values
AcinarVolumeSTEPanizer = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
SurfaceDensity = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
AbsoluteSurface = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
TotalTestPoints = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
NonParenchymalPoints= [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
ParenchymalPoints = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
ParenchymalVolume = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
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
		if verbose:
			print WhichRat + Rat[Sample] + ': Acinus', "%02d" % (Acinus),\
				'|',"%03d" % (TotalSlices),'Files',\
				'|',AllAreaTestPoints_Vol,'test points',\
				'|',"%04d" % (Interceptions),'Interceptions',\
				'|',"%03d" % (AcinusTestPoints),'Points in Acinus',\
				'|',"%03d" % (NonParenchymalPoints[Sample][Acinus]),'Nonparenchymal Points'
					
		# Volume = AcinusTestPoints * Area_Vol * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize		
		AcinarVolumeSTEPanizer[Sample][Acinus] = ((( AcinusTestPoints * Area_Vol * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize ) / ShrinkageFactor ) / 1e12 ) # scaling volume to cm^3: http://is.gd/wbZ81O	
						
		# Total of all points = AllAreaTestPoints_Vol (in file) * Total of slices
		TotalTestPoints[Sample][Acinus] = AllAreaTestPoints_Vol * TotalSlices
		
		# Parenchym-Volume = Parenchymal points * Area_Vol * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize. Parenchymal points are total points minus NParenchymalpoints (which are all points outside the sample and in non-parenchyma)
		ParenchymalPoints[Sample][Acinus] = TotalTestPoints[Sample][Acinus] - NonParenchymalPoints[Sample][Acinus]
		ParenchymalVolume[Sample][Acinus] = ParenchymalPoints[Sample][Acinus] * Area_Vol * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize
		
		# Surface Density = 2 * Interceptions / Length
		Length = ( LinePointLength_Vol * AcinusTestPoints ) / 1e4 # Linepointlength has to be calculated for the acinar reference space, /1e4 to scale from um to cm, http://is.gd/7wu1UD
		SurfaceDensity[Sample][Acinus] = 2 * Interceptions / Length
		
		# Absolute Surface = Surface density * acinar volume
		AbsoluteSurface[Sample][Acinus] = SurfaceDensity[Sample][Acinus] * AcinarVolumeSTEPanizer[Sample][Acinus]
	
print 'Reading STEPanizer, number of alveoli (Evelyne)'
# Read data from each STEPanizer .csv-file and calculate the desired values
AlveolarFraction = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
NumberOfAlveoli = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	for CurrentFile in sorted(CSVFileAlveoli[Sample][:]):
		Acinus = int(CurrentFile[CurrentFile.find('acinus')+len('acinus'):CurrentFile.find('acinus')+len('acinus')+2])
		DisectorThickness = float(CurrentFile[CurrentFile.find('Thickness-')+len('Thickness-'):CurrentFile.find('Thickness-')+len('Thickness-')+3]) # from filename, in um
		TotalSlices = len(glob.glob(os.path.join((os.path.join(SamplePath[Sample],'acinus'+ str("%02d" % (Acinus)),STEPanizerAlveoliDir)),'*.jpg')))
		FileData = csv.reader(open(CurrentFile,'rb'),dialect=csv.excel_tab)
		for line in FileData:
			if len(line) > 0:
				if line[0] == 'Num 1':
					Counts = int(line[3])
				if line[0] == 'Pixel size:':
					STEPanizerPixelSize_Alveoli = double(line[1])				
				if line[0] == 'a(p):':
					Area_Alveoli = double(line[1])*STEPanizerPixelSize_Alveoli**2
				if line[0] == 'l(p):':
					LinePointLength_Alveoli = double(line[1])*STEPanizerPixelSize_Alveoli
		
		# Counts are *all* counted bridges, (from a to b and from b to a). According 
		# to Stefan, we thus have to double the disector volume. This is then the 
		# volume density of the counts in said acinus.
		# Evelyne did count every second image, but this is really only relevant for
		# the volume of the acini. From her numbers we get the Counts per volume (=AlveolarFraction).
		# This is then multiplied by the volume of the acinus to get the number
		# of alveoli in each acinus. The volume is taken from Davids Cavaglieri estimation
		# above (AcinarVolumeSTEPanizer).
		AlveolarFraction[Sample][Acinus] = Counts / ( ( Area_Alveoli * ( DisectorThickness / ShrinkageFactor ) ) * 2 ) * 1e12 # Counts/cm^3
		# DisectorThickness = um, Area_Alveoli = um^2 -> 10^12 um^3 = 1 cm^3: http://is.gd/Cr6kUL
							
		NumberOfAlveoli[Sample][Acinus] = \
			AlveolarFraction[Sample][Acinus] * AcinarVolumeSTEPanizer[Sample][Acinus]
			
		# Hsiah2010 p. 407:
		# Counting the number of entrance rings in paired sections by the disector
		# technique allows estimation of total number of alveoli in the lung N(a,L) (112, 113).
		# N(a,L) is the product of the number of alveolar openings per unit parenchyma
		# volume (Sn/Vp) with the volume density of parenchyma per unit lung volume VV(p,L)
		# and the absolute lung volume:
		# N(a,L,) = (Sn/Vp) * VV(p,L) * V(L) (Formula 17)				

		# Give out counted/assessed data if desired
		if verbose:
			print WhichRat + Rat[Sample] + ': Acinus', "%02d" % (Acinus),\
				'| ' + str("%03d" % (TotalSlices)) + ' Files' +\
				' | ' + str(Counts) + ' counts' +\
				' | ' + str(int(AlveolarFraction[Sample][Acinus])) + ' AlveolarFraction' +\
				' | ' + str(int(NumberOfAlveoli[Sample][Acinus])) + ' alveoli'
				
########################################################################
#                                                                      #
#	CALCULATIONG DATA                                                  #
#                                                                      #
########################################################################
	
print 'Mean acinar volume MeVisLab'
AcinarVolumeMeanMeVisLab = [np.nan for Sample in range(len(Rat))]
AcinarVolumeSTDMeVisLab = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not scanned'
	else:
		AcinarVolumeMeanMeVisLab[Sample] = np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample]))
		AcinarVolumeSTDMeVisLab[Sample] = np.std(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample]))
		print WhichRat + Rat[Sample] + ':', AcinarVolumeMeanMeVisLab[Sample], 'cm^3'
print 'Mean acinar volume for all samples:', np.mean(np.ma.masked_invalid(AcinarVolumeMeanMeVisLab)), 'cm^3'
print 'Standard deviation of the mean acinar volume for all samples:', np.std(np.ma.masked_invalid(AcinarVolumeMeanMeVisLab))
print '________________________________________________________________________________'
print

print 'Mean acinar volume STEPanizer'
AcinarVolumeMeanSTEPanizer = [np.nan for Sample in range(len(Rat))]
AcinarVolumeSTDSTEPanizer = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not scanned'
	else:
		AcinarVolumeMeanSTEPanizer[Sample] = np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample]))
		AcinarVolumeSTDSTEPanizer[Sample] = np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample]))
		print WhichRat + Rat[Sample] + ':',AcinarVolumeMeanSTEPanizer[Sample],'cm^3'
print 'Mean acinar volume (mean of *all* acini):',np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer)),'cm^3, Standard deviation:',np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer))
print 'Mean acinar volume (mean of means of each sample):', np.mean(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer)),'cm^3, Standard deviation:',np.std(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer))
print

NumberOfOutliers = [ 0 for Sample in range(len(Rat))]
BiggerThan = 2
# Remove Outlier if > (mean +- BiggerThan*STD)
if RemoveOutliers:
	print 'Removing outliers if larger or smaller than mean +-',BiggerThan,'*STD'
	for Sample in range(len(Rat)):
		for Acinus in range(MaximumAcini):
			if not np.isnan(AcinarVolumeSTEPanizer[Sample][Acinus]):
				if not (
					AcinarVolumeMeanSTEPanizer[Sample] - ( BiggerThan * AcinarVolumeSTDSTEPanizer[Sample] )
					<=
					AcinarVolumeSTEPanizer[Sample][Acinus]
					<=
					AcinarVolumeMeanSTEPanizer[Sample] + ( BiggerThan * AcinarVolumeSTDSTEPanizer[Sample] )
					):
					print 'STEPanizer'
					print WhichRat + Rat[Sample] + ', Acinus',Acinus,\
						'volume:',np.round(AcinarVolumeSTEPanizer[Sample][Acinus],decimals=5),\
						'is larger than',\
						np.round(AcinarVolumeMeanSTEPanizer[Sample] + ( BiggerThan * AcinarVolumeSTDSTEPanizer[Sample] ),decimals=5),\
						'or smaller than',\
						np.round(AcinarVolumeMeanSTEPanizer[Sample] - ( BiggerThan * AcinarVolumeSTDSTEPanizer[Sample] ),decimals=5),\
						'thus setting to NaN'
					AcinarVolumeSTEPanizer[Sample][Acinus] = np.nan
					NumberOfOutliers[Sample] += 1
				if not (
					AcinarVolumeMeanMeVisLab[Sample] - ( BiggerThan * AcinarVolumeSTDMeVisLab[Sample] )
					<=
					AcinarVolumeMeVisLab[Sample][Acinus]
					<=
					AcinarVolumeMeanMeVisLab[Sample] + ( BiggerThan * AcinarVolumeSTDMeVisLab[Sample] )
					):
					print 'MeVisLab'
					print WhichRat + Rat[Sample] + ', Acinus',Acinus,\
						'volume:',np.round(AcinarVolumeMeVisLab[Sample][Acinus],decimals=5),\
						'is larger than',\
						np.round(AcinarVolumeMeanMeVisLab[Sample] + ( BiggerThan * AcinarVolumeSTDMeVisLab[Sample] ),decimals=5),\
						'or smaller than',\
						np.round(AcinarVolumeMeanMeVisLab[Sample] - ( BiggerThan * AcinarVolumeSTDMeVisLab[Sample] ),decimals=5),\
						'thus setting to NaN'
					AcinarVolumeMeVisLab[Sample][Acinus] = np.nan
					NumberOfOutliers[Sample] += 1
	for Sample in range(len(Rat)):
		if NumberOfOutliers[Sample] > 0:
			print WhichRat + Rat[Sample] + ': Removed',NumberOfOutliers[Sample],'Outlier(s) from the volume data'
	print '________________________________________________________________________________'		

	print 'After removing the outliers, we now have a'
	print 'Mean acinar volume MeVisLab'
	AcinarVolumeMeanMeVisLab = [np.nan for Sample in range(len(Rat))]
	for Sample in range(len(Rat)):
		if Beamtime[Sample] == '':
			print WhichRat + Rat[Sample] + ': not scanned'
		else:
			AcinarVolumeMeanMeVisLab[Sample] = np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample]))
			print WhichRat + Rat[Sample] + ':', AcinarVolumeMeanMeVisLab[Sample], 'cm^3'
	print 'Mean acinar volume for all samples:', np.mean(np.ma.masked_invalid(AcinarVolumeMeanMeVisLab)), 'cm^3'
	print 'Standard deviation of the mean acinar volume for all samples:', np.std(np.ma.masked_invalid(AcinarVolumeMeanMeVisLab))
	print '________________________________________________________________________________'
	print

print 'Mean acinar volume STEPanizer'
AcinarVolumeMeanSTEPanizer = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not scanned'
	else:
		AcinarVolumeMeanSTEPanizer[Sample] = np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample]))
		print WhichRat + Rat[Sample] + ':',AcinarVolumeMeanSTEPanizer[Sample],'cm^3'
print 'Mean acinar volume (mean of *all* acini):',np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer)),'cm^3, Standard deviation:',np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer))
print 'Mean acinar volume (mean of means of each sample):', np.mean(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer)),'cm^3, Standard deviation:',np.std(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer))
print

# Absolute airspace Volume / mean acinar volume = Number of Acini
print 'Number of acini (= absolute airspace volume from stefan / mean acinar volume)'
NumberOfAcini = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ' was not scanned'
	else:
		NumberOfAcini[Sample] = AbsoluteParenchymalVolume[Sample] / AcinarVolumeMeanSTEPanizer[Sample]
		print WhichRat + Rat[Sample],'contains',int(np.round(NumberOfAcini[Sample])),'acini'
print 'Rodriguez1987 (page 146) states a total of',NumberOfAciniRodriguez,' acini for the whole rat lungs.'
print 'We have a mean of',int(np.round(np.mean(np.ma.masked_invalid(NumberOfAcini)))),'acini for the whole rat lungs.'
print

for Sample in range(len(Rat)):
	if Beamtime[Sample]:
		print WhichRat + Rat[Sample] + ': -',np.round(np.mean(np.ma.masked_invalid(AlveolarFraction[Sample]))),'Alveoli per cm^3 (Evelyne	)'
		print (len(WhichRat + Rat[Sample])+1) * ' ','-',RULVolume[Sample],'cm^3 Volume (Stefan)'
		print (len(WhichRat + Rat[Sample])+1) * ' ','-',np.round(np.mean(np.ma.masked_invalid(AlveolarFraction[Sample])) * RULVolume[Sample] / 1e6 ),' * 10^6 estimated alveoli for the whole lung.'
print

print 'Mean number of alveoli per acinus'
MeanAlveolarNumber = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not scanned'
	else:
		MeanAlveolarNumber[Sample] = np.mean(np.ma.masked_invalid(NumberOfAlveoli[Sample]))
		print WhichRat + Rat[Sample] + ':', MeanAlveolarNumber[Sample]
print 'Mean number of alveoli (mean of all samples):', np.round(np.mean(np.ma.masked_invalid(NumberOfAlveoli))),'Standard deviation:',np.std(np.ma.masked_invalid(NumberOfAlveoli))
print 'Mean number of alveoli (mean of means of samples):', np.round(np.mean(np.ma.masked_invalid(MeanAlveolarNumber))),'Standard deviation:',np.std(np.ma.masked_invalid(MeanAlveolarNumber))
print 'From estimation of Alveoli per Lung (Lilian/Stefan) and # of acini (Rodriguez) we calculate an estimated',int(round(NumberOfAlveoliPerAcinusEstimated)),'alveoli per acinus.'
print 'Our number is thus',np.round(NumberOfAlveoliPerAcinusEstimated/np.mean(np.ma.masked_invalid(NumberOfAlveoli)),decimals=3),'times smaller.'
print

print 'Mean acinar surface'
MeanAcinarSurface = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ': not scanned'
	else:
		MeanAcinarSurface[Sample] = np.mean(np.ma.masked_invalid(AbsoluteSurface[Sample]))
		print WhichRat + Rat[Sample] + ':', MeanAcinarSurface[Sample], 'cm^2'
print 'Mean acinar surface for all samples:',np.round(np.mean(np.ma.masked_invalid(MeanAcinarSurface)),decimals=3), 'cm^2'
print 'Standard deviation of the mean acinar surface for all samples:', np.std(np.ma.masked_invalid(MeanAcinarSurface))
print

# Number of Acini * Mean acinar Surface = Mean airspace surface
print 'Diffusion surface (=number of acini * mean acinar surface)'
DiffusionSurface = [np.nan for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] == '':
		print WhichRat + Rat[Sample] + ' was not scanned'
	else:
		DiffusionSurface[Sample] = NumberOfAcini[Sample] * MeanAcinarSurface[Sample]
		print WhichRat + Rat[Sample],'contains',np.round(DiffusionSurface[Sample],3),'cm^2 of mean airspace surface'
print

print 'Stefan measured the absolute airspace surface with EM and got'
for Sample in range(len(Rat)):
	print WhichRat + Rat[Sample] +':',np.round(AbsoluteAirspaceSurface[Sample],decimals=3),'cm^2'
print
print 'Stefans mean absolute airspace surface is',np.round(np.mean(AbsoluteAirspaceSurface),decimals=3),'cm^2.'
print 'Our mean airspace surface is',np.round(np.mean(np.ma.masked_invalid(DiffusionSurface)),decimals=3),\
	'cm^2, approximately',np.round(np.mean(AbsoluteAirspaceSurface) / np.mean(np.ma.masked_invalid(DiffusionSurface)),decimals=3),'x smaller.'
print

# MeVisLab volume compared to STEPanizer volume (STEPanizer/MeVisLab)
STEPanizerMeVisLabDifference = [[np.nan for Acinus in range(MaximumAcini)] for Sample in range(len(Rat))]
for Sample in range(len(Rat)):
	if Beamtime[Sample] != '':
		for Acinus in range(MaximumAcini):
			if isnan(AcinarVolumeMeVisLab[Sample][Acinus]) == False:
				STEPanizerMeVisLabDifference[Sample][Acinus] = np.round((AcinarVolumeSTEPanizer[Sample][Acinus] / AcinarVolumeMeVisLab[Sample][Acinus]),decimals=3)
				if verbose:
					print WhichRat + Rat[Sample],'Acinus',str(Acinus) + ':',STEPanizerMeVisLabDifference[Sample][Acinus]
		if verbose:
			print '(a mean of',np.round(np.mean(np.ma.masked_invalid(STEPanizerMeVisLabDifference[Sample])),decimals=3),')'
	else:
		if verbose:
			print WhichRat + Rat[Sample],'was not scanned'

print '________________________________________________________________________________'
print
print 'MeVisLab to STEPanizer-comparison'
print 'Mean acinar volume'
for Sample in range(len(Rat)):
	if Beamtime[Sample] != '':
		print WhichRat + Rat[Sample] +': - mean acinar volume =',np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample])),'cm^3'
		print (len(WhichRat+Rat[Sample])+2) * ' ' + '- mean MeVisLab volume =',np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample])),'cm^3'
		print (len(WhichRat+Rat[Sample])+2) * ' ' + '- mean Differences = ',np.round(np.mean(np.ma.masked_invalid(STEPanizerMeVisLabDifference[Sample])),decimals=3),'times'
		
print 'The mean of each and every difference (STEPanizer-/MeVisLab-volume) is',np.round(np.mean(np.ma.masked_invalid(STEPanizerMeVisLabDifference)),decimals=3),'times bigger'


### Calculation for variant
NumberOfAciniVariant = [np.nan for Sample in range(len(Rat))]
TotalNumberOfAlveoli = [np.nan for Sample in range(len(Rat))]

TotalNumberOfAlveoli[1] = 39373856. # From DatenblattStefan.xls
TotalNumberOfAlveoli[3] = 56599446.
TotalNumberOfAlveoli[4] = 54728409.

for Sample in range(len(Rat)):
	NumberOfAciniVariant[Sample] = TotalNumberOfAlveoli[Sample] / np.mean(np.ma.masked_invalid(NumberOfAlveoli[Sample]))

print NumberOfAciniVariant
print np.std(np.ma.masked_invalid(NumberOfAciniVariant))
### Calculation for variant

# Write the data as variables directly to acinus.tex, so we don't have to copy-paste it all the time...
# The replacement stuff comes from http://is.gd/LIYXmb
for line in fileinput.FileInput('p:\\doc\\#Docs\\AcinusPaper\\acinus.tex',inplace=1):
	if '\\newcommand{\\shrinkagefactor}{' in line:
		print '\\newcommand{\\shrinkagefactor}{' + str(ShrinkageFactor) + '\\xspace} % Shrinkagefactor used for the calculation' 
	elif '\\newcommand{\\numberofacini}{' in line:
		print '\\newcommand{\\numberofacini}{' + str(TotalAssessedAcini) + '\\xspace}'
	elif '\\newcommand{\\numberofaciniB}{' in line:
		print '\\newcommand{\\numberofaciniB}{' + str(AssessedAcini[1]) + '\\xspace}'
	elif '\\newcommand{\\numberofaciniD}{' in line:
		print '\\newcommand{\\numberofaciniD}{' + str(AssessedAcini[3]) + '\\xspace}'
	elif '\\newcommand{\\numberofaciniE}{' in line:
		print '\\newcommand{\\numberofaciniE}{' + str(AssessedAcini[4]) + '\\xspace}'
	elif '\\newcommand{\\meantotalnumberofacini}{' in line:
		print '\\newcommand{\\meantotalnumberofacini}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAcini))))) + '\\xspace}'
	elif '\\newcommand{\\meantotalnumberofaciniSTD}{' in line:
		print '\\newcommand{\\meantotalnumberofaciniSTD}{' + str(int(np.round(np.std(np.ma.masked_invalid(NumberOfAcini))))) + '\\xspace} % add "ddof=1" to get the same STD as with "=STDEV()" in Excel'
	elif '\\newcommand{\\meantotalnumberofacinirounded}{' in line:
		print '\\newcommand{\\meantotalnumberofacinirounded}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAcini))/100)*100)) + '\\xspace}'
	elif '\\newcommand{\\totalnumberofaciniB}{' in line:
		print '\\newcommand{\\totalnumberofaciniB}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAcini[1]))))) + '\\xspace}'
	elif '\\newcommand{\\totalnumberofaciniD}{' in line:
		print '\\newcommand{\\totalnumberofaciniD}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAcini[3]))))) + '\\xspace}'
	elif '\\newcommand{\\totalnumberofaciniE}{' in line:
		print '\\newcommand{\\totalnumberofaciniE}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAcini[4]))))) + '\\xspace}'
	elif '\\newcommand{\\meantotalnumberofaciniVariant}{' in line:
		print '\\newcommand{\\meantotalnumberofaciniVariant}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAciniVariant))))) + '\\xspace}'
	elif '\\newcommand{\\meantotalnumberofaciniSTDVariant}{' in line:
		print '\\newcommand{\\meantotalnumberofaciniSTDVariant}{' + str(int(np.round(np.std(np.ma.masked_invalid(NumberOfAciniVariant))))) + '\\xspace} % add "ddof=1" to get the same STD as with "=STDEV()" in Excel'
	elif '\\newcommand{\\totalnumberofaciniBVariant}{' in line:
		print '\\newcommand{\\totalnumberofaciniBVariant}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAciniVariant[1]))))) + '\\xspace}'
	elif '\\newcommand{\\totalnumberofaciniDVariant}{' in line:
		print '\\newcommand{\\totalnumberofaciniDVariant}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAciniVariant[3]))))) + '\\xspace}'
	elif '\\newcommand{\\totalnumberofaciniEVariant}{' in line:
		print '\\newcommand{\\totalnumberofaciniEVariant}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAciniVariant[4]))))) + '\\xspace}'
	elif '\\newcommand{\\acinarvolumeB}{' in line:
		print '\\newcommand{\\acinarvolumeB}{' + str('%.3f' % (np.mean(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer[1])) * 1000 )) + '\\xspace} % mm^3 (changed from cm^3 to mm^3 in ReadVolumeSurfaceAndAlveaolarNumber.py) http://is.gd/Z3fUjF'
	elif '\\newcommand{\\acinarvolumeD}{' in line:
		print '\\newcommand{\\acinarvolumeD}{' + str('%.3f' % (np.mean(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer[3])) * 1000 )) + '\\xspace} % mm^3'
	elif '\\newcommand{\\acinarvolumeE}{' in line:
		print '\\newcommand{\\acinarvolumeE}{' + str('%.3f' % (np.mean(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer[4])) * 1000 )) + '\\xspace} % mm^3'
	elif '\\newcommand{\\meanacinarvolume}{' in line:
		print '\\newcommand{\\meanacinarvolume}{' + str('%.3f' % (np.mean(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer)) * 1000 )) + '} % mm^3, (mean acinar volume)'
	elif '\\newcommand{\\meanacinarvolumeSTD}{' in line:
		print '\\newcommand{\\meanacinarvolumeSTD}{' + str('%.3f' % (np.std(np.ma.masked_invalid(AcinarVolumeMeanSTEPanizer)) * 1000 )) + '} % (Standard deviation of acinar volumes), add "ddof=1" to get the same STD as with "=STDEV()" in Excel'
	elif '\\newcommand{\\meannumberofalveoli}{' in line:
		print '\\newcommand{\\meannumberofalveoli}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAlveoli))))) + '\\xspace} % (Mean number of alveoli per acinus)'
	elif '\\newcommand{\\meannumberofalveoliSTD}{' in line:
		print '\\newcommand{\\meannumberofalveoliSTD}{' + str(int(np.round(np.std(np.ma.masked_invalid(NumberOfAlveoli))))) + '\\xspace}'
	elif '\\newcommand{\\numberofalveoliB}{' in line:
		print '\\newcommand{\\numberofalveoliB}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAlveoli[1]))))) + '\\xspace}'
	elif '\\newcommand{\\numberofalveoliD}{' in line:
		print '\\newcommand{\\numberofalveoliD}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAlveoli[3]))))) + '\\xspace}'
	elif '\\newcommand{\\numberofalveoliE}{' in line:
		print '\\newcommand{\\numberofalveoliE}{' + str(int(np.round(np.mean(np.ma.masked_invalid(NumberOfAlveoli[4]))))) + '\\xspace}'
	elif '\\newcommand{\\difference}{' in line:
		print '\\newcommand{\\difference}{' + str(np.round(np.mean(np.ma.masked_invalid(STEPanizerMeVisLabDifference)),decimals=3)) + '\\xspace} % X times bigger (acinar volumes STEPanizer/MeVisLab-volumes)'
	elif '\\newcommand{\\meanacinarsurface}{' in line:
		print '\\newcommand{\\meanacinarsurface}{' + str(np.round(np.mean(np.ma.masked_invalid(MeanAcinarSurface)),decimals=3) * 100) + '} % mm^2 (changed from cm^2 to mm^2 in ReadVolumeSurfaceAndAlveaolarNumber.py) http://is.gd/99Sa9v'
	elif '\\newcommand{\\meanacinarsurfaceSTD}{' in line:
		print '\\newcommand{\\meanacinarsurfaceSTD}{' + str(np.round(np.std(np.ma.masked_invalid(MeanAcinarSurface)),decimals=3) * 100) + '}'
	elif '\\newcommand{\\acinarsurfaceB}{' in line:
		print '\\newcommand{\\acinarsurfaceB}{' + str(np.round(np.mean(np.ma.masked_invalid(MeanAcinarSurface[1])),decimals=3) * 100) + '\\xspace} % mm^2'
	elif '\\newcommand{\\acinarsurfaceD}{' in line:
		print '\\newcommand{\\acinarsurfaceD}{' + str(np.round(np.mean(np.ma.masked_invalid(MeanAcinarSurface[3])),decimals=3) * 100) + '\\xspace} % mm^2'
	elif '\\newcommand{\\acinarsurfaceE}{' in line:
		print '\\newcommand{\\acinarsurfaceE}{' + str(np.round(np.mean(np.ma.masked_invalid(MeanAcinarSurface[4])),decimals=3) * 100) + '\\xspace} % mm^2'
	elif '\\newcommand{\\meanairspacesurface}{' in line:
		print '\\newcommand{\\meanairspacesurface}{' + str(int(np.round(np.mean(np.ma.masked_invalid(DiffusionSurface))))) + '} % cm^2'
	elif '\\newcommand{\\airspacesurfaceB}{' in line:
		print '\\newcommand{\\airspacesurfaceB}{' + str(int(np.round(np.mean(np.ma.masked_invalid(DiffusionSurface[1]))))) + '\\xspace} % cm^2'
	elif '\\newcommand{\\airspacesurfaceD}{' in line:
		print '\\newcommand{\\airspacesurfaceD}{' + str(int(np.round(np.mean(np.ma.masked_invalid(DiffusionSurface[3]))))) + '\\xspace} % cm^2'
	elif '\\newcommand{\\airspacesurfaceE}{' in line:
		print '\\newcommand{\\airspacesurfaceE}{' + str(int(np.round(np.mean(np.ma.masked_invalid(DiffusionSurface[4]))))) + '\\xspace} % cm^2'
	elif '\\newcommand{\\airspacedifference}{' in line:
		print '\\newcommand{\\airspacedifference}{' + str(np.round(np.mean(AbsoluteAirspaceSurface) / np.mean(np.ma.masked_invalid(DiffusionSurface)),decimals=3)) + '\\xspace} % times'
	elif '\\newcommand{\\numberofoutliers}{' in line:
		print '\\newcommand{\\numberofoutliers}{' + str(sum(NumberOfOutliers)) + '\\xspace} % Number of Outliers removed'
	elif '\\newcommand{\\biggerthan}{' in line:
		print '\\newcommand{\\biggerthan}{' + str(BiggerThan) + '\\xspace} % Outliers bigger/smaller than mean +- BiggerThan * STD have been removed'
	elif '\\newcommand{\\bodyweightB}{' in line:
		print '\\newcommand{\\bodyweightB}{' + str(int(round(np.mean(np.ma.masked_invalid(BodyWeight[1]))))) + '\\xspace} % g'
	elif '\\newcommand{\\bodyweightD}{' in line:
		print '\\newcommand{\\bodyweightD}{' + str(int(round(np.mean(np.ma.masked_invalid(BodyWeight[3]))))) + '\\xspace} % g'
	elif '\\newcommand{\\bodyweightE}{' in line:
		print '\\newcommand{\\bodyweightE}{' + str(int(round(np.mean(np.ma.masked_invalid(BodyWeight[4]))))) + '\\xspace} % g'
	elif '\\newcommand{\\meanbodyweight}{' in line:
		print '\\newcommand{\\meanbodyweight}{' + str(int(round(np.mean(np.ma.masked_invalid(BodyWeight))))) + '\\xspace} % g'
	elif '\\newcommand{\\meanbodyweightSTD}{' in line:
		print '\\newcommand{\\meanbodyweightSTD}{' + str('%.3f' % (np.std(np.ma.masked_invalid(BodyWeight)))) + '\\xspace}'
	elif '\\newcommand{\\totallungvolumeB}{' in line:
		print '\\newcommand{\\totallungvolumeB}{' + str(int(np.mean(np.ma.masked_invalid(TotalLungVolume[1]) * 1000 ))) + '\\xspace} % mm^3'
	elif '\\newcommand{\\totallungvolumeD}{' in line:
		print '\\newcommand{\\totallungvolumeD}{' + str(int(np.mean(np.ma.masked_invalid(TotalLungVolume[3]) * 1000 ))) + '\\xspace} % mm^3'
	elif '\\newcommand{\\totallungvolumeE}{' in line:
		print '\\newcommand{\\totallungvolumeE}{' + str(int(np.mean(np.ma.masked_invalid(TotalLungVolume[4]) * 1000 ))) + '\\xspace} % mm^3'
	elif '\\newcommand{\\meantotallungvolume}{' in line:
		print '\\newcommand{\\meantotallungvolume}{' + str(int(round(np.mean(np.ma.masked_invalid(TotalLungVolume) * 1000 )))) + '} % mm^3'
	elif '\\newcommand{\\meantotallungvolumeSTD}{' in line:
		print '\\newcommand{\\meantotallungvolumeSTD}{' + str('%.3f' % (np.std(np.ma.masked_invalid(TotalLungVolume) * 1000 ))) + '\\xspace}'
	elif '\\newcommand{\\parenchymalvolumeB}{' in line:
		print '\\newcommand{\\parenchymalvolumeB}{' + str(int(round(np.mean(np.ma.masked_invalid(AbsoluteParenchymalVolume[1]) * 1000 )))) + '\\xspace} % mm^3'
	elif '\\newcommand{\\parenchymalvolumeD}{' in line:
		print '\\newcommand{\\parenchymalvolumeD}{' + str(int(round(np.mean(np.ma.masked_invalid(AbsoluteParenchymalVolume[3]) * 1000 )))) + '\\xspace} % mm^3'
	elif '\\newcommand{\\parenchymalvolumeE}{' in line:
		print '\\newcommand{\\parenchymalvolumeE}{' + str(int(round(np.mean(np.ma.masked_invalid(AbsoluteParenchymalVolume[4]) * 1000 )))) + '\\xspace} % mm^3'
	elif '\\newcommand{\\meanparenchymalvolume}{' in line:
		print '\\newcommand{\\meanparenchymalvolume}{' + str(int(round(np.mean(np.ma.masked_invalid(AbsoluteParenchymalVolume) * 1000 )))) + '\\xspace} % mm^3'
	elif '\\newcommand{\\meanparenchymalvolumeSTD}{' in line:
		print '\\newcommand{\\meanparenchymalvolumeSTD}{' + str('%.3f' % (np.std(np.ma.masked_invalid(AbsoluteParenchymalVolume) * 1000 ))) + '\\xspace}'
	else:
		print line, # the ',' at the end prevents a newline

if ShrinkageFactor != 1:
	print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	print '!                                                        !'
	print '!   All is calculated with a Shrinkagefactor of',100*ShrinkageFactor,'%   !'
	print '!                                                        !'	
	print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print

if PlotTheData == False:
	exit()

if TikZTheData:
	print "Using 'plt.plot' instead of 'plt.scatter', since 'tikz_save('file.tikz')' doesn't work otherwise"
	print "just add 'only marks' to the TikZ-code (in the \\addplot-command)"

PlotMeVisLabVolume = []
PlotSTEPAnizerVolume = []
for Sample in range(len(Rat)):
	for Acinus in range(MaximumAcini):
		if isnan(AcinarVolumeMeVisLab[Sample][Acinus]) == False:
			PlotMeVisLabVolume.append(AcinarVolumeMeVisLab[Sample][Acinus])
			PlotSTEPAnizerVolume.append(AcinarVolumeSTEPanizer[Sample][Acinus])
			
PlotNormalizedMeVisLabVolume = [float(i)/max(PlotMeVisLabVolume) for i in PlotMeVisLabVolume]
PlotNormalizedSTEPanizerVolume = [float(i)/max(PlotSTEPAnizerVolume) for i in PlotSTEPAnizerVolume]

# Plot volumes for all samples
plt.figure(num=None,figsize=(16,9))
for Sample in range(len(Rat)):
	if Beamtime[Sample]:
		# MeVisLab
		plt.subplot(2,len(Rat),Sample+1)
		plt.scatter(range(MaximumAcini),AcinarVolumeMeVisLab[Sample],color=color[Sample])
		# Plot mean MeVisLab-Volume
		plt.axhline(np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample])),color='c',linestyle='dashed',linewidth=2,label='MeVisLab mean +- 3xSTD') # no x-range necessary
		plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample]))+(3*np.std(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample])))),color='c',linestyle='dotted',alpha=0.5)
		plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample]))-(3*np.std(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample])))),color='c',linestyle='dotted',alpha=0.5)
		plt.title(WhichRat+Rat[Sample] + '\nMeVisLab (' + str(AssessedAcini[Sample]) +')')			
		if Sample == 0:
			plt.ylabel('Volume [cm^3]')
		#~ plt.xticks(arange(TotalAssessedAcini))	
		#~ plt.xlabel('Acinus')
		plt.xlim([0,None])
		# STEPanizer
		plt.subplot(2,len(Rat),Sample+len(Rat)+1)
		plt.scatter(range(MaximumAcini),AcinarVolumeSTEPanizer[Sample],color=color[Sample])
		# Plot mean STEPanizer-Volume
		plt.axhline(np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample])),color='c',linestyle='dashed',linewidth=2,label='STEPanizer mean +- 3xSTD') # no x-range necessary
		plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample]))+(3*np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample])))),color='c',linestyle='dotted',alpha=0.5)
		plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample]))-(3*np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample])))),color='c',linestyle='dotted',alpha=0.5)
		plt.title('STEPanizer (' + str(AssessedAcini[Sample]) +')')
		if Sample == 0:
			plt.ylabel('Volume [cm^3]')	
		#~ plt.xticks(arange(TotalAssessedAcini))
		plt.xlabel('Acinus')
		plt.xlim([0,None])
if SaveFigures:
	plt.savefig('plot_mevislab_stepanizer_volumes.png')
if TikZTheData:
	tikz_save('plot_mevislab_stepanizer_volumes.tikz')

sys.exit()
# Plot volumes for single samples
for Sample in range(len(Rat)):
	if Beamtime[Sample]:
		# MeVisLab
		plt.figure(num=None,figsize=(16,9))
		plt.plot(range(MaximumAcini),AcinarVolumeMeVisLab[Sample],color=color[Sample])
		# Plot mean MeVisLab-Volume
		plt.axhline(np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample])),color='c',linestyle='dashed',linewidth=2,label='MeVisLab mean & 3xSTD') # no x-range necessary
		plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample]))+(3*np.std(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample])))),color='c',linestyle='dashed',alpha=0.5)
		plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample]))-(3*np.std(np.ma.masked_invalid(AcinarVolumeMeVisLab[Sample])))),color='c',linestyle='dashed',alpha=0.5)
		plt.title(WhichRat+Rat[Sample] + ' MeVisLab (' + str(AssessedAcini[Sample]) +')')			
		if Sample == 0:
			plt.ylabel('Volume [mm^3]')
		#~ plt.xticks(arange(TotalAssessedAcini))	
		#~ plt.xlabel('Acinus')
		plt.xlim([0,None])
		if SaveFigures:
			plt.savefig('plot_mevisvolume_' + str(WhichRat+Rat[Sample]) + '.png')
		if TikZTheData:
			tikz_save('plot_mevisvolume_' + str(WhichRat+Rat[Sample]) + '.tikz')
		# STEPanizer
		plt.figure(num=None,figsize=(16,9))
		plt.plot(range(MaximumAcini),AcinarVolumeSTEPanizer[Sample],color=color[Sample])
		# Plot mean STEPanizer-Volume
		plt.axhline(np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample])),color='c',linestyle='dotted',linewidth=2,label='STEPanizer mean & 3xSTD') # no x-range necessary
		plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample]))+(3*np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample])))),color='c',linestyle='dotted',alpha=0.5)
		plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample]))-(3*np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer[Sample])))),color='c',linestyle='dotted',alpha=0.5)
		plt.title('STEPanizer (' + str(AssessedAcini[Sample]) +')')
		if Sample == 0:
			plt.ylabel('Volume [mm^3]')	
		#~ plt.xticks(arange(TotalAssessedAcini))
		plt.xlabel('Acinus')
		plt.xlim([0,None])
		if SaveFigures:
			plt.savefig('plot_stepanizervolume_' + str(WhichRat+Rat[Sample]) + '.png')
		if TikZTheData:
			tikz_save('plot_stepanizervolume_' + str(WhichRat+Rat[Sample]) + '.tikz',encoding='utf8')
			
plt.show()

# plot both MeVisLab- and STEPanizer-Volumes in one plot
plt.figure(num=None,figsize=(16,9))
for Sample in range(len(Rat)):
	plt.plot(range(MaximumAcini*Sample,
		MaximumAcini*(Sample+1)),
		AcinarVolumeSTEPanizer[Sample],
		color[Sample],
		marker='o',
		label=str(WhichRat+Rat[Sample]) + ': STEPanizer')	
	plt.plot(range(MaximumAcini*Sample,MaximumAcini*(Sample+1)),
		AcinarVolumeMeVisLab[Sample],
		color[Sample]+'--',
		marker='*',
		label=str(WhichRat+Rat[Sample]) + ': MeVisLab')
# Plot mean STEPanizer-Volume
plt.axhline(np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer)),color='c',linestyle='dotted',linewidth=2,label='STEPanizer mean +- 3xSTD') # no x-range necessary
plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer))+(3*np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer)))),color='c',linestyle='dotted',alpha=0.5)
plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeSTEPanizer))-(3*np.std(np.ma.masked_invalid(AcinarVolumeSTEPanizer)))),color='c',linestyle='dotted',alpha=0.5)
	
# Plot mean MeVisLab-Volume
plt.axhline(np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab)),color='c',linestyle='dashed',linewidth=2,label='MeVisLab mean +- 3xSTD') # no x-range necessary
plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab))+(3*np.std(np.ma.masked_invalid(AcinarVolumeMeVisLab)))),color='c',linestyle='dashed',alpha=0.5)
plt.axhline((np.mean(np.ma.masked_invalid(AcinarVolumeMeVisLab))-(3*np.std(np.ma.masked_invalid(AcinarVolumeMeVisLab)))),color='c',linestyle='dashed',alpha=0.5)
plt.title('Volumes')
plt.legend(loc='best')
if SaveFigures:
	plt.savefig('plot_stepanizervolumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_stepanizervolumes.tikz')
plt.draw()
plt.show()

import win32api
win32api.MessageBox(0, 'Select one Acinus, I will then show you a slice of the selected acinus afterwards', 'Acinus Selection', 0x00001000)
tellme('Select acinus to look at a slice\nStop with middle mouse button')

AcinusToLookAtCoordinates = plt.ginput(1,timeout=0)
AcinusToLookAt = int(round(AcinusToLookAtCoordinates[0][0]))

plt.plot(AcinusToLookAtCoordinates[0][0],AcinusToLookAtCoordinates[0][1],'go',markersize=20.0,alpha=0.5)
plt.draw()

for Sample in range(len(Rat)):
	if (AcinusToLookAt >= MaximumAcini*Sample) and (AcinusToLookAt < MaximumAcini*(Sample+1)):
		print 'The selected Acinus',AcinusToLookAt,'corresponds to acinus #',\
			AcinusToLookAt-MaximumAcini*Sample,'of the sorted list of acini of sample',WhichRat+Rat[Sample],\
			'which is',\
			os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample])[\
			os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample]).find('acinus'):\
			os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample]).find('acinus')+len('acinus')+2]
		#~ print CSVFileVolume[Sample]
		print 'This acinus can be found in the folder ',os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample]),\
			' which contains',\
			len(glob.glob(os.path.join(os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample]),'*.jpg'))),\
			'.jpg files'
		print 'We are opening',\
			glob.glob(os.path.join(os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample]),'*.jpg'))[int(len(glob.glob(os.path.join(os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample]),'*.jpg')))/2)],\
			'to look at.'
		os.startfile(glob.glob(os.path.join(os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample]),'*.jpg'))[int(len(glob.glob(os.path.join(os.path.dirname(CSVFileVolume[Sample][AcinusToLookAt-MaximumAcini*Sample]),'*.jpg')))/2)])

print '________________________________________________________________________________'
print

# Plot the interesting stuff
# Plot MeVisLab- against STEPanizer-volumes
plt.figure(num=None,figsize=(16,9))
for Sample in range(len(Rat)):
	plt.subplot(len(Rat),1,Sample+1)
	if Beamtime[Sample] != '':
		if TikZTheData:
			plt.plot(range(MaximumAcini),
				AcinarVolumeMeVisLab[Sample],
				marker='o',
				c='g')
			plt.plot(range(MaximumAcini),
				AcinarVolumeSTEPanizer[Sample],
				marker='o',
				c='b')
		else:
			plt.scatter(range(MaximumAcini),AcinarVolumeMeVisLab[Sample],c='g')
			plt.scatter(range(MaximumAcini),AcinarVolumeSTEPanizer[Sample],c='b')
			plt.legend(['MeVisLab','STEPanizer'],loc='best')
	if Beamtime[Sample] == '':
		plt.title('Volumes of acini of ' + WhichRat + Rat[Sample] + ' were not assessed')
	else:
		plt.title('Volumes of acini of ' + Beamtime[Sample] + '\\' + WhichRat+Rat[Sample])
	plt.xlim([0,MaximumAcini])
	#~ plt.ylim([0,None])
	plt.xticks(arange(MaximumAcini))
plt.tight_layout()
if SaveFigures:
	plt.savefig('plot_mevis_vs_stepanizer_volumes.png',transparent=False)
if TikZTheData:
	tikz_save('plot_mevis_vs_stepanizer_volumes.tikz')

#~ plt.figure(num=None,figsize=(16,9))
#~ for Sample in range(len(Rat)):
	#~ plt.subplot(len(Rat),1,Sample+1)
	#~ if TikZTheData:
		#~ plt.plot(range(MaximumAcini),
			#~ STEPanizerMeVisLabDifference[Sample],
			#~ marker='o',
			#~ c=color[Sample])
	#~ else:
		#~ plt.scatter(range(MaximumAcini),STEPanizerMeVisLabDifference[Sample],c=color[Sample])		
	#~ plt.xlim([0,MaximumAcini])
	#~ plt.ylim([0,None])
	#~ if Beamtime[Sample] == '':
		#~ plt.title(WhichRat + Rat[Sample])
	#~ else:
		#~ plt.title(Beamtime[Sample] + '\\' + WhichRat + Rat[Sample] + ': STEPanizer/MeVisLab-Difference')
	#~ plt.xticks(arange(MaximumAcini))
#~ plt.tight_layout()
#~ if SaveFigures:
	#~ plt.savefig('plot_mevis_vs_stepanizer_volumeratio.png',transparent=False)
#~ if TikZTheData:
	#~ tikz_save('plot_mevis_vs_stepanizer_volumeratio.tikz')
#~ 
#~ # Boxplot of Volumes
#~ plt.figure(num=None,figsize=(16,9))
#~ for Sample in range(len(Rat)):
	#~ plt.xlim([0,2])
	#~ plt.subplot(3,len(Rat),Sample+1)#,axisbg=color[Sample])
	#~ if Beamtime[Sample] != '':
		#~ plt.boxplot([x for x in AcinarVolumeSTEPanizer[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	#~ plt.title(str(WhichRat) + str(Rat[Sample])+'\nRUL-Volume='+str(RULVolume[Sample])+' cm^3')
	#~ plt.ylabel('Acinar Volume')
	#~ plt.subplot(3,len(Rat),Sample+1+len(Rat))#,axisbg=color[Sample])
	#~ if Beamtime[Sample] != '':
		#~ plt.boxplot([x for x in SurfaceDensity[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	#~ plt.ylabel('Surface Density')
	#~ plt.subplot(3,len(Rat),Sample+1+2*len(Rat))#,axisbg=color[Sample])
	#~ if Beamtime[Sample] != '':
		#~ plt.boxplot([x for x in AbsoluteSurface[Sample] if not math.isnan(x)],1) # http://is.gd/Ywxpiz remove all np.Nan for boxplotting
	#~ plt.ylabel('Absolute Surface')
#~ plt.tight_layout()
#~ if SaveFigures:
	#~ plt.savefig('plot_boxplot_of_volumes.png',transparent=False)
#~ if TikZTheData:
	#~ tikz_save('plot_boxplot_of_volumes.tikz')
#~ 
#~ # Plotting Volumes, Surface Density and Absolute Surfaces
#~ plt.figure(num=None,figsize=(16,9))
#~ for Sample in range(len(Rat)):
	#~ plt.subplot(311)
	#~ plt.xlim([0,MaximumAcini])
	#~ if Beamtime[Sample] != '':
		#~ if TikZTheData:
			#~ plt.plot(range(MaximumAcini),
				#~ AcinarVolumeSTEPanizer[Sample],
				#~ marker='o',
				#~ c=color[Sample])
		#~ else:
			#~ plt.scatter(range(MaximumAcini),AcinarVolumeSTEPanizer[Sample],c=color[Sample])
	#~ plt.title('Volume (AcinusTestPoints * Area * STEPanizerPixelSize_Vol * SliceNumber * TOMCATVoxelSize)')
	#~ plt.xlim([0,MaximumAcini])
	#~ plt.ylim([0,max(max(AcinarVolumeSTEPanizer))])
	#~ plt.subplot(312)
	#~ plt.xlim([0,MaximumAcini])
	#~ if Beamtime[Sample] != '':
		#~ if TikZTheData:
			#~ plt.plot(range(MaximumAcini),
				#~ SurfaceDensity[Sample],
				#~ marker='o',
				#~ c=color[Sample])
		#~ else:
			#~ plt.scatter(range(MaximumAcini),SurfaceDensity[Sample],c=color[Sample])
	#~ plt.title('Surface Density (2 * Int / Length)')
	#~ plt.subplot(313)
	#~ plt.xlim([0,MaximumAcini])
	#~ plt.ylim([0,None])
	#~ if Beamtime[Sample] != '':
		#~ if TikZTheData:
			#~ plt.plot(range(MaximumAcini),
				#~ AbsoluteSurface[Sample],
				#~ marker='o',
				#~ c=color[Sample])
		#~ else:
			#~ plt.scatter(range(MaximumAcini),AbsoluteSurface[Sample],c=color[Sample])
	#~ plt.title('Absolute Surface (SurfaceDensity * AcinarVolumeSTEPanizer)')
	#~ plt.legend([Beamtime[1] + "\\" + WhichRat + Rat[1],Beamtime[3] + '\\' + WhichRat + Rat[3],Beamtime[4] + '\\' + WhichRat + Rat[4]],loc='best')#loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3)
	#~ plt.xlim([0,MaximumAcini])
	#~ plt.ylim([0,None])
	#~ plt.tight_layout()
#~ if SaveFigures:
	#~ plt.savefig('plot_acinarvolume_surfacedensity_absolutesurface.png',transparent=False)
#~ if TikZTheData:
	#~ tikz_save('plot_acinarvolume_surfacedensity_absolutesurface.tikz')
#~ 
#~ # Plot AlveolarFraction vs. Volume (Evelyne vs. David)
#~ plt.figure(num=None,figsize=(16,9))
#~ plt.subplots_adjust(hspace=1)
#~ for Sample in range(len(Rat)):
	#~ ax = plt.subplot(len(Rat),1,Sample+1)
	#~ plt.scatter(range(MaximumAcini),AlveolarFraction[Sample],c='r')
	#~ plt.scatter(range(MaximumAcini),AcinarVolumeSTEPanizer[Sample],c='b')
	#~ plt.title('AlveolarFraction ' + WhichRat + Rat[Sample] +': ' +\
		#~ str(np.sum(np.ma.masked_invalid(AlveolarFraction[Sample]))) +\
		#~ ' (total)')
	#~ # Shink plot to make space for the legend: http://stackoverflow.com/a/4701285/323100
	#~ box = ax.get_position()
	#~ ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])		
	#~ plt.legend(['AlveolarFraction (Evelyne) ','Volume (David)'],loc='best')
	#~ plt.xlim([0,MaximumAcini])
	#~ plt.xticks(arange(MaximumAcini))
	#~ plt.ylim([0,None])
	#~ plt.tight_layout()
#~ if SaveFigures:
	#~ plt.savefig('plot_AlveolarFraction_vs_volume.png',transparent=False)
#~ if TikZTheData:
	#~ tikz_save('plot_AlveolarFraction_vs_volume.tikz')
	#~ 
#~ # Plot AlveolarFraction vs. Volume (Evelyne)
#~ plt.figure(num=None,figsize=(16,9))
#~ for Sample in range(len(Rat)):
	#~ plt.subplot(2,len(Rat),Sample)
	#~ plt.plot(AlveolarFraction[Sample],marker='o')
	#~ plt.title('AlveolarFraction ' + str(WhichRat) + str(Rat[Sample]))
	#~ plt.subplot(2,len(Rat),Sample+len(Rat))
	#~ plt.plot(NumberOfAlveoli[Sample],marker='o')
	#~ plt.title('# of Alveoli ' + str(WhichRat) + str(Rat[Sample]))
#~ if SaveFigures:
	#~ plt.savefig('plot_AlveolarFraction_and_number_of_acini.png',transparent=False)
#~ if TikZTheData:
	#~ tikz_save('plot_AlveolarFraction_and_number_of_acini.tikz')
	#~ 

# http://is.gd/JWhkjn
# -*- noplot -*-
# t = np.arange(10)
# plt.figure(num=None,figsize=(16,9))
# plt.plot(t, np.sin(t))
# print "Please click"
# x = plt.ginput(3)
# print "clicked",x[0][0]
# plt.show()

plt.show()
