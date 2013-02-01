# Script to find *all* instances of R108CXX-Samples, where XX = 04, 10, 21, 36, 60
# Made for having a list to go through, since TotalCommander was not useable for this.
# It's probably best to call this Python-script in the MS-DOS console with
# python FindAllSamples.py > Output.txt to have the output in a .txt-File

TeraStationsList = ['R:','S:','T:','V:']
RootDirectory = 'SLS'
SampleList = ['04','10','21','36','60']
		
import os
import fnmatch

matches = []
for TeraStation in TeraStationsList:
	SearchPath = os.path.join(TeraStation,os.sep,RootDirectory)
	print SearchPath
	for root, directories, files in os.walk(SearchPath):
		print root		
		for Sample in SampleList:
			for hit in fnmatch.filter(directories,'R108C' + Sample + '*'):
				print '---> found:',os.path.join(SearchPath,root,hit)
				matches.append(os.path.join(SearchPath,root,hit))
 
print 
print

print 'I found the following directories'
for hits in list(set(matches)):
	print hits
