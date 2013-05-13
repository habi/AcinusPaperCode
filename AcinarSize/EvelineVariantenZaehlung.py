# -*- coding: utf8 -*-

import os
import glob
import csv
import numpy as np

print 80 * '_'

Drive = 'R:\SLS'

if os.path.exists(Drive) == False:
    print 'Cannot read ' + str(Drive) + '. Exiting!'
    exit()

Experiment = 'R108C'
Day = 4
Animal = 'A'
Beamtime = '2010a'

Disector = [3, 5, 7]

Acinus = [25, 28, 30]

# Find all .xls-files in the folders for each acinus ($Acinus) and each
# disector thickness ($Disector, while making shure that we have two numbers
# after the comma).
# This gives us a list in a list, which we then covert to a simple list with
# "''.join(biglist)".
ExportFile = [[''.join(glob.glob(os.path.join(Drive, Beamtime, 'mrg',
                                              str(Experiment +\
                                                  '%02d' % Day +\
                                                  Animal + 't-mrg'),
                                              'acinus' + str(WhichAcinus),
                                              'voxelsize1.48-every10slice-' +\
                                              'DisectorThickness-' +\
                                              str('%0.2f' %
                                                  (thickness * 1.48)) +\
                                              'um-or' + str(thickness) +\
                                              'slices', '*.xls')))
    for thickness in Disector]
    for WhichAcinus in Acinus]

Counts = [[0 for i in range(3)] for k in range(3)]
PixelSize = [[0 for i in range(3)] for k in range(3)]
CountingArea = [[0 for i in range(3)] for k in range(3)]
DisectorThickness = [[0 for i in range(3)] for k in range(3)]

# Go into the files and read the relevant data
for i in range(3):
    for k in range(3):
        print 'Reading data from', ExportFile[i][k]

# Read and show the data
print 80 * '_'
for i in range(3):
    print 'For Acinus', Acinus[i]
    for k in range(3):
        FileData = csv.reader(open(ExportFile[i][k], 'rb'),
                              dialect=csv.excel_tab)
        for line in FileData:
            if len(line) > 0:
                if any('SUM' in item for item in line):
                    Counts[i][k] = int(line[3])
                if any('Pixel size:' in item for item in line):
                    PixelSize[i][k] = float(line[1])
                if any('a(p):' in item for item in line):
                    CountingArea[i][k] = int(line[1])
        DisectorThickness[i][k] = float('%0.2f' % (Disector[k] * 1.48))
        print '  Disector distance:', Disector[k], 'slices (' +\
            str(DisectorThickness[i][k]) + ' um):'
        print '    -', Counts[i][k], 'counts at'
        print '    -', PixelSize[i][k], 'um pixel size'
        print '    - and thus a countig area of', CountingArea[i][k],\
            'px or', int(round(CountingArea[i][k] * PixelSize[i][k] ** 2)),\
            'um^2'

print 80 * '_'

AlveolarFraction = [[0 for i in range(3)] for k in range(3)]

# Calculation from ReadVolumeSurfaceAndAlveaolarNumber.py, line 245ff
print 'We can thus calculate the counts per volume (alveolar fracton) for',\
    'each acinus, which is:'
for i in range(3):
    for k in range(3):
        AlveolarFraction[i][k] = Counts[i][k] / ((CountingArea[i][k] *
                                                  DisectorThickness[i][k]) *
                                                  2) * 1e12
            
for i in range(3):
    print 'For acinus', Acinus[i], 'we found a mean of',\
        str(int(round(np.mean(np.mean(AlveolarFraction[i]))))),'counts/cm^3'
    for k in range(3):
        print '  At a disector distance of', Disector[k],\
            'slices'
        print '    -', int(round(AlveolarFraction[i][k])), 'Counts/cm^3.'
        print '    - This is the mean plus', \
            round(((np.mean(AlveolarFraction[i]) - AlveolarFraction[i][k])
                    / np.mean(AlveolarFraction[i]))*100,2), '%'

print 80 * '_'
