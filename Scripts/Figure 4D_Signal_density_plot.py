#!/usr/bin/python
# required files : deeptools

import os, sys, time
from time import asctime

if len(sys.argv) != 6:
    print "Usage : ToolName.py <threads> <bigWigs> <BedFile> <StartEnd> <outputName>"
    print "MultiFile = comma divided"
    exit(1)
inFileThreads = sys.argv[1]
inFileBigWigs = sys.argv[2]
inFileBedFile = sys.argv[3]
inFileStartEnd = sys.argv[4]
outFileName = sys.argv[5]

# Function ----------------------------------------------------------
def timeNow ():
	localtime = time.asctime( time.localtime(time.time()) )
	print "###", localtime

def comma2space (inString): # Used multi bigWigs & StartEnd 
	tmp = inString.strip().split(",")

	spaceString = ""
	for splitNum in range(len(tmp)):
		if splitNum > 0:
			spaceString += " "
		spaceString += tmp[splitNum]
	return spaceString

def dotSplit (inString):
	tmp = inString.strip().split(",")
	for splitNum in range(len(tmp)):
		tmp[splitNum] = tmp[splitNum].split(".")[0]

	spaceString = ""
	for splitNum in range(len(tmp)):
		if splitNum > 0:
			spaceString += " "
		spaceString += tmp[splitNum]
	return spaceString
#-------------------------------------------------------------------
#Deeptools input arrangementSS
bigWigName = comma2space(inFileBigWigs)
StartEnd = inFileStartEnd.split(",")
labelName = dotSplit(inFileBigWigs)


### Deeptools : computeMatrix
print("# Deeptools : computeMatrix")
timeNow ()
print("### Code : computeMatrix reference-point --referencePoint center --skipZeros -p %s -S %s -R %s -a %s -b %s --samplesLabel %s -out %s.gz" % \
	(inFileThreads, bigWigName, inFileBedFile, StartEnd[0], StartEnd[1], labelName, outFileName))
os.system("computeMatrix reference-point --referencePoint center --skipZeros -p %s -S %s -R %s -a %s -b %s --samplesLabel %s -out %s.gz" % \
	(inFileThreads, bigWigName, inFileBedFile, StartEnd[0], StartEnd[1], labelName, outFileName))


## Deeptools : plotHeatmap
print("# Deeptools : plotHeatmap")
timeNow ()

os.system("plotHeatmap --refPointLabel center --heatmapHeight 7 --heatmapWidth 4 -m output/%s.gz -out output/Dpt_1out_comMarix_%s.pdf --yMin 0 -z Binding -x Distance" % \
	(outFileName, outFileName))
