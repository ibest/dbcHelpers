from modules.IdentifiedBarcodes import IdentifiedBarcodes
from modules.SampleSheet import SampleSheet
from modules.PrimerSheet import PrimerSheet
from modules.BarcodeLookup import BarcodeLookup
import modules.Validation as Validation
from datetime import datetime
import sys

def HammingDist(p1, p2, s1, s2):
	ss = s1 if len(s1) < len(s2) else s2
	sl = s1 if len(s1) >= len(s2) else s2
	ps = p1 if len(s1) < len(s2) else p2
	pl = p1 if len(s1) >= len(s2) else p2
	lenDiff = len(sl) - len(ss)
	minDist = float('inf');
	for i in range(4):
		diffs = i;
		if (lenDiff < i):
			diffs += i - lenDiff
		matches = []
		j = 0
		for k in range(i,len(sl)):
			if j >= len(ss):
				break
			if (Validation.CheckNucleotidesDifferent(ss[j],sl[k])):
				diffs += 1
			else:
				matches.append((j,k))
			j += 1
		minDist = min(diffs,minDist)
		diffs = i
		if (lenDiff < i):
			diffs += i - lenDiff
		matches = []
		j = 0
		for k in range(i,len(sl)):
			if k >= len(ss):
				break
			if (Validation.CheckNucleotidesDifferent(ss[k],sl[j])):
				diffs += 1
			else:
				matches.append((k,j))
			j += 1
		minDist = min(diffs,minDist)
	return minDist


sampleSheetName = None
primerSheetName = None
barcodeLookupName = None
identifiedBarcodesName = None
supressBarcodeWarning = False

if (len(sys.argv) > 0):
	count = 0
	while count < len(sys.argv):
		if sys.argv[count] == '-h' or sys.argv[count] == '-H' or sys.argv[count] == '--help':
			print("-B FILENAME, --barcodes_file FILENAME\tfile with barcodes")
			print("-P FILENAME, --primer_file FILENAME\tfile with primers")
			print("-S FILENAME, --sample_metadata FILENAME\tfile with sample metadata")
			print("-I FILENAME, --identified_barcodes FILENAME\tfile with identified barcodes")
			print("-E\tSuppress warnings for barcodes not in the sample sheet.")
			exit(0)
		if sys.argv[count] == "--sample_metadata" or sys.argv[count] == '-S':
			count += 1
			if (count < len(sys.argv)):
				sampleSheetName = sys.argv[count]
		elif sys.argv[count] == "--primer_file" or sys.argv[count] == '-P':
			count += 1
			if (count < len(sys.argv)):
				primerSheetName = sys.argv[count]
		elif sys.argv[count] == "--barcodes_file" or sys.argv[count] == '-B':
			count += 1
			if (count < len(sys.argv)):
				barcodeLookupName = sys.argv[count]
		elif sys.argv[count] == "--identified_barcodes" or sys.argv[count] == '-I':
			count += 1
			if (count < len(sys.argv)):
				identifiedBarcodesName = sys.argv[count]
		elif (sys.argv[count] == "-E"):
			supressBarcodeWarning = True
			count += 1
		count += 1

if sampleSheetName == None:
	sampleSheetName = input("Sample sheet name: ")
#if primerSheetName == None:
#	primerSheetName = input("Primer sheet name: ")
if barcodeLookupName == None:
	barcodeLookupName = input("Barcode lookup sheet name: ")
if identifiedBarcodesName == None:
	identifiedBarcodesName = input("Identified barcodes name: ")

identifiedBarcodes = IdentifiedBarcodes()
sampleSheet = SampleSheet()
barcodeLookup = BarcodeLookup()

identifiedBarcodes.ReadFile(identifiedBarcodesName)
sampleSheet.ReadFile(sampleSheetName)
barcodeLookup.ReadFile(barcodeLookupName)

spreadSheet = []
row = ["Project","Samples","Barcode"]
row[len(row):] = identifiedBarcodes.primers
spreadSheet.append(row)
for barcode in barcodeLookup.barcodeNames:
	row = [barcode]
	project = ""
	samples = []
	read = []
	if barcode in identifiedBarcodes.barcodes:
		barcodeIndex = identifiedBarcodes.barcodes.index(barcode)
		read = identifiedBarcodes.reads[barcodeIndex]
	else:
		read = [0]*len(identifiedBarcodes.reads[0])
	if barcode in sampleSheet.barcodeIDs:
		barcodeIndices = [i for i, x in enumerate(sampleSheet.barcodeIDs) if x == barcode]
		for barcodeIndex in barcodeIndices:
			sampleIndex = sampleSheet.sampleBarcodes.index(barcodeIndex)
			sample = sampleSheet.sampleIDs[sampleIndex]
			samples.append(sample)
			project = sampleSheet.projectIDs[sampleSheet.sampleProject[sampleIndex]]
	row.insert(0,samples)
	row.insert(0,project)
	row[len(row):] = read
	spreadSheet.append(row)
row = ["","","None"]
row[len(row):] = identifiedBarcodes.reads[-1]
spreadSheet.append(row)

for i in range(1,len(spreadSheet)-1):
	readSum = 0
	barcode = spreadSheet[i][2]
	primers = set()
	if barcode in sampleSheet.barcodeIDs:
		barcodeIndices = [i for i, x in enumerate(sampleSheet.barcodeIDs) if x == barcode]
		for barcodeIndex in barcodeIndices:
			sampleIndex = sampleSheet.sampleBarcodes.index(barcodeIndex)
			primers = primers | set(sampleSheet.GetSamplePrimers(sampleIndex))
	for j in range(3,len(spreadSheet[i])-1):
		primer = spreadSheet[0][j]
		readSum += int(spreadSheet[i][j])
		if (len(primers) > 0 and primer not in primers and int(spreadSheet[i][j]) > 10):
			print("Warning: Barcode %s and primer %s not paired but have %s reads." %(barcode,primer,spreadSheet[i][j]))
	readSum += int(spreadSheet[i][-1])
	if (readSum > 10):
		percent = int(spreadSheet[i][-1]) / readSum
		if (len(primers) == 0):
			if (supressBarcodeWarning == False):
				print("Warning: Barcode %s has %d reads, %d%% of which match no primers, but does not appear in the sample sheet. " %(barcode, readSum, percent * 100), end='')
				closest = None
				l1 = None
				l2 = None
				dist = float('inf')
				for i in range(len(barcodeLookup.barcodes)):
					if i == barcodeIndex or barcodeLookup.barcodeNames[i] not in sampleSheet.barcodeIDs:
						continue;
					check1 = HammingDist(barcode,barcodeLookup.barcodeNames[i],barcodeLookup.barcodes[barcodeIndex][0],barcodeLookup.barcodes[i][0])
					check2 = HammingDist(barcode,barcodeLookup.barcodeNames[i],barcodeLookup.barcodes[barcodeIndex][1],barcodeLookup.barcodes[i][1])
					#print("%d %d" %(check1,check2))
					thisDist = check1 + check2
					if thisDist < dist:
						dist = thisDist
						closest = barcodeLookup.barcodeNames[i]
						l1 = barcodeLookup.barcodes[i][0]
						l2 = barcodeLookup.barcodes[i][1]
				if (dist < 3):
					print("Possible barcode bleed: l1:%s, l2:%s. Possible source: %s l1:%s, l2:%s" %(barcodeLookup.barcodes[barcodeIndex][0],barcodeLookup.barcodes[barcodeIndex][1],closest,l1,l2))
				else:
					print("Possible barcode contamination: l1:%s, l2:%s. Possible source: %s l1:%s, l2:%s" %(barcodeLookup.barcodes[barcodeIndex][0],barcodeLookup.barcodes[barcodeIndex][1],closest,l1,l2))
		else:
			if percent >= 0.1:
				print("Warning: %d%% of %d read(s) for %s did not match a primer. " %(percent * 100,readSum,barcode), end='')
				closest = None
				l1 = None
				l2 = None
				dist = float('inf')
				for i in range(len(barcodeLookup.barcodes)):
					if i == barcodeIndex or barcodeLookup.barcodeNames[i] not in sampleSheet.barcodeIDs:
						continue;
					check1 = HammingDist(barcode,barcodeLookup.barcodeNames[i],barcodeLookup.barcodes[barcodeIndex][0],barcodeLookup.barcodes[i][0])
					check2 = HammingDist(barcode,barcodeLookup.barcodeNames[i],barcodeLookup.barcodes[barcodeIndex][1],barcodeLookup.barcodes[i][1])
					#print("%d %d" %(check1,check2))
					thisDist = check1 + check2
					if thisDist < dist:
						dist = thisDist
						closest = barcodeLookup.barcodeNames[i]
						l1 = barcodeLookup.barcodes[i][0]
						l2 = barcodeLookup.barcodes[i][1]
				if (dist < 3):
					print("Possible barcode bleed: l1:%s, l2:%s. Possible source: %s l1:%s, l2:%s" %(barcodeLookup.barcodes[barcodeIndex][0],barcodeLookup.barcodes[barcodeIndex][1],closest,l1,l2))
				else:
					print("Possible barcode contamination: l1:%s, l2:%s. Possible source: %s l1:%s, l2:%s" %(barcodeLookup.barcodes[barcodeIndex][0],barcodeLookup.barcodes[barcodeIndex][1],closest,l1,l2))
	#else:
	#	print("Barcode %s is fine." %barcode)
	
for project in sampleSheet.projectIDs:
	file = open("%s_%s.txt" %(datetime.now().strftime("%Y%m%d"), project),"w")
	projectRows = [i for i, row in enumerate(spreadSheet) if row[0] == project]
	projectIndex = sampleSheet.projectIDs.index(project)
	sampleIndices = [i for i, x in enumerate(sampleSheet.sampleProject) if x == projectIndex]
	primers = set()
	for i in sampleIndices:
		primers = primers | set(sampleSheet.GetSamplePrimers(i))
	header = [x for i, x in enumerate(spreadSheet[0]) if x in primers]
	header[:0] = spreadSheet[0][1:3]
	file.write("%s\n" %("\t".join(header)))
	for row in projectRows:
		data = [x for i, x in enumerate(spreadSheet[row]) if spreadSheet[0][i] in primers]
		data[:0] = spreadSheet[row][1:3]
		data[0] = ",".join(data[0])
		file.write("%s\n" %("\t".join(data)))
	file.close()

