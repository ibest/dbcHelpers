from modules.IdentifiedBarcodes import IdentifiedBarcodes
from modules.SampleSheet import SampleSheet
from modules.PrimerSheet import PrimerSheet
from modules.BarcodeLookup import BarcodeLookup
from datetime import datetime
import sys

sampleSheetName = None
primerSheetName = None
barcodeLookupName = None
identifiedBarcodesName = None

if (len(sys.argv) > 0):
	count = 0
	while count < len(sys.argv):
		if sys.argv[count] == '-h' or sys.argv[count] == '-H' or sys.argv[count] == '--help':
			print("-B FILENAME, --barcodes_file FILENAME\tfile with barcodes")
			print("-P FILENAME, --primer_file FILENAME\tfile with primers")
			print("-S FILENAME, --sample_metadata FILENAME\tfile with sample metadata")
			print("-I FILENAME, --identified_barcodes FILENAME\tfile with identified barcodes")
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
			print("Warning: Barcode %s has %d reads, %d%% of which match no primers, but does not appear in the sample sheet." %(barcode, readSum, percent * 100))
		else:
			if percent >= 0.1:
				print("Warning: %d%% of %d read(s) for %s did not match a primer." %(percent * 100,readSum,barcode))
	
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

