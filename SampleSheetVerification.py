from modules.SampleSheet import SampleSheet
from modules.PrimerSheet import PrimerSheet
from modules.BarcodeLookup import BarcodeLookup
import modules.Validation as Validation
import sys
from datetime import datetime

sampleSheetName = None
primerSheetName = None
barcodeLookupName = None
newSampleSheetName = None

if (len(sys.argv) > 0):
	count = 0
	while count < len(sys.argv):
		if sys.argv[count] == '-h' or sys.argv[count] == '-H' or sys.argv[count] == '--help':
			print("-B FILENAME, --barcodes_file FILENAME\tfile with barcodes")
			print("-P FILENAME, --primer_file FILENAME\tfile with primers")
			print("-S FILENAME, --sample_metadata FILENAME\tfile with sample metadata")
			print("-n FILENAME\tCustom output samplesheet name. Default is {YearMonthDay}SampleSheet.txt")
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
		elif sys.argv[count] == "--output_name" or sys.argv[count] == '-N':
			count += 1
			if (count < len(sys.argv)):
				newSampleSheetName = sys.argv[count]
		count += 1

if sampleSheetName == None:
	sampleSheetName = input("Sample sheet name: ")
if primerSheetName == None:
	primerSheetName = input("Primer sheet name: ")
if barcodeLookupName == None:
	barcodeLookupName = input("Barcode lookup sheet name: ")
if newSampleSheetName == None:
	newSampleSheetName = "%s%s" %(datetime.now().strftime("%Y%m%d"), "SampleSheet.txt")

sampleSheet = SampleSheet()
primerSheet = PrimerSheet()
barcodeLookup = BarcodeLookup()

sampleSheet.ReadFile(sampleSheetName)
primerSheet.ReadFile(primerSheetName)
barcodeLookup.ReadFile(barcodeLookupName)


print("Sanitizing names...")
sampleSheet.sampleIDs = Validation.SanitizeNames(sampleSheet.sampleIDs, "SID")
sampleSheet.barcodeIDs = Validation.SanitizeNames(sampleSheet.barcodeIDs, "BID")
sampleSheet.projectIDs = Validation.SanitizeNames(sampleSheet.projectIDs, "PID")
sampleSheet.primerPairIDs = Validation.SanitizeNames(sampleSheet.primerPairIDs, "PPID")
primerSheet.pairIDs = Validation.SanitizeNames(primerSheet.pairIDs, "PPID")
primerSheet.p5PrimersIDs = Validation.SanitizeNames(primerSheet.p5PrimersIDs, "P5ID")
primerSheet.p7PrimerIDs = Validation.SanitizeNames(primerSheet.p7PrimerIDs, "P7ID")
barcodeLookup.barcodeNames = Validation.SanitizeNames(barcodeLookup.barcodeNames, "BID")

print("Checking for primer collisions...")
primerCollisions = Validation.CheckPrimerCollision(primerSheet.p5PrimersIDs, primerSheet.p5Sequences)
primerCollisions[len(primerCollisions):] = Validation.CheckPrimerCollision(primerSheet.p7PrimerIDs, primerSheet.p7Sequences)

print("Checking for primers not represented in the sample sheet...")
primersRepresented = Validation.CheckPrimersRepresented(sampleSheet.primerPairIDs, primerSheet.pairIDs)

print("Checking that all p5 and p7 primers have a matching pair...")
primerPairsValid = Validation.CheckPrimerPairs(primerSheet.p5Pairs, primerSheet.p7Pairs, primerSheet.pairIDs)

print("Checking for duplicate barcodes...")
duplicateBarcodes = Validation.CheckDuplicateBarcodes(sampleSheet.duplicateBarcodes, sampleSheet.barcodeIDs, sampleSheet.primerPairIDs, sampleSheet.sampleBarcodes, sampleSheet.samplePrimers, primerCollisions)

if len(primerCollisions) > 0 or not primersRepresented or not primerPairsValid or not duplicateBarcodes:
	response = input("Print sample sheet? (Yes/No): ")
	if (len(response) > 0 and (response[0] != 'y' and response[0] != 'Y')):
		exit(1)

print("Writing sample sheet...")
sampleSheet.WriteSampleSheet(newSampleSheetName)