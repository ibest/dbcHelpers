class SampleSheet:
	sampleIDs = []
	barcodeIDs = []
	projectIDs = []
	primerPairIDs = []

	sampleProject = []
	samplePrimers = []
	sampleBarcodes = []

	duplicateBarcodes = []

	def ReadFile(self, name):
		file = open(name)
		file.readline()
		for line in file:
			cols = line.split("\t")
			cols[-1] = cols[-1].strip()
			#check num cols
			sample = self.InsertSample(cols[0])
			primers = self.InsertPrimers(cols[1].split(','))
			self.SetSamplePrimers(sample,primers)
			barcode = self.InsertBarcode(cols[2])
			self.SetSampleBarcode(sample,barcode)
			project = self.InsertProject(cols[3])
			self.SetSampleProject(sample,project)
		file.close()

	def GetSampleProject(self, sampleIndex):
		return self.projectIDs[self.sampleProject[sampleIndex]]

	def GetSampleBarcode(self, sampleIndex):
		return self.barcodeIDs[self.sampleBarcodes[sampleIndex]]

	def GetSamplePrimers(self, sampleIndex):
		primers = []
		for primer in self.samplePrimers[sampleIndex]:
			primers.append(self.primerPairIDs[primer])
		return primers

	def InsertSample(self, sample):
		if sample not in self.sampleIDs:
			self.sampleIDs.append(sample)
			return len(self.sampleIDs) - 1
		else:
			print("Duplicate sample %s" %sample)
			return self.sampleIDs.index(sample)

	def InsertPrimers(self, primers):
		pid = []
		for primer in primers:
			if primer not in self.primerPairIDs:
				self.primerPairIDs.append(primer)
				pid.append(len(self.primerPairIDs) - 1)
			else:
				pid.append(self.primerPairIDs.index(primer))
		return pid

	def SetSamplePrimers(self, sampleIndex, primers):
		if (len(self.samplePrimers) == sampleIndex):
			self.samplePrimers.append(primers)
		else:
			print("Mismatched sample primers length")

	def InsertBarcode(self, barcode):
		if barcode not in self.barcodeIDs:
			self.barcodeIDs.append(barcode)
			return len(self.barcodeIDs) - 1
		else:
			print("Duplicate barcode %s" %barcode)
			i = self.barcodeIDs.index(barcode)
			self.duplicateBarcodes.append(i)
			return i

	def SetSampleBarcode(self, sampleIndex, barcodeIndex):
		if (len(self.sampleBarcodes) == sampleIndex):
			self.sampleBarcodes.append(barcodeIndex)
		else:
			print("Mismatched barcodes length %d" %sampleIndex)

	def InsertProject(self, project):
		if project not in self.projectIDs:
			self.projectIDs.append(project)
			return len(self.projectIDs) - 1
		else:
			return self.projectIDs.index(project)

	def SetSampleProject(self, sampleIndex, projectIndex):
		if (len(self.sampleProject) == sampleIndex):
			self.sampleProject.append(projectIndex)
		else:
			if (sampleIndex < len(self.sampleProject) and self.sampleProject[sampleIndex] != projectIndex):
				print("Duplicate sample %s with different project. Adding sample." %self.sampleIDs[sampleIndex])

	def WriteSampleSheet(self, fileName):
		file = open(fileName, 'w')
		file.write("SampleID\tPrimerPairID\tBarcodeID\tProjectID\n")
		for i in range(len(self.sampleIDs)):
			file.write("%s\t%s\t%s\t%s\n" %(self.sampleIDs[i],",".join(self.GetSamplePrimers(i)),self.GetSampleBarcode(i),self.GetSampleProject(i)))
		file.close()