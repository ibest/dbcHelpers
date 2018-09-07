class BarcodeLookup:
	barcodeNames = []
	barcodes = []

	def ReadFile(self, name):
		file = open(name)
		file.readline()
		for line in file:
			cols = line.split("\t")
			#check num cols
			cols[-1] = cols[-1].strip();
			self.barcodeNames.append(cols[0])
			self.barcodes.append(cols[1:])
		file.close()