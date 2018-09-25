class IdentifiedBarcodes:
	barcodes = []
	primers = []
	reads = []

	def ReadFile(self, name):
		file = open(name)
		self.InsertPrimers(file.readline())
		for line in file:
			self.ReadLine(line)
		file.close()

	def InsertPrimers(self, line):
		cols = line.split("\t")
		cols[-1] = cols[-1].strip()
		cols = cols[1:]
		for col in cols:
			self.primers.append(col)

	def ReadLine(self, line):
		cols = line.split("\t")
		cols[-1] = cols[-1].strip()
		self.barcodes.append(cols[0])
		self.reads.append(cols[1:])
