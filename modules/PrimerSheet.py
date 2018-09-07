class PrimerSheet:
	pairIDs = []
	p5PrimersIDs = []
	p5Sequences = []
	p7PrimerIDs = []
	p7Sequences = []

	p5Pairs = []
	p7Pairs = []

	def ReadFile(self, name):
		file = open(name)
		file.readline()
		for line in file:
			cols = line.split("\t")
			cols[-1] = cols[-1].strip()
			#check num cols
			readPair = cols[0]
			pair = self.InsertPairID(cols[1])
			primer = self.InsertPrimer(readPair, cols[2])
			self.SetPrimerPair(readPair,primer,pair)
			self.InsertSequence(readPair, cols[3])
		file.close()

	def InsertPairID(self, pair):
		if pair not in self.pairIDs:
			self.pairIDs.append(pair)
			return len(self.pairIDs) - 1
		else:
			return self.pairIDs.index(pair)

	def InsertPrimer(self, readPair, primer):
		pairs = []
		if (readPair == "P5"):
			pairs = self.p5PrimersIDs;
		elif (readPair == "P7"):
			pairs = self.p7PrimerIDs;
		else:
			print("Invalid read pair")
			return;
		#if primer not in pairs:
		pairs.append(primer)
		return len(pairs) - 1
		#else:
		#	print("Duplicate primer ID %s" %primer)
		#	return pairs.index(primer)

	def SetPrimerPair(self, readPair, primerIndex, pairIndex):
		pairs = []
		if (readPair == "P5"):
			pairs = self.p5Pairs;
		elif (readPair == "P7"):
			pairs = self.p7Pairs;
		else:
			print("Invalid read pair")
			return;
		if (len(pairs) == primerIndex):
			pairs.append(pairIndex)
		else:
			print("Mismatched primer pair length")

	def InsertSequence(self, readPair, sequence):
		seq = []
		if (readPair == "P5"):
			seq = self.p5Sequences;
		elif (readPair == "P7"):
			seq = self.p7Sequences;
		else:
			print("Invalid read pair")
			return;
		seq.append(sequence)
		return len(seq) - 1

