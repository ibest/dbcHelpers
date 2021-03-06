import string
import re
import sys

def HammingDist(p1, p2, s1, s2):
	ss = s1 if len(s1) < len(s2) else s2
	sl = s1 if len(s1) >= len(s2) else s2
	ps = p1 if len(s1) < len(s2) else p2
	pl = p1 if len(s1) >= len(s2) else p2
	lenDiff = len(sl) - len(ss)
	for i in range(5):
		diffs = i;
		if (lenDiff < i):
			diffs += i - lenDiff
		matches = []
		j = 0
		for k in range(i,len(sl)):
			if j >= len(ss):
				break
			if (CheckNucleotidesDifferent(ss[j],sl[k])):
				diffs += 1
			else:
				matches.append((j,k))
			j += 1
		if (diffs <= 4):
			PrintPrimerCollision(ps,pl,ss,sl,matches)
			return True
		diffs = i
		if (lenDiff < i):
			diffs += i - lenDiff
		matches = []
		j = 0
		for k in range(i,len(sl)):
			if k >= len(ss):
				break
			if (CheckNucleotidesDifferent(ss[k],sl[j])):
				diffs += 1
			else:
				matches.append((k,j))
			j += 1
		if (diffs <= 4):
			PrintPrimerCollision(ps,pl,ss,sl,matches)
			return True
	return False

def PrintPrimerCollision(ps, pl, ss, sl, matches):
	print("Error: %s collides with %s" %(pl,ps))
	ssWidth = len(ss)
	slWidth = len(sl)
	diff = matches[0][0] - matches[0][1]
	if diff < 0:
		ssWidth += abs(diff)
	elif diff > 0:
		slWidth += abs(diff)
	matchIndices = []
	for match in matches:
		if diff == 0:
			matchIndices.append(match[0])
		elif diff > 0:
			matchIndices.append(match[0])
		else:
			matchIndices.append(match[1])
	print('{:>{width}}'.format(sl, width=slWidth))
	for i in range(slWidth):
		if i in matchIndices:
			print("|", end='')
		else:
			print(" ", end='')
	print()
	print('{:>{width}}'.format(ss, width=ssWidth))
	print("\n")

def CheckNucleotidesDifferent(c1, c2):
	c1Set = CreateNucleotideSet(c1)
	c2Set = CreateNucleotideSet(c2)
	if (len(c1Set & c2Set) == 0):
		return True
	return False
	

def CreateNucleotideSet(c):
	if (c == 'A'):
		return set(['A'])
	if (c == 'G'):
		return set(['G'])
	if (c == 'C'):
		return set(['C'])
	if (c == 'T'):
		return set(['T'])
	if (c == 'U'):
		return set(['U'])
	if (c == 'R'):
		return set(['A','G'])
	if (c == 'Y'):
		return set(['C','T'])
	if (c == 'N'):
		return set(['A','G','C','T'])
	if (c == 'W'):
		return set(['A','T'])
	if (c == 'S'):
		return set(['G','C'])
	if (c == 'M'):
		return set(['A','C'])
	if (c == 'K'):
		return set(['G','T'])
	if (c == 'B'):
		return set(['G','C','T'])
	if (c == 'H'):
		return set(['A','C','T'])
	if (c == 'D'):
		return set(['A','G','T'])
	if (c == 'V'):
		return set(['A','G','C'])

def CheckPrimerCollision(primers, sequences):
	collisions = []
	for i in range(len(primers)):
		for j in range(i,len(primers)):
			if (primers[i] == primers[j]):
				continue
			if HammingDist(primers[i],primers[j],sequences[i],sequences[j]):
				collisions.append((i,j))
	return collisions

def CheckPrimersRepresented(ssPrimers, psPrimers):
	notRepresented = set(psPrimers) - set(ssPrimers)
	if (len(notRepresented) > 0):
		print("Warning: primers in primer sheet but not represented in sample sheet: %s\n" %(', '.join(notRepresented)))
		return False
	return True

def CheckPrimerPairs(p5Pairs, p7Pairs, pairIDs):
	notSymmetric = set(p5Pairs) ^ set(p7Pairs)
	if (len(notSymmetric) > 0):
		pairNames = []
		for pair in notSymmetric:
			pairNames.append(pairIDs[pair]);
		print("Warning: primer pair IDs not used in p5 and p7 primer groups: %s\n" %', '.join(pairNames))
		return False
	return True

def CheckDuplicateBarcodes(duplicates, barcodeIDs, primerIDs, sampleBarcodes, samplePrimers, primerCollisions):
	error = False
	for duplicate in duplicates:
		primers = []
		for i in range(len(sampleBarcodes)):
			if (sampleBarcodes[i] == duplicate):
				primers.append(samplePrimers[i])
		if (len(primers) > 0):
			for p in range(0,len(primers)):
				collisions = set(primers[p]);
				c = []
				for primerCollision in primerCollisions:
					if primerCollision[0] in collisions and primerCollisions[1] in collisions:
						c.append(primerIDs[primerCollision[0]])
						c.append(primerIDs[primerCollision[1]])
				if (len(c) > 0):
					print("Error: Barcode %s used with colliding primers: %s" %(barcodeIDs[duplicate], ", ",join(c)))
			collisions = set(primers[0])
			for p in range(1,len(primers)):
				collisions = collisions & set(primers[p])
			if (len(collisions) == 0):
				print("Warning: Barcode %s used multiple times, but with different primers." %barcodeIDs[duplicate])
				error = True
			else:
				print("Error: Barcode %s used multiple time with the same primers: " %barcodeIDs[duplicate], end = '')
				collisionNames = []
				for collision in collisions:
					collisionNames.append(primerIDs[collision])
				print(", ".join(collisionNames))
				error = True
	return error


def SanitizeNames(names, field):
	newNames = []
	for name in names:
		newName = name
		if (field == "SID" and not newName[0].isalpha()):
			newName = "s%s" %newName
		newName = re.sub("[^_.A-Za-z0-9]+",".",newName)
		if newName != name:
			print("Warning: %s: %s renamed %s\n" %(field,name,newName))
		newNames.append(newName)
	return newNames