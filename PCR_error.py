#!/usr/bin/env python

#setup interpreter env
import os
import random
os.chdir('/home/vince/code/PCR_error/PCR_error_primary/50000_primary_singleton_filtered')

templateFile = 'filtered.fna'
outFna = '3_cycle_primary_singleton_filtered_HFtaq.fna'
ID = 'KIM8002'
maxCycle = 3
#taq polymerase error rate = 2.28E-5
invErrRate = int(1/4.4E-7)

#examine template dimensions for array
lines = 0
with open(templateFile,'rU') as template:
	for seq in template:
		if '>' in seq:
			continue
		else:
			lines += 1

#create array according to above dimensions
substrateArray = [i for i in range(lines)]
subArrayIndex = 0
#read template sequences line by line
with open(templateFile,'rU') as template:
	#fill array with template
	for seq in template:
		if '>' in seq:
			continue
		else:
			seq = str.split(seq)[0]
			substrateArray[subArrayIndex] = [seq,0]
			subArrayIndex += 1

#cycle
n = 1
#starting template size
n0 = lines
#count errors
errCount = 0

#begin cycle; run until specified max cycle
while n < maxCycle:
	cycleCoeff = 2**n
	#productArray size = double substrateArray length for given cycle
	productArray = [i for i in range(n0*cycleCoeff)]
	#"a" indexes substrate carryover
	for a in range(len(substrateArray)):
		# substrateArray[0] = seq, substrateArray[1] = # of sustained mutations
		productArray[a] = substrateArray[a]
	#switch a for j; "j" indexes synthesized product
	for j in range(len(substrateArray)):
		synthStr = ['',substrateArray[j][1]]
		for k in range(len(substrateArray[j][0])):
			mut = random.randrange(0,invErrRate)
			base = substrateArray[j][0][k]
			otherBases = set(['A','G','T','C']) - set(base)
			if mut == 5:
				#replace w/ 'X' or any other marker if needed
				#random.sample returns list w/ 1 element
				synthStr[0] += random.sample(otherBases,1)[0]
				synthStr[1] += 1
				errCount += 1
			else:
				synthStr[0] += substrateArray[j][0][k]
		l = j + len(substrateArray)
		productArray[l] = synthStr
	# #advance cycle count
	n += 1
	#convert products to substrate
	substrateArray = productArray
	#end cycle
#write modified filtered.fna
with open(outFna,'w') as out:
	for z in range(0,len(productArray)):
		out.write('>' + ID + '_' + str(z + 1) + '\n')
		out.write(productArray[z][0] + '\n')

