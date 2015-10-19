#!/usr/bin/env python

#setup interpreter env
import os
os.chdir('/home/vince/code/PCR_error')

templateFile = 'template.txt'
maxCycle = 30

#examine template dimensions for array
lines = 0
with open(templateFile,'rU') as template:
	for seq in template:
		lines += 1

#read template sequences line by line
substrateArray = [i for i in range(lines)]
subArrayIndex = 0
with open(templateFile,'rU') as template:
	#build dictionary of unique substrate to be amplified; template data should be formatted as seq\tcount
	for seq in template:
		seq = str.split(seq)[0]
		substrateArray[subArrayIndex] = seq
		subArrayIndex += 1


