#!/usr/bin/env python
from multiprocessing import Process, Pool
import threading
import itertools
import argparse
import random
import time
'''
PCR Error Simulation
'''

RANDOM_BASE = ["A","C","G","T"]


def read_fasta_tuples(fasta):
	sequence_records = []
	with open(fasta, 'r') as fasta_handle:
		records_list = fasta_handle.read().split('>')
		for record in records_list[1::]:
			id, seq = record.rstrip('\n').split('\n')
			sequence_records.append((id,seq))
	return sequence_records

def make_output(records, outfile):
	with open(outfile, 'w') as handle:
		for record in records:
			id = ">{}".format(record[0])
			seq = record[1]
			handle.write("{}\n{}\n".format(id,seq))

if __name__ == '__main__':
	#  Argument Parser
	parser = argparse.ArgumentParser(description='PCR Error Simulation')

	# Input file
	parser.add_argument('-i','--input',dest='input', required=True, help='The input fasta')
	# Output file
	parser.add_argument('-o','--output',dest='output', required=True, help='The output file')
	# Cycles
	parser.add_argument('-c','--cycles',dest='cycles', default=1, help='Number of PCR cycles to simulate')
	# Error Rate
	parser.add_argument('-e','--error',dest='error', default=2272727, help='The inverse error rate')
	# Seed
	parser.add_argument('-s','--seed',dest='seed', default=time.time(), help='The seed for the PRNG')
	# Processes
	parser.add_argument('-t','--threads',dest='threads', default=8, help='The number of threads to start')

	# Parse arguments
	args = parser.parse_args()
	infile = args.input
	outfile = args.output
	cycles = int(args.cycles)
	error = int(args.error)
	seed = args.seed
	threads = int(args.threads)

	records = read_fasta_tuples(infile)

	def pcr_seq(seqrecord):
		pcr_product = []
		# Parent Seq
		#Copy original seq
		pcr_product.append(seqrecord)
		# PCR operation
		pcr_id = seqrecord[0]
		pcr_seq = list(seqrecord[1])
		# Introduce Error 
		mutations = 0
		mutation_list = []
		for index,base in enumerate(pcr_seq):
			if random.randint(1,error) == 1:
				#Shuffle
				newbase = random.choice(list(set(RANDOM_BASE)-set(base))[0])
				pcr_seq[index] = newbase
				mutations += 1
				mutation_list.append((index,base,newbase))
			else:
				pass
		pcr_id_mutations =  "{} M={} {}".format(pcr_id, mutations, mutation_list)
		pcr_seq_formatted = "".join(c for c in pcr_seq)
		pcr_product.append((pcr_id_mutations,pcr_seq_formatted))
		return pcr_product

	cycle = 1
	pcr_results = []
	while cycle <= cycles:
		pool = Pool(processes=threads)
		print "Cycle: {}".format(cycle)
		if cycle == 1:
			# Create Jobs
			map_results = pool.map(pcr_seq, records) # 2d array need to flatten
			pcr_results = list(itertools.chain(*map_results))
		else:
			# Create Jobs
			map_results = pool.map(pcr_seq, pcr_results)
			pcr_results = list(itertools.chain(*map_results))
		print len(pcr_results)
		cycle += 1

	#print pcr_results


	make_output(pcr_results,outfile)

