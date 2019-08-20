#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.7'

'''
Dependence: pysam (module).
Sorted BAM with index.
There must be barcode at the end of QNAME column in BAM/SAM files and seperated by ":".
e.g. NB501707:7:HNFCWAFXX:1:21101:5445:18780:AGGAAAGAAGCC
Correct barcode with 2 bases mis-match or 1 base indels.
Indels in read can cause quality decreas but do not influence sequence result.
And merge same barcode reads together.
	1. previous cluter output
'''

BARCODE_LEN = 12  # barcode length is 12bp
HELF_BC_LEN = int(BARCODE_LEN/2)  # BARCODE_LEN/2
MERGE_CUTOFF = 5  # same barcode reads merge cutoff
EDIT_DIST_CUTOFF = 2  # edit distance cutoff, 2 mismatch, 1 indel
SLIDE_CUTOFF = 3  # read coordinate slide distance cutoff
BARCODE_SEPARATOR = ':'  # barcode in qname and the sepatator

import pysam
import gzip

class BcRead(object):
	'''light weight class to store both reads and barcode.'''
	__slots__ = '__read', '__barcode'
	def __init__(self, read, barcode):
		self.__read = read
		self.__barcode = barcode
	def get_read(self):
		return self.__read
	def get_barcode(self):
		return self.__barcode
	def set_barcode(self, barcode):
		self.__barcode = barcode

def main(input_file, output_file, region_chrom, region_start, region_end):
	input_frag, result_BC, input_BC = {}, {}, {}  # input ,result and total fragment
	cor_read = 0  # corrected reads and fragments
	stat = {}  # state of same barcode loaded reads number for read 1
	import os.path
	if os.path.splitext(input_file)[1] == '.bam':  # BAM file
		bamfile = pysam.AlignmentFile(input_file, "rb")
	else:
		raise IOError('Please choose a sorted BAM/SAM file.')
	read_PL = PositionalList()  # a positional list for store reads information
	cluster_chrom, cluster_start, cluster_end = None, 0, 0
	with pysam.AlignmentFile('%s.filtered.bam' % output_file, 'wb', header=bamfile.header) as no_pass_bam, pysam.AlignmentFile('%s.M%d.bam' % (output_file, MERGE_CUTOFF), "wb", header=bamfile.header) as pass_bam, gzip.open('%s.BC_correct.txt.gz' % output_file, 'wt') as BC_correct_txt:
		for read in bamfile.fetch(contig=region_chrom, start=region_start, stop=region_end):
			if read.is_unmapped or read.cigarstring is None:  # * cigar
				continue
			qname = read.query_name
			qname = qname.split(' ')[0]  # in case of ' ' seperation in qname
			barcode = qname.split(BARCODE_SEPARATOR)[-1]  # barcodo in qname: NB501707:7:HNFCWAFXX:1:21101:5445:18780:AGGAAAGAAGCC
			chrom = read.reference_name
			input_frag[qname] = input_frag.get(qname, 0) + 1  # same paired reads in same qname
			input_BC[barcode] = input_BC.get(barcode, 0) + 1  # total input barcode number
			ref_pos_start, ref_pos_end = read.reference_start, read.reference_end
			# R1/R2 determine
			if read.is_read1:
				read_1_2 = 'r1'
			elif read.is_read2:
				read_1_2 = 'r2'
			else:
				read_1_2 = 'r0'
			barcode += ':%s' % read_1_2  # add R1/R2 information in barcode
			if cluster_chrom is not None:  # not first read in BAM
				if chrom != cluster_chrom or abs(ref_pos_start - cluster_start) >= SLIDE_CUTOFF:  # at least one end is same
					output_BAM(read_PL, no_pass_bam, pass_bam, chrom, cluster_chrom, cluster_start, cluster_end, stat, result_BC)  # clear the list of reads
					if len(read_PL) >= MERGE_CUTOFF:  # if not directly clear the list
						cor_read += barcode_correct(read_PL, BC_correct_txt, pass_bam, no_pass_bam, stat, result_BC)  # find barcode family and correct error barcode
					cluster_chrom, cluster_start, cluster_end = chrom, ref_pos_start, ref_pos_end  # new cluster
			else:  # first read in BAM
				cluster_chrom, cluster_start, cluster_end = chrom, ref_pos_start, ref_pos_end
			read_PL.add_last(BcRead(read, barcode))  # add read information to linked list
		cor_read += barcode_correct(read_PL, BC_correct_txt, pass_bam, no_pass_bam, stat, result_BC)  # fininsh last reads
		output_BAM(read_PL, no_pass_bam, pass_bam, chrom, 'end', cluster_start, cluster_end, stat, result_BC)  # clear the list of reads
	bamfile.close()
	with open('%s.cluster.log' % output_file, 'w') as log_txt:
		result_member, result_fam_read, filtered_fam_read = 0, 0, 0
		for k,v in stat.items():  # k is member, v is count number
			if k >= MERGE_CUTOFF:
				result_fam_read += v  # family(reads) of >= M5 reads
				result_member += k*v  # member of >= M5 reads
			else:
				filtered_fam_read += v  # family(reads) of < M5 reads
		input_hist = {}  # input family member hist
		for k,v in input_BC.items():  # k is member, v is count number
			input_hist[v] = input_hist.get(v, 0) + 1

		input_member = sum(input_frag.values())  # reads
		input_frag = len(input_frag)  # fragment = single reads + pair reads (qname)
		input_fam = len(input_BC)  # family = multi fragment group (barcode)
		result_read = sum(result_BC.values())  #
		result_fam = len(result_BC)  # barcode of >= M5 reads 
		'''
		input_fam, result_fam family are sum of barcode types.
		input_member, result_member and filtered_fam_read member are sum of reads.
		result_read and result_fam_read are famliy in reads. 
		'''
		log_txt.write('========== Raw Data ==========\n')
		log_txt.write('Input reads    :\t%d\n' % (input_member))
		log_txt.write('Input fragments:\t%d\n' % (input_frag))
		log_txt.write('Input family   :\t%d\n' % (input_fam))
		log_txt.write('Mean members per fam(reads):\t%.2f\n' % (input_member/input_fam))
		log_txt.write('========== Process  ==========\n')
		log_txt.write('Filtered reads :\t{:d}\t{:.2%}\n'.format(input_member-result_member, (input_member-result_member)/input_member))
		log_txt.write('Corrected reads:\t{:d}\t{:.2%}\n'.format(cor_read, cor_read/input_member))
		log_txt.write('Merged reads   :\t{:d}\t{:.2%}\n'.format(result_member, result_member/input_member))
		log_txt.write('==========  Result  ==========\n')
		log_txt.write('Filtered M1-M4 family(reads) :\t%d\n' % (filtered_fam_read))
		log_txt.write('Total M5 family(reads)       :\t%d\n' % (result_fam_read))
		log_txt.write('Total M5 family(fragments)   :\t%d\n' % (result_fam))
		log_txt.write('Mean reads per fam(reads)    :\t%.2f\n' % (result_member/result_fam_read))
		log_txt.write('Mean reads per fam(fragments):\t%.2f\n' % (result_member/result_fam))
		log_txt.write('======== Input & output =======\n')
		log_txt.write('#Member\tInput_family\tOutput_reads\n')
		for k,v in sorted(input_hist.items()):  # k is member, v is count number
			out_read = stat[k] if k in stat else 0
			log_txt.write('%d\t%d\t%d\n' % (k, v, out_read))
	return True

def barcode_correct(read_PL, BC_correct_txt, pass_bam, no_pass_bam, stat, result_BC):
	'''correct barcode in same coordinate cluster.'''
	barcode_seq = {}  # key is barcode and value is corrected barcode
	barcode_fre = {}  # key is barcode and value is barcode frequency
	cor_read = 0  # corrected reads/fragments number
	# store barcode seq and barcode frequency in dict
	cursor = read_PL.first()
	while cursor is not None:
		this_barcode = cursor.get_element().get_barcode()
		barcode_seq[this_barcode] = this_barcode  # key is barcode and value will correct barcode
		barcode_fre[this_barcode] = barcode_fre.get(this_barcode, 0) + 1  #key is barcode and value is barcode frequency
		cursor = read_PL.after(cursor)
	# compare barcode among two barcode, compute complexity: O(n^2)
	srt_barcode_ar = sorted(barcode_seq.keys())  # sorted barcode array
	head_dt, tail_dt = {}, {}  # head 6bp and tail 5bp dict, key is 6/5bp and value is 12bp whole barcode
	for barcode in barcode_seq.keys():  # key is barcode and value will correct barcode
		head_6, tail_5 = barcode[:HELF_BC_LEN], barcode[HELF_BC_LEN:]  # head 6bp(6=12/2) and tail 6bp include R1/R2 information
		for kmer, kmer_dt in ((head_6, head_dt), (tail_5, tail_dt)):
			if barcode_fre[barcode] == 0:  # already corrected 
				continue
			if kmer in kmer_dt:  # same head 6 bp to previous barcode
				prev_barcode, this_barcode = kmer_dt[kmer], barcode  # origin barcode
				if prev_barcode.split(':')[-1] == this_barcode.split(':')[-1] and Levenshtein(prev_barcode, this_barcode) <= EDIT_DIST_CUTOFF:  # edit distance < 2, .split(':')[-1] is R1/R2
					if barcode_fre[prev_barcode] > barcode_fre[this_barcode]:  # choose a barcode in higher frequency
						right_barcode = prev_barcode
						error_barcode = this_barcode
					else:
						right_barcode = this_barcode
						error_barcode = prev_barcode
					barcode_seq[error_barcode] = right_barcode  # value of barcode_seq is corrected barcode
					barcode_fre[right_barcode] += barcode_fre[error_barcode]  # key is barcode and value is barcode frequency
					barcode_fre[error_barcode] = 0
					kmer_dt[kmer] = right_barcode
			else:
				kmer_dt[kmer] = barcode
	# reset corrected barcode to BcRead class
	output = {}  # corrected barcode frequncy > MERGE_CUTOFF
	cursor = read_PL.first()
	while cursor is not None:
		old_cursor = cursor
		cursor = read_PL.after(cursor)
		this_read = old_cursor.get_element().get_read()
		this_barcode = old_cursor.get_element().get_barcode()
		this_cor_barcode = barcode_seq[this_barcode]
		if this_cor_barcode != this_barcode:  # correct barcode
			BC_correct_txt.write('%s=>%s\n' % (this_cor_barcode, this_barcode))
			cor_read += 1  # corrected reads number +1
			old_cursor.get_element().set_barcode(this_cor_barcode)  # reset correct barcode
		if barcode_fre[this_cor_barcode] >= MERGE_CUTOFF:  # corrected barcode frequncy > MERGE_CUTOFF
			if this_cor_barcode not in output:
				output.setdefault(this_cor_barcode, [])
			output[this_cor_barcode].append(read_PL.delete(old_cursor).get_read())  # pysam alignment object
	write_BAM(output, stat, pass_bam, no_pass_bam, result_BC)
	return cor_read

def output_BAM(read_PL, no_pass_bam, pass_bam, chrom, cluster_chrom, cluster_start, cluster_end, stat, result_BC):
	'''clear positional list of reads and output flitered reads.'''
	output = {}  # key is barcode and value is reads
	if cluster_chrom != chrom:  # fininsh one chromosome, output all cache reads.
		cursor = read_PL.first()
		while cursor is not None:
			old_cursor = cursor
			cursor = read_PL.after(cursor)
			this_read = old_cursor.get_element().get_read()
			this_barcode = old_cursor.get_element().get_barcode()
			output.setdefault(this_barcode, [])
			output[this_barcode].append(read_PL.delete(old_cursor).get_read())  # pysam alignment object
	else:  # output previous cluster reads
		cursor = read_PL.first()
		while cursor is not None:
			old_cursor = cursor
			cursor = read_PL.after(cursor)
			if old_cursor.get_element().get_read().reference_start < cluster_start-SLIDE_CUTOFF:
				this_read = old_cursor.get_element().get_read()
				this_barcode = old_cursor.get_element().get_barcode()
				output.setdefault(this_barcode, [])
				output[this_barcode].append(read_PL.delete(old_cursor).get_read())  # pysam alignment object
	write_BAM(output, stat, pass_bam, no_pass_bam, result_BC)
	return True

def write_BAM(output, stat, pass_bam, no_pass_bam, result_BC):
	for barcode in output:  # i is array of pysam read in same barcode 
		i = output[barcode]
		stat[len(i)] = stat.get(len(i), 0) + 1  # k is barcode member number
		if len(i) < MERGE_CUTOFF:  # 1,2,3,4
			for j in i:  # j is pysam read in same barcode
				no_pass_bam.write(j)
		else:  # >= 5
			merged_read = consensus_calling(barcode, i)
			result_BC[barcode] = result_BC.get(barcode, 0) + 1  # result fragment will in same barcode
			pass_bam.write(merged_read)
	return True

def consensus_calling(barcode, read_ar):  # read_ar is BcRead object array
	'''call consensus and return a bam read'''
	base2index = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
	index2base = ['A', 'C', 'G', 'T', 'N']
	# construct a BAM-read
	merged_read = pysam.AlignedSegment()
	merged_read.query_name = '%s:%d %s:%d' % (barcode, len(read_ar), read_ar[0].reference_name, read_ar[0].reference_start)
	merged_read.reference_id = read_ar[0].reference_id
	merged_read.mapping_quality = 60
	#merged_read.flag = 16 if read_ar[0].is_reverse is True else 0  # as singel ends
	merged_read.flag = read_ar[0].flag  # as singel ends
	# call correct start, end , cigar
	seq_ar, start_ar = [], []  # sequence and reference positional start of reads
	start_dt, end_dt, len_dt, cigar_dt = {}, {}, {}, {}  # start, end
	for read in read_ar:  # to find correct start, end and cigar
		end_dt[read.reference_end] = end_dt.get(read.reference_end, 0) + 1
		ref_pos_start, ref_pos_end = read.reference_end, read.reference_start
		read_cigar = read.cigartuples
		if read_cigar[0][0] == 4:  # soft clipping
			read_cigar.pop(0)
		if read_cigar[-1][0] == 4:
			read_cigar.pop()
		read_start, read_end = 0, 0  # read position
		for x in read.get_aligned_pairs(matches_only=True, with_seq=True):
			read_pos, ref_pos, ref_seq = x
			if ref_pos is not None and read_pos is not None and ref_pos is not None:
				if ref_pos_start > ref_pos:
					ref_pos_start, read_start = ref_pos, read_pos
				if ref_pos_end < ref_pos:
					ref_pos_end, read_end = ref_pos, read_pos
		if ref_pos_start == ref_pos_end:  # all softclipping reads
			continue
		seq = read.query_sequence[read_start : read_end+1]
		read_cigar.insert(0, len(seq))  # insert ref pos of read start to front of array
		read_cigar.insert(0, ref_pos_start)  # insert ref pos of read start to front of array
		read_cigar = tuple(read_cigar)
		cigar_dt[read_cigar] = cigar_dt.get(read_cigar, 0) + 1
		seq_ar.append(seq)
		start_ar.append(ref_pos_start)
	con_cigar = max(cigar_dt, key=cigar_dt.get)  # consensus cigar
	con_start, con_len, con_cigar = con_cigar[0], con_cigar[1], con_cigar[2:]
	merged_read.reference_start = con_start
	merged_read.cigartuples = con_cigar
	# call consensus sequence
	base_ar = [ [0, 0, 0, 0, 0] for x in range(0, con_len)]  # item in array is position
	con_seq, con_qual = '', []  # result sequence and quality of the reads in same barcode and coordinate
	agv_member = 0
	while seq_ar: 
		read_seq = seq_ar.pop()
		read_start = start_ar.pop()
		agv_member += 1
		base_ar_start = 0  #
		if read_start < con_start:  #read slide to left
			read_seq = read_seq[con_start-read_start:]
		if read_start > con_start:  #read slide to left
			base_ar_start = read_start-con_start  # start from some base later
		for x in range(0, min(len(read_seq), con_len)):
			base = read_seq[x]
			index = base2index[base]
			if x+base_ar_start < min(len(read_seq), con_len):  # avoid index out of range
				base_ar[x+base_ar_start][index] += 1
	for x in range(0, con_len):
		if sum(base_ar[x]) == 0:  # always the end of read is empty
			continue
		base = index2base[base_ar[x].index(max(base_ar[x]))]  # ATCG
		seq_error = (sum(base_ar[x])-max(base_ar[x]))/sum(base_ar[x])
		qual = base_q(seq_error)  # ascii
		con_seq += base
		con_qual.append(qual)
	merged_read.query_sequence = con_seq
	merged_read.query_qualities = con_qual
	merged_read.template_length = len(con_seq)  # as single ends
	#merged_read.tags = read.tags
	#merged_read.next_reference_id = None  # in same chromosome
	#merged_read.next_reference_start = None
	return merged_read

def base_q(seq_error):
	import math
	if seq_error == 0:
		result = 40
	else:
		result = min(int(-10*math.log10(seq_error)), 40)
	return result

def Levenshtein(seq1, seq2):
	'''Levenshtein/Edit Distance'''
	rows = len(seq1)+1
	cols = len(seq2)+1
	dist = [[0 for x in range(cols)] for x in range(rows)]
	for i in range(1, rows):
		dist[i][0] = i
	for i in range(1, cols):
		dist[0][i] = i

	for col in range(1, cols):
		for row in range(1, rows):
			if seq1[row-1] == seq2[col-1]:
				cost = 0
			else:
				cost = 1
			dist[row][col] = min(dist[row-1][col] + 1,       # deletion
								 dist[row][col-1] + 1,       # insertion
								 dist[row-1][col-1] + cost)  # substitution
	return dist[row][col]

class PositionalList(object):
	'''A sequential container of elements allowing positional access.'''

	##### Position class#####
	class Position(object):
		'''An abstraction representing the location of a single element.'''
		def __init__(self, container, node):
			'''Constructor should not be invoked by user.'''
			self.__container = container  # instance of PositionList class
			self.__node = node  # instance of _Node class

		def get_container(self):
			return self.__container

		def get_node(self):
			return self.__node

		def get_element(self):
			'''Return the element stored at this Position.'''
			return self.get_node().get_element()

		def __eq__(self, other):
			'''Return True if other is a Position represeting the same location.'''
			return type(other) is type(self) and other.get_node() is self.get_node()

		def __ne__(self, other):
			'''Retrun True if other does not represent the same loaction.'''
			return not (self == other)

	##### utility method  #####
	def __validate(self, p):
		'''Return position's node, or raise approprate error if invalid.'''
		if not isinstance(p, self.Position):
			raise TypeError('p must be proper Position type')
		if p.get_container() is not self:
			raise ValueError('p does not belong to this container')
		if p.get_node().get_next() is None:
			raise ValueError('p is no longer valid')
		return p.get_node()

	def __make_position(self, node):
		'''Return Position instance for given node (or None if sentinel).'''
		if node is self.__header or node is self.__trailer:
			return None
		return self.Position(self, node)

	##### _Node class  #####
	class _Node(object):
		'''Lightweigth, nonpublic class for storing a double linked node.'''
		__slots__ = '__element', '__prev', '__next'

		def __init__(self, e, p, n):
			self.__element = e
			self.__prev = p
			self.__next = n

		def get_prev(self):
			return self.__prev

		def get_next(self):
			return self.__next

		def get_element(self):
			return self.__element

		def set_prev(self, p):
			self.__prev = p

		def set_next(self, n):
			self.__next = n

		def set_element(self, e):
			self.__element = e

	##### Positional list class  #####
	def __init__(self):
		'''Creat an empty list'''
		self.__header = self._Node(None, None, None)
		self.__trailer = self._Node(None, None, None)
		self.__header.set_next(self.__trailer)
		self.__trailer.set_prev(self.__header)
		self.__size = 0

	def __len__(self):
		'''Return the number of elements in the list.'''
		return self.__size

	def is_empty(self):
		'''Return True if the list is empty.'''
		return self.__size == 0

	##### accessors  #####
	def first(self):
		'''Return the first Position in the list (or None if list is empty).'''
		return self.__make_position(self.__header.get_next())

	def last(self):
		'''Return the first Position in the list (or None if list is empty).'''
		return self.__make_position(self.__trailer.get_prev())

	def before(self, p):
		'''Return the Position just before Position p (or None if p is first).'''
		node = self.__validate(p)
		return self.__make_position(node.get_prev())

	def after(self, p):
		'''Return the Position just after Position p (or None if p is last).'''
		node = self.__validate(p)
		return self.__make_position(node.get_next())

	def __iter__(self):
		'''Generatea forward iteration of the elements of the list.'''
		cursor = self.first()
		while cursor is not None:
			yield cursor.get_element()
			cursor = self.after(cursor)

	##### mutators  #####
	def __insert_between(self, e, predecessor, successor):
		'''Add element e between two existing nodes and return new node.'''
		newest = self._Node(e, predecessor, successor)
		predecessor.set_next(newest)
		successor.set_prev(newest)
		self.__size += 1
		return self.__make_position(newest)

	def __delete_node(self, node):
		'''Delete nonsentinel node from the list and returen its element.'''
		predecessor = node.get_prev()
		successor = node.get_next()
		predecessor.set_next(successor)
		successor.set_prev(predecessor)
		self.__size -= 1
		element = node.get_element()
		node.set_prev(None)
		node.set_next(None)
		node.set_element(None)
		return element

	def add_first(self, e):
		'''Insert element e at the font  of the list and return new Postion.'''
		return self.__insert_between(e, self.__header, self.__header.get_next())

	def add_last(self, e):
		'''Insert element e at the back of the list and return new position.'''
		return self.__insert_between(e, self.__trailer.get_prev(), self.__trailer)

	def add_before(self, p, e):
		'''Insert element e into list after Positon p and return new Postion.'''
		original = self.__validate(p)
		return self.__insert_between(e, original.get_prev(), original)

	def add_after(self, p, e):
		'''Insert element e into list after Position pand return new Position.'''
		original = self.__validate(p)
		return self.__insert_between(e, original, original.get_next())

	def delete(self, p):
		'''Remove and return the elemet at Position p.'''
		original = self.__validate(p)
		return self.__delete_node(original)

	def replace(self, p, e):
		'''
		Replase the element at Position p.
		Retrun the element formerly at Position p.
		'''
		original = self.__validate(p)
		old_value = orginal.get_element()
		original.set_element(e)
		return old_value

if __name__ == '__main__':
	import sys
	if len(sys.argv) < 3:
		raise ValueError('Usage: python3 make_consensus_bam.py <input BAM/SAM> <outprefix> (optional)<chr:start-end>')
		sys.exit()
	else:
		input_file = sys.argv[1]
		output_file = sys.argv[2]
		region_chrom, region_start, region_end = None, None, None  # init of chromosome start and end of assigned region
		if len(sys.argv) == 4:
			region_chrom, coor = sys.argv[3].split(':')
			region_start, region_end = coor.split('-')
			region_start, region_end = int(region_start), int(region_end)
		main(input_file, output_file, region_chrom, region_start, region_end)
