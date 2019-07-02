#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.1'

'''
Input bam of bwa-meth.
'''

import sys
import pysam

def main(vcf, bam, chrom, output_file):
	with pysam.VariantFile(vcf) as f,\
			pysam.AlignmentFile(bam, 'rb') as samfile,\
			open(output_file, 'w') as out:
		for rec in f.fetch(chrom):
			#print(rec.chrom, rec.pos, rec.ref, rec.alts)
			if len(rec.ref) == 1:
				in_vcf, out_vcf, tot_cov = 0, 0, 0
				frag_start_set, frag_end_set = set(), set()
				for read in samfile.fetch(chrom, rec.pos-1, rec.pos):
					if read.reference_start is not None and read.next_reference_start is not None and\
							read.infer_read_length() is not None:
						frag_start = min(read.reference_start, read.next_reference_start)
						frag_end = max(read.reference_start, read.next_reference_start) + read.infer_read_length()
						frag_start_set.add(frag_start)
						frag_end_set.add(frag_end)
					for cycle,ref_pos in read.get_aligned_pairs(matches_only=True):
						if ref_pos + 1 == rec.pos:
							alt_base = read.query_sequence[cycle].upper()
							tot_cov += 1
							if alt_base in rec.alts:
								in_vcf += 1
							else:
								out_vcf += 1
							break
				if out_vcf > 2 and out_vcf/tot_cov > 0.2:
					out.write('{}\t{}\t{}:{}-{}\t{}\t{}\n'.format(chrom, rec.pos, chrom,
						max(min(frag_start_set), rec.pos-200),
						min(max(frag_end_set), rec.pos+200),
						out_vcf,
						tot_cov))

if __name__ == '__main__':
	vcf, bam, chrom, output_file = sys.argv[1:]
	main(vcf, bam, chrom, output_file)
