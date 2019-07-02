#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.1'

import glob

for bed_file in glob.glob('./*_raw.txt'):
	data = {}
	sam = bed_file.split('/')[-1].split('.')[0]
	with open(bed_file, 'r') as f, open('{}.srt'.format(bed_file), 'w') as out:
		for line in f:
			info = line.rstrip().split()
			chrom, start, end, umi_n, val = info[0], int(info[1]), int(info[2]), int(info[3]), int(info[4])
			data.setdefault(chrom, {})
			data[chrom].setdefault(start, [0, 0])
			if umi_n == 1:
				data[chrom][start][0] += val
			else:
				data[chrom][start][1] += val
	result = {}
	for chrom, vc in sorted(data.items()):
		for start, (umi_1, umi_2) in sorted(vc.items()):
			result.setdefault(chrom, {})
			result[chrom].setdefault(start, umi_2/umi_1)
			out.write('{}\t{}\t{}\n'.format(chrom,start, umi_2/umi_1))
