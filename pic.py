#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.1'

import os
import glob

for coor_file in glob.glob('./find_over_kill/*txt'):
	sam = coor_file.split('/')[-1].split('.')[0]
	with open(coor_file, 'r') as f:
		for line in f:
			info = line.rstrip().split('\t')
			chrom, pos, coor = info[0], info[1], info[2]
			#print('./umi_vis.py vis	{} pic/{}.{}.{} ./raw_bam/{}.chrE.vcf.gz ./raw_bam/{}.mapped.bam'.format(coor, sam, chrom, pos, sam, sam))
			print(info)
			os.system('./umi_vis.py vis	{} pic/{}.{}.{} ./raw_bam/{}.chrE.vcf.gz ./raw_bam/{}.mapped.bam ./raw_bam/{}.CS.mapped.bam'.format(coor, sam, chrom, pos, sam, sam, sam))
