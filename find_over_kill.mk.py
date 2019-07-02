#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.1'

import os
import sys
import glob

input_file = sys.argv[1] if len(sys.argv) > 1 else './sample.info'

# write make file
idx, all_name = 0, ''
with open(os.path.basename(__file__)[:-3], 'w') as out, open('./sample.info', 'r') as f:
	for line in f:
		sam, vcf, raw_bam, cs_C = line.rstrip().split('\t')
		chrom = 'chrE'
		out.write('find_over_kill/%s.txt: %s %s\n' % (sam, vcf, raw_bam,))
		out.write('\tpython3 find_over_kill.py %s %s %s find_over_kill/%s.txt' %
				(vcf, raw_bam, chrom, sam))
		out.write('\n')
		all_name += ' find_over_kill/%s.txt' % (sam)

	out.write('all:%s\n' % (all_name))
