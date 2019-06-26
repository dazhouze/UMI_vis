#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'ZHOU Ze <dazhouze@link.cuhk.edu.hk>'
__version__ = '0.1'

import sys
import pysam
import threading
import genomeview

class ThreadWithReturn(threading.Thread):
	def __init__(self, group=None, target=None, name=None, args=(), kwargs={}, Verbose=None):
		threading.Thread.__init__(self, group, target, name, args, kwargs)
		self._return = None

	def run(self):
		if self._target is not None:
			self._return = self._target(*self._args, **self._kwargs)

	def join(self, *args):
		threading.Thread.join(self, *args)
		return self._return

def umi_visualization(bams, chrom, start, end, ouput_pic):
	pass

def umi_distribution(bams, chrom, start, end, ouput_pic):
	pass

def usage():
	result = '\nProgram: {}\n'.format(__file__)
	result += 'Version: {}\n'.format(__version__)
	result += 'Contact: {}\n'.format(__author__)
	result += '\nUsage:\n'
	result += '\tUMI_vis.py <dis/vis> <chr:start-end> <output> <bam1> [bam2 ...]\n'
	result += '\nCommands:\n'
	result += '\tdis\tcheck the umi_distribution in given coordinate\n'
	result += '\tvis\tvisualize the UMI in given coordinate. PDF and JPG format is support.\n'
	result += '\nExample:\n'
	result += '\tUMI_vis.py dis chrX:1000-2000 ouput.txt sample1.bam sample2.bam\n'
	result += '\tUMI_vis.py vis chrX:1000-2000 ouput.pdf sample1.bam sample2.bam\n'
	return result

if __name__ == '__main__':
	# parameter check
	if len(sys.argv) < 5 or sys.argv[1] not in ('dis', 'vis'):  # untraced paramters
		sys.exit('\n*** Incorrect parameter ***\n{}'.format(usage()))
	method, coor, ouput_pic, bams, = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]

	try:
		info = coor.split(':')
		if len(info) == 1:
			chrom, start, end = info[0], None, None
		else:
			chrom, co = info
			info = co.split('-')
			if len(info) == 1:
				start = int(info[0])
				end = start + 1 
			else:
				start, end = int(info[0]), int(info[1])
	except:
		sys.exit('\n*** Incorrect coordinate {} ***\n{}'.format(coor, usage()))
	print(method, chrom, start, end)
	if method == 'dis':
		umi_distribution(bams, chrom, start, end, ouput_pic)
	else:
		umi_visualization(bams, chrom, start, end, ouput_pic)
