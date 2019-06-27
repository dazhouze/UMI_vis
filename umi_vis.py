#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'ZHOU Ze <dazhouze@link.cuhk.edu.hk>'
__version__ = '0.1'

'''
UMI is in BAM RX tag.
'''

import sys
import pysam
import threading
import numpy as np
import matplotlib.colors as mc

WIN_SIZE = 1000 # 1k
REF_FA = '/lustre/zhouze/Database/hg19_HBV_EBV.fa'

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

##### visalization ######
class ColorIter(object):
	def __init__(self, num):
		self.colors = ['darkblue', 'darkcyan', 'darkgoldenrod', 'darkgreen', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkturquoise', 'darkviolet', 'lightcoral', 'lightgreen', 'lightpink', 'lightsalmon', 'lightseagreen', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred']
		self.color = None
		self.pointer = 0
		self.len = len(self.colors)
	def next_color(self): 
		self.color = self.colors[self.pointer % self.len]
		self.pointer += 2
		return self.color
	def this_color(self):
		return self.color

def color_by_umi(interval):
	interval_len = interval.read.infer_read_length()
	frag_start = min(interval.read.reference_start, interval.read.next_reference_start)
	start_l, start_r = frag_start, max(interval.read.reference_start, interval.read.next_reference_start)
	umi = '{} {} {}'.format(interval.read.get_tag('RX'), start_l, start_r)
	if umi in count:
		return count[umi]
	return 'lightgrey'

def filter_by_umi(interval):
	interval_len = interval.infer_read_length()
	if interval.reference_start is None or\
			interval.next_reference_start is None or\
			interval_len is None:
		return False
	frag_start = min(interval.reference_start, interval.next_reference_start)
	frag_end = max(interval.reference_start, interval.next_reference_start)+interval_len
	#start_l, start_r = frag_start, max(interval.reference_start, interval.next_reference_start)
	#umi = '{} {} {}'.format(interval.get_tag('RX'), start_l, start_r)
	#if umi in count and\
	if		frag_start >= start and\
			frag_end <= end:
		return True
	return False

def stats_umi(bam, chrom, start, end):
	count = {}
	with pysam.AlignmentFile(bam, 'rb') as samfile:
		for read in samfile.fetch(chrom, start, end):
			if read.template_length > 0:
				continue
			read_len = read.infer_read_length()
			if read.reference_start is None or\
					read.next_reference_start is None or\
					read_len is None:
				continue
			frag_start = min(read.reference_start, read.next_reference_start)
			frag_end = max(read.reference_start, read.next_reference_start)+read_len
			if frag_start >= start and frag_end <= end:
				start_l, start_r = frag_start, max(read.reference_start, read.next_reference_start)
				umi = '{} {} {}'.format(read.get_tag('RX'), start_l, start_r)
				count[umi] = count.get(umi, 0) + 1
	frag_draw, umi_draw = 0, 0
	for umi, umi_n in count.items():
		frag_draw += umi_n
		umi_draw += 1
	return count, frag_draw, umi_draw

def umi_visualization(bams, chrom, start, end, output):
	import genomeview
	doc = genomeview.Document(1000)
	# genome
	source = genomeview.FastaGenomeSource(REF_FA)
	gv = genomeview.GenomeView(chrom, start, end, "+", source)
	axis = genomeview.Axis()
	gv.add_track(axis)
	for bam in bams:
		# BAM track
		name = bam.split('/')[-1]
		bam_track = genomeview.PairedEndBAMTrack(bam, name=name)
		gv.add_track(bam_track)
		label_track = genomeview.track.TrackLabel('{}:{}-{}'.format(chrom, start, end))
		gv.tracks.insert(0, label_track)
		# format
		global count, colors
		count, frag_draw, umi_draw = stats_umi(bam, chrom, start, end)
		colors = ColorIter(umi_draw)  # color generater
		umi_ar = list(count.keys())
		for umi in umi_ar:
			umi_n = count[umi]
			if umi_n < 2:
				del count[umi]
			else:
				count[umi] = colors.next_color()
		print('{}:{}-{} UMI:{} Fragments:{}'.format(chrom, start, end, umi_draw, frag_draw))
		bam_track.color_fn = color_by_umi
		bam_track.include_read_fn = filter_by_umi  # exculde reads out of
	doc.elements.append(gv)
	genomeview.save(doc, '{}.svg'.format(output))


##### distribution ######
def count_umi(bam, chrom, start, end):
	count = {}
	with pysam.AlignmentFile(bam, 'rb') as samfile:
		for read in samfile.fetch(chrom, start, end):
			if read.template_length > 0:
				continue
			if read.reference_start is None or read.reference_end is None:
				continue
			if not start <= (read.reference_start+read.reference_end)/2 < end:
				continue
			umi = read.get_tag('RX')
			count[umi] = count.get(umi, 0) + 1
	dist = {}
	for umi, umi_n in count.items():
		dist[umi_n] = dist.get(umi_n, 0) + 1
	return count, dist

def umi_distribution(bam, chrom, start, end, output):
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.backends.backend_pdf
	import seaborn as sns
	dist = {}  # umi count
	threads_list = []
	if end is None:
		with pysam.AlignmentFile(bam, 'rb') as samfile:
			start, end = 0, samfile.get_reference_length(chrom) - 1
	for l in range(int(start/WIN_SIZE)*WIN_SIZE, int(end/WIN_SIZE)*WIN_SIZE+1, WIN_SIZE):
		s, e = l, l+WIN_SIZE 
		s = start if s < start else s
		e = end if e > end else e
		if s >= e:
			continue
		t = ThreadWithReturn(target=count_umi,
				args=(bam, chrom, s, e),
				name='{}\t{}\t{}'.format(chrom, s, e))
		threads_list.append(t)
		t.start()
	with open('{}.txt'.format(output), 'w') as out:
		for t in threads_list:
			count, result = t.join()
			for umi_n, val in sorted(result.items()):
				dist[umi_n] = dist.get(umi_n, 0) + val
				out.write('{}\t{}\t{}\n'.format(t.name, umi_n, val))

	with matplotlib.backends.backend_pdf.PdfPages('{}.pdf'.format(output)) as pdf_all:
		fig_x, fig_y = [], []
		for umi_n, val in sorted(dist.items()):
			if umi_n > 10:
				break
			fig_x.append(umi_n)
			fig_y.append(val)
		fig = plt.figure()
		tot_frag = sum(dist.values())
		fig_y = np.array(fig_y) / tot_frag
		ax = sns.barplot(fig_x, fig_y,
				color = 'b',
			)
		plt.xticks(fontsize=15)
		plt.yticks(fontsize=15)
		title = output.split('/')[-1].split('.')[0]
		coor = chrom
		if start is not None:
			coor += ':{}'.format(start)
		if end is not None:
			coor += '-{}'.format(end)
		plt.title('{}\n{} #{}'.format(title, coor, tot_frag), fontsize=20)
		plt.ylabel('Frequency (%)', fontsize=20)
		plt.xlabel('UMI ', fontsize=20)
		plt.tight_layout()
		pdf_all.savefig()

def usage():
	result = '\nProgram: {}\n'.format(__file__)
	result += 'Version: {}\n'.format(__version__)
	result += 'Contact: {}\n'.format(__author__)
	result += '\nUsage:\n'
	result += '\tumi_vis.py <dis/vis> <chr:start-end> <output> <bam1> [bam2 ...]\n'
	result += '\nCommands:\n'
	result += '\tdis\tcheck the umi_distribution in given coordinate\n'
	result += '\tvis\tvisualize the UMI in given coordinate. PDF and JPG format is support.\n'
	result += '\nExample:\n'
	result += '\tumi_vis.py dis chrX:1000-2000 ouput.txt sample1.bam\n'
	result += '\tumi_vis.py vis chrX:1000-2000 ouput.pdf sample1.bam sample2.bam\n'
	return result

if __name__ == '__main__':
	# parameter check
	if len(sys.argv) < 5 or sys.argv[1] not in ('dis', 'vis'):  # untraced paramters
		sys.exit('\n*** Incorrect parameter ***\n{}'.format(usage()))
	method, coor, output, bams, = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]

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
	if method == 'dis':
		if len(bams) > 1:
			sys.exit('\n*** Incorrect bam file number (allow 1) {} ***\n{}'.format(bams, usage()))
		umi_distribution(bams[0], chrom, start, end, output)
	else:
		umi_visualization(bams, chrom, start, end, output)
