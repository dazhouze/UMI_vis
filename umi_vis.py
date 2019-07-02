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
import genomeview
import numpy as np
from random import shuffle
import matplotlib.cm as cm
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

# from genomeview
def color_by_nuc(base):
	colors = {"A":"blue", "C":"organge", "G":"green", "T":"black", "N":"gray"}
	#return colors.get(str(interval.variant.alts[0]).upper(), "gray")
	return colors.get(base.upper(), "gray")

class VCFTrack(genomeview.IntervalTrack):
	def __init__(self, vcf_path, name=None, chrom=None, start=None, end=None):
		super().__init__([], name=name)
		self.vcf = pysam.VariantFile(vcf_path)
		self.intervals = self
		self.color_fn = color_by_nuc
		self.row_height = 20
		self.min_variant_pixel_width = 2
		self.chrom = chrom
		self.start, self.end = start, end
	
	def __iter__(self):
		if self.start is None:
			start, end = self.scale.start, self.scale.end
			chrom = genomeview.match_chrom_format(self.scale.chrom, self.vcf.header.contigs)
		else:
			chrom = self.chrom
			start, end = self.start, self.end
		for variant in self.vcf.fetch(chrom, start, end):
			#print(dir(variant))
			interval = genomeview.Interval(
					variant.id if variant.id is not None else '{}:{}:{}'.format(variant.chrom, variant.start, variant.alts),
					variant.chrom,
					variant.start,
					variant.stop,
					None)
			interval.variant = variant
			yield interval
	
	def draw_interval(self, renderer, interval):
		# overriding this method isn't really necessary - we're just going to make
		# sure that every variant is at least several screen pixels wide, even
		# if we're zoomed pretty far out
		start = self.scale.topixels(interval.start)
		end = self.scale.topixels(interval.end)
		
		if end - start < self.min_variant_pixel_width:
			mid = (end + start) / 2
			start = mid - self.min_variant_pixel_width/2
			end = mid + self.min_variant_pixel_width/2
		
		#row = self.intervals_to_rows[interval.id]
		#row = 1  # fix to second line
		for idx,alt in enumerate(interval.variant.alts):
			row = self.intervals_to_rows[interval.id] + idx
			top = row * (self.row_height + self.margin_y)
			#color = self.color_fn(interval)
			color = self.color_fn(alt)
			yield from renderer.rect(start, top, end-start, self.row_height, fill=color,
    		                         **{"stroke":"none"})

##### visalization ######
def color_iter(num):
	num = int(num/3)+1
	x = np.arange(num)
	ys = [i+x+(i*x)**2 for i in range(num)]
	one, two, three = cm.spring(np.linspace(0, 0.6, len(ys))), cm.plasma(np.linspace(0.1, 0.8, len(ys))), cm.cool(np.linspace(0, 0.7, len(ys)))
	result = []
	for o,t,h in zip(one, two, three):
		result.append(o)
		result.append(t)
		result.append(h)
	return iter(result)

class ColorIter(object):
	def __init__(self, num):
		self.colors = color_iter(num)
		self.color = None
	def next_color(self): 
		c = next(self.colors)
		self.color = mc.to_hex(c[:3])
		return self.color
	def this_color(self):
		return self.color

def color_by_umi(interval):
	interval_len = interval.read.infer_read_length()
	frag_start = min(interval.read.reference_start, interval.read.next_reference_start)
	start_l, start_r = frag_start, max(interval.read.reference_start, interval.read.next_reference_start)
	umi = '{}'.format(interval.read.get_tag('RX'), start_l, start_r)
	if umi in umi_colors:
		return umi_colors[umi]
	return 'lightgrey'

def filter_by_umi(interval):
	interval_len = interval.infer_read_length()
	if interval.reference_start is None or\
			interval.next_reference_start is None or\
			interval_len is None:
		return False
	frag_start = min(interval.reference_start, interval.next_reference_start)
	frag_end = max(interval.reference_start, interval.next_reference_start)+interval_len
	if	frag_start >= start and\
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
				umi = '{}'.format(read.get_tag('RX'), start_l, start_r)
				count[umi] = count.get(umi, 0) + 1
	frag_draw, umi_draw, umi_all = 0, 0, 0
	for umi, umi_n in count.items():
		frag_draw += umi_n
		umi_all += 1
		if umi_n > 1:
			umi_draw += 1
	return count, frag_draw, umi_draw, umi_all

def umi_visualization(bams, chrom, start, end, output):
	doc = genomeview.Document(1000)
	global umi_colors
	umi_colors = {}
	# genome
	source = genomeview.FastaGenomeSource(REF_FA)
	gv = genomeview.GenomeView(chrom, max(0, start-150), end+150, "+", source)
	axis = genomeview.Axis()
	gv.add_track(axis)
	label_track = genomeview.track.TrackLabel('{}:{}-{}'.format(chrom, start, end))
	gv.tracks.insert(0, label_track)
	for bam in bams:
		# VCF/BAM track
		name = bam.split('/')[-1]
		if bam[-7:] == '.vcf.gz':   # VCF track
			variant_track = VCFTrack(bam, name, chrom, start, end)
			gv.add_track(variant_track)
		else:  # bam track
			track = genomeview.PairedEndBAMTrack(bam, name=name)
			gv.add_track(track)
			# format
			global count, colors
			count, frag_draw, umi_draw, umi_all = stats_umi(bam, chrom, start, end)
			colors = ColorIter(umi_draw)  # color generater
			umi_ar = list(count.keys())
			for umi in umi_ar:
				umi_n = count[umi]
				if umi_n > 1 and umi not in umi_colors:
					umi_colors[umi] = colors.next_color()
			print('{}.svg {}:{}-{} UMI:{} Fragments:{}'.
					format(output, chrom, start, end, umi_all, frag_draw))
			track.color_fn = color_by_umi
			track.include_read_fn = filter_by_umi  # exculde reads out of
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
		fig_y = np.array(fig_y) / tot_frag * 100
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
	result += '\tumi_vis.py <dis/vis> <chr:start-end> <output> [vcf.gz] <bam1> [bam2 ...]\n'
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
