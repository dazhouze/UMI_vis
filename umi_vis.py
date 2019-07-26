#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'ZHOU Ze <dazhouze@link.cuhk.edu.hk>'
__version__ = '0.2'

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
	colors = {"A":"blue", "C":"orange", "G":"green", "T":"black", "N":"gray"}
	return colors.get(base.upper(), "gray")

class VCFTrack(genomeview.IntervalTrack):
	def __init__(self, vcf_path, name=None, chrom=None, start=None, end=None):
		super().__init__([], name=name)
		self.vcf = pysam.VariantFile(vcf_path)
		self.intervals = self
		self.color_fn = color_by_nuc
		self.row_height = 15
		self.min_variant_pixel_width = 3
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
			row = idx
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
	colormaps = cm.spring(np.linspace(0, 0.6, len(ys))), cm.plasma(np.linspace(0.35, 0.8, len(ys))), cm.cool(np.linspace(0.2, 0.9, len(ys)))
	result = []
	for colormap in colormaps:
		for c in colormap:
			result.append(c)
	shuffle(result)
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

def filter_by_coor(interval):
	'''
	read_start = interval.reference_start
	read_end = interval.reference_end
	if	read_start >= start and\
			read_end <= end:
		return True
	return False
	'''
	return True

def stats_umi(bam, chrom, start, end):
	count, any_umi = {}, set()
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
			any_umi.add(read.get_tag('RX'))
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
	return count, frag_draw, umi_draw, umi_all, any_umi

def umi_visualization(bams, chrom, start, end, output, coor):
	doc = genomeview.Document(1000)
	global umi_colors
	umi_colors = {}
	# genome
	source = genomeview.FastaGenomeSource(REF_FA)
	gv = genomeview.GenomeView(chrom, max(0, start), end, "+", source)
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
			#track = genomeview.SingleEndBAMTrack(bam, name=name)
			gv.add_track(track)
			track.nuc_colors = {"A":"blue", "C":"orange", "G":"green", "T":"black", "N":"gray"}
			track.quick_consensus = False
			# format
			global count, colors
			count, frag_draw, umi_draw, umi_all, any_umi = stats_umi(bam, chrom, start, end)

			if len(coor.split(':')) == 3:
				colors = ColorIter(umi_draw)  # color generater
				umi_ar = list(count.keys())
				for umi in umi_ar:
					umi_n = count[umi]
					if umi_n > 1 and umi not in umi_colors:
						umi_colors[umi] = colors.next_color()
				track.color_fn = color_by_umi
				track.include_read_fn = filter_by_umi  # exculde reads out of
			else:
				colors = ColorIter(len(any_umi))  # color generater
				for umi in any_umi:
					umi_colors[umi] = colors.next_color()
				track.include_read_fn = filter_by_coor  # exculde reads out of
				track.color_fn = color_by_umi
				#track.color_fn = lambda x: "lightgray"
	doc.elements.append(gv)
	genomeview.save(doc, '{}.svg'.format(output))
	#genomeview.save(doc, '{}.pdf'.format(output), outformat='pdf')

##### distribution ######
def count_umi(bam, chrom, start, end):
	count_umi = {}  # same umi reads
	count_coor = {}  # same coor reads
	with pysam.AlignmentFile(bam, 'rb') as samfile:
		for read in samfile.fetch(chrom, start, end):
			if read.template_length > 0:  # half of reads
				continue
			if read.reference_start is None or\
					read.reference_end is None or\
					read.infer_read_length() is None:
				continue
			mid_pos = (read.reference_start + read.reference_end)/2
			if not start <= mid_pos < end:
				continue
			frag_start = min(read.reference_start, read.next_reference_start)
			frag_end = max(read.reference_start, read.next_reference_start)+read.infer_read_length()
			umi = read.get_tag('RX')
			umi = '{}:{}:{}:{}'.format(read.reference_name, frag_start, frag_end, umi)
			count_umi[umi] = count_umi.get(umi, 0) + 1
			coor = '{}:{}:{}'.format(read.reference_name, frag_start, frag_end)
			count_coor[coor] = count_coor.get(coor, 0) + 1
	dist_umi = {}  # umi distribution
	for umi, umi_n in count_umi.items():
		dist_umi[umi_n] = dist_umi.get(umi_n, 0) + 1
	dist_coor = {}  # coor distribution
	for coor, coor_n in count_coor.items():
		dist_coor[coor_n] = dist_coor.get(coor_n, 0) + 1
	return dist_umi, dist_coor

def umi_distribution(bam, chrom, start, end, output, coor):
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.backends.backend_pdf
	import seaborn as sns
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
	tot_dist_umi, tot_dist_coor = {}, {}  # umi count
	for t in threads_list:
		dist_umi, dist_coor = t.join()
		for umi_n, val in sorted(dist_umi.items()):
			tot_dist_umi[umi_n] = tot_dist_umi.get(umi_n, 0) + val
		for coor_n, val in sorted(dist_coor.items()):
			tot_dist_coor[coor_n] = tot_dist_coor.get(coor_n, 0) + val

	xlim, ylim = [], []
	for dt in (tot_dist_umi, tot_dist_coor):
		xmax, ymax = max(dt.keys()), max(dt.values()) / sum(dt.values()) * 100
		xlim.append(xmax)
		ylim.append(ymax)
	xlim, ylim = min(xlim), max(ylim)*1.1

	with matplotlib.backends.backend_pdf.PdfPages('{}.pdf'.format(output)) as pdf_all,\
			open('{}.txt'.format(output), 'w') as out:
		# umi distribution
		fig_x, fig_y = [], []
		for umi_n, val in sorted(tot_dist_umi.items()):
			fig_x.append(umi_n)
			fig_y.append(val/sum(tot_dist_umi.values())*100)
		tot_family = sum(tot_dist_umi.values())
		tot_reads = 0
		for coor_n, family_n in tot_dist_umi.items():
			tot_reads += coor_n * family_n

		fig = plt.figure()
		plt.plot(fig_x, fig_y, color='b')
		plt.xticks(fontsize=15)
		plt.yticks(fontsize=15)
		title = output.split('/')[-1].split('.')[0]
		plt.title('{} ({})\nFragments:{:,}\nFamily:{:,} Mean:{:.2f}'.
				format(title, coor, tot_reads, tot_family, tot_reads/tot_family), fontsize=18)
		out.write('umi\t{}\t{}\t{}\n'.format(tot_reads, tot_family, tot_reads/tot_family))
		plt.ylabel('Frequency (%)', fontsize=20)
		plt.xlabel('The no. of fragments per UMI cluster', fontsize=20)
		plt.tight_layout()
		pdf_all.savefig()
		plt.ylim(ymax = ylim)
		plt.xlim(-5, xlim)
		pdf_all.savefig()

		# same distribution
		fig_x, fig_y = [], []
		for coor_n, val in sorted(tot_dist_coor.items()):
			fig_x.append(coor_n)
			fig_y.append(val/sum(tot_dist_coor.values())*100)
		tot_family = sum(tot_dist_coor.values())
		tot_reads = 0
		for coor_n, family_n in tot_dist_coor.items():
			tot_reads += coor_n * family_n

		fig = plt.figure()
		plt.plot(fig_x, fig_y, color='b')
		plt.xticks(fontsize=15)
		plt.yticks(fontsize=15)
		title = output.split('/')[-1].split('.')[0]
		plt.title('{} ({})\nFragments:{:,}\nFamily:{:,} Mean:{:.2f}'.
				format(title, coor, tot_reads, tot_family, tot_reads/tot_family), fontsize=18)
		out.write('rmdup\t{}\t{}\t{}\n'.format(tot_reads, tot_family, tot_reads/tot_family))
		plt.ylabel('Frequency (%)', fontsize=20)
		plt.xlabel('The frequency of duplicated fragments', fontsize=20)
		plt.tight_layout()
		pdf_all.savefig()
		plt.ylim(ymax = ylim)
		plt.xlim(-5, xlim)
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
	result += '\tumi_vis.py dis chrX:1000-2000 ouput sample.bam\n'
	result += '\tumi_vis.py vis chrX:1000-2000 ouput sample.vcf sample.bam1 sample.bam2\n'
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
				start = max(0, start-50)
				end = end + 50
			else:
				start, end = int(info[0]), int(info[1])
	except:
		sys.exit('\n*** Incorrect coordinate {} ***\n{}'.format(coor, usage()))
	if method == 'dis':
		if len(bams) > 1:
			sys.exit('\n*** Incorrect bam file number (allow 1) {} ***\n{}'.format(bams, usage()))
		umi_distribution(bams[0], chrom, start, end, output, coor)
	else:
		umi_visualization(bams, chrom, start, end, output, coor)
