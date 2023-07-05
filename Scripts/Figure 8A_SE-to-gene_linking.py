#!/usr/bin/python3

#######################################################################
# 
# Made by KeunHong
# Revised by Jeong-Woon
#
# Purpose : gene-peak linking using correlation between expression and Epigenomic signal
# Dep Programs: 1.bedtools
#
# Updated : v1.3/ 200128 - changed tss info (Ensembl v99)
#           v2.0/ 200309 - applied conservative null model (TCGA ATAC strategy)
#                    added threads function
#			v2.1/ 200317 - list to numpy & pandas
#
#			v2.2/ 200322 - fail to use multithread in p-value calculation
#			v2.3/ 200323 - eventually used for / list to calculate p-value
#			v3.0/ 210531 - add distance column & update tss files (v102/2020-08)
#
#######################################################################

### Change this part --------------------------------------------------
tss_gene = "tss_info/tss_geneID_CanFam31v102.st.uq.bed"
tss_trans = "tss_info/tss_transcriptID_CanFam31v102.bed"

### import  -----------------------------------------------------------
import os, sys, random, math
import numpy as np
import pandas as pd
from pandas import DataFrame
import scipy as sp
import scipy.stats as stat
from time import time
from optparse import OptionParser
from concurrent.futures import ProcessPoolExecutor
# FDR package
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

#import rpy2.rinterface
#base = importr('base')
#print(base._libPaths())

stats = importr('stats')

### Analysis code -----------------------------------------------------
def read_input (check_gene_type, in_exp, in_signal):
	print("1. Reading expression and signal tables\n")
	# check gene type
	global in_tss_file
	tss = ""
	if check_gene_type == "gene":
		in_tss_file = tss_gene
	elif check_gene_type == "transcript":
		in_tss_file = tss_trans
	# make tss gene list
	pd_in_tss = pd.read_csv(in_tss_file, sep = '\t', header = None)
	pd_in_tss.index = pd_in_tss.iloc[:,3]
	## read RNA expression table
	pd_in_exp = pd.read_csv(in_exp, sep = '\t', header = 0, index_col = 0)
	## read ATAC signal table
	pd_in_signal = pd.read_csv(in_signal, sep = '\t', header = 0)
	pd_in_signal.index = pd_in_signal.iloc[:,0]

	return pd_in_tss, pd_in_exp, pd_in_signal

def check_sample_genes_overlap (pd_in_tss, pd_in_exp, pd_in_signal):
	print("2-1. check overlapped samples (RNA vs ATAC)")
	# 1. Check overlapped genes between tss and RNA expression table
	pd_in_tss_gene_rmdu = pd_in_tss.index.drop_duplicates()
	ls_idx_tss = pd_in_tss_gene_rmdu.values.tolist()
	ls_idx_exp = pd_in_exp.index.values.tolist()
	overlap_tss_exp = list(set(ls_idx_exp) & set(ls_idx_tss))
	print(f"# of genes: overlap: {len(overlap_tss_exp)} [exp: {len(set(ls_idx_exp)) - len(overlap_tss_exp)}/ tss: {len(set(ls_idx_tss)) - len(overlap_tss_exp)}]")
	
	pd_in_tss_st = pd_in_tss.loc[overlap_tss_exp]
	pd_in_exp_st_1st = pd_in_exp.loc[overlap_tss_exp]

	# 2. Check overlapped samples between RNA and ATAC expression table
	print("2-2. check overlapped genes (tss vs RNA)")
	sample_exp = pd_in_exp_st_1st.columns
	sample_signal = pd_in_signal.columns
	ls_overlap = list(set(sample_exp.tolist()) & set(sample_signal.tolist()))
	print(f"# of samples: overlap: {len(ls_overlap)} [excepted exp: {len(sample_exp) - len(ls_overlap)}/ signal: {len(sample_signal) - len(ls_overlap)}]\n")
	
	pd_in_exp_st_2nd = pd_in_exp_st_1st[ls_overlap]
	pd_in_signal_st = pd_in_signal[ls_overlap]

	return pd_in_tss_st, pd_in_exp_st_2nd, pd_in_signal_st

def signal_filter_log2_variance (pd_in_tss, pd_in_exp, pd_in_signal, pd_in_signal_bed, zero_line_cutoff, add_RNA_exp, add_ATAC_signal, variance_cutoff, out_suffix):
	
	print("3-1. Noise filtering")
	len_pd_in_exp, len_pd_in_signal = len(pd_in_exp.index), len(pd_in_signal.index)
	len_pd_in_exp_cut, len_pd_in_signal_cut = round(len_pd_in_exp * variance_cutoff), round(len_pd_in_signal * variance_cutoff)
	print(f"# of cut and total line: RNA {len_pd_in_exp_cut}/{len_pd_in_exp}, ATAC {len_pd_in_signal_cut}/{len_pd_in_signal}")

	num_zero_cut = round(len(pd_in_exp.columns) * (1 - zero_line_cutoff))
	# Check number of zero in each line
	ls_exp_nonzero_idx = []
	ls_signl_nonzero_idx = []
	for idx in pd_in_exp.index:
		exp = pd_in_exp.loc[idx]
		exp_zero = exp[exp == 0]
		if len(exp_zero) < num_zero_cut:
			ls_exp_nonzero_idx.append(idx)
	for idx in pd_in_signal.index:
			signal = pd_in_signal.loc[idx]
			signal_zero = signal[signal == 0]
			if len(signal_zero) < num_zero_cut:
				ls_signl_nonzero_idx.append(idx)
	pd_in_exp_zero_ft = pd_in_exp.loc[ls_exp_nonzero_idx]
	pd_in_signal_zero_ft = pd_in_signal.loc[ls_signl_nonzero_idx]

	print("3-2. add score and log2 conversion")
	pd_in_exp_log2p = np.log2(pd_in_exp_zero_ft + add_RNA_exp)
	pd_in_signal_log2p = np.log2(pd_in_signal_zero_ft + add_ATAC_signal)

	print("3-3. Low variance filtering\n")
	# RNA/ coefficient of variation, CV || relative standard deviation, RSD
	pd_in_exp_log2p['cov'] = pd_in_exp_log2p.std(axis = 1) / pd_in_exp_log2p.mean(axis = 1)
	pd_in_exp_log2p.fillna(0, inplace = True) #NaN to 0
	pd_in_exp_log2p_st = pd_in_exp_log2p.sort_values(by = ['cov'], ascending = False)

	# ATAC/ coefficient of variation, CV || relative standard deviation, RSD
	pd_in_signal_log2p['cov'] = pd_in_signal_log2p.std(axis = 1) / pd_in_signal_log2p.mean(axis = 1)
	pd_in_signal_log2p.fillna(0, inplace = True) #NaN to 0
	pd_in_signal_log2p_st = pd_in_signal_log2p.sort_values(by = ['cov'], ascending = False)

	# Variance filtering
	if len(pd_in_exp_log2p_st.index) > len_pd_in_exp_cut:
		pd_exp_percent = pd_in_exp_log2p_st.iloc[:(len_pd_in_exp_cut)]
	else:
		pd_exp_percent = pd_in_exp_log2p_st

	if len(pd_in_signal_log2p_st.index) > len_pd_in_signal_cut:
		pd_signal_percent = pd_in_signal_log2p_st.iloc[:(len_pd_in_signal_cut)]
	else:
		pd_signal_percent = pd_in_signal_log2p_st
	# Remove 'cov' column
	del pd_exp_percent['cov']
	del pd_signal_percent['cov']

	# Write tss file including only filtered gene list
	pd_in_tss_st = pd_in_tss.loc[pd_exp_percent.index]
	pd_in_tss_st.to_csv(f'tmp_tss_{out_suffix}.bed', sep = '\t', header = False, index = False)

	# make peak bed file in filtered peaks list
	pd_bed = pd_in_signal_bed.loc[pd_signal_percent.index]
	pd_bed.to_csv(f'tmp_peak_{out_suffix}.bed', sep = '\t', header = False, index = False)

	return pd_exp_percent, pd_signal_percent


def gene_peak_network (range_length, out_suffix):
	print(f"4. Making network between genes and peaks (range = {range_length})")
	# make connect file using bedtools
	command = f"bedtools window -a tmp_tss_{out_suffix}.bed -b tmp_peak_{out_suffix}.bed -w {range_length} > tmp_connect_{out_suffix}.txt"
	print(command)
	os.system(command)
	print("bedtools window end")
	# read "tmp_connect.txt" file and extract network info

	pd_connect_tmp = pd.read_csv(f'tmp_connect_{out_suffix}.txt', sep = '\t', header = None)
	pd_connect_rmdu = pd_connect_tmp[[3, 9]].drop_duplicates()
	pd_connect_rmdu.columns = ['gene', 'signal']

	print(f"# of uniq and total links: {len(pd_connect_rmdu.index)} / {len(pd_connect_tmp.index)}")
	print(f"# of genes and peaks in links: {len(pd_connect_rmdu['gene'].drop_duplicates())} / {len(pd_connect_rmdu['signal'].drop_duplicates())}\n")

	return pd_connect_rmdu

def prepare_null_model (pd_exp_in, pd_signal_in, ran_genes):
	pd_in_genes = pd.Series(pd_exp_in.index)
	print("5. Calculation of mean and sd of correlation(consevative null model)\n")
	rt_pr_sp = pool.map(map_null_model, pd_in_genes) #[[geneA, mean, sd], [geneB, mean, sd], .... ]
	return list(rt_pr_sp)

def map_null_model(gene):
	np_corr_pr = np.array([])
	np_corr_sp = np.array([])

	ls_signal_idx_cp = ls_signal_idx.copy()
	## remove peaks within linking region
	for con in ls_connect:
		if con[0] == gene:
			if con[1] in ls_signal_idx_cp:
				ls_signal_idx_cp.remove(con[1])

	## random sampling
	list_random = random.sample(ls_signal_idx_cp, num_random_genes)

	for ran in list_random:
		pd_exp_gene = pd.Series(pd_exp_st.loc[gene], dtype=np.float64)
		pd_signal_peak = pd.Series(pd_signal_st.loc[ran], dtype=np.float64)
		out_corr_pr = pd.Series.corr(pd_exp_gene, pd_signal_peak, method='pearson')
		out_corr_sp = pd.Series.corr(pd_exp_gene, pd_signal_peak, method='spearman')
		np_corr_pr = np.append(np_corr_pr, out_corr_pr)
		np_corr_sp = np.append(np_corr_sp, out_corr_sp)

	out_pr_sp = [gene, str(round(np.mean(np_corr_pr), 5)), str(round(np.std(np_corr_pr), 5)), str(round(np.mean(np_corr_sp), 5)), str(round(np.std(np_corr_sp), 5))]
	return out_pr_sp

def make_results (pd_exp_in, pd_signal_in):
	print("6. calculate gene-peak correlation and p-value\n")
	# calculate correlation and significance p-value
	for con in ls_connect:
		pd_exp_gene = pd.Series(pd_exp_in.loc[con[0]], dtype=np.float64)
		pd_signal_peak = pd.Series(pd_signal_in.loc[con[1]], dtype=np.float64)
		
		out_corr_pr = pd.Series.corr(pd_exp_gene, pd_signal_peak, method='pearson')
		out_corr_sp = pd.Series.corr(pd_exp_gene, pd_signal_peak, method='spearman')

		pr_mean = pd_null_model.loc[con[0], 'pr_mean']
		pr_std = pd_null_model.loc[con[0], 'pr_std']
		sp_mean = pd_null_model.loc[con[0], 'sp_mean']
		sp_std = pd_null_model.loc[con[0], 'sp_std']

		null_pr = sp.stats.norm(loc = float(pr_mean), scale = float(pr_std))
		null_sp = sp.stats.norm(loc = float(sp_mean), scale = float(sp_std))

		pr_pv = pvalue(null_pr.cdf(float(out_corr_pr)))
		sp_pv = pvalue(null_sp.cdf(float(out_corr_sp)))

		con.append(out_corr_pr)
		con.append(out_corr_sp)
		con.append(pr_pv)
		con.append(sp_pv)

def pvalue (value):
	if value >= 0.5:
		return (1 - value)*2
	else:
		return value*2

def fdr ():
	print("7. Calculation of adjust.pvalue\n")
	pd_pvalue = pd.DataFrame(ls_connect, columns = ['gene', 'peak', 'pr_cor', 'sp_cor', 'pr_p', 'sp_p'])

	pr_q = pd.DataFrame(list(stats.p_adjust(list(pd_pvalue['pr_p']), method = 'BH')))
	sp_q = pd.DataFrame(list(stats.p_adjust(list(pd_pvalue['sp_p']), method = 'BH')))

	pd_pvalue['pr_q'], pd_pvalue['sp_q'] = pr_q, sp_q
	pd_pvalue_rd = round(pd_pvalue, 7)
	return pd_pvalue_rd

def cal_distance (pd_tss_in, pd_signal_in, pd_pvalue_in):
	#ls_tss = []
	#ls_peak = []
	ls_dist = []
	#type_dist=[]
	for gene, peak in zip(pd_pvalue_in['gene'], pd_pvalue_in['peak']):
		pos_tss = (pd_tss_in.loc[gene][1]+pd_tss_in.loc[gene][2])/2
		pos_peak = (pd_signal_in.loc[peak][1]+pd_signal_in.loc[peak][2])/2
		#ls_tss.append(gene)
		#ls_peak.append(peak)
		dist = abs(pos_tss - pos_peak)
		#type_dist.append(type(dist))

		if type(dist) != np.float64:
			ls_dist.append(min(dist.values.tolist()))
		else:
			ls_dist.append(dist)

	#print(set(type_dist))
	#pd_pvalue_in['tss_test'] = pd.DataFrame(list(ls_tss))
	#pd_pvalue_in['peak_test'] = pd.DataFrame(list(ls_peak))
	pd_pvalue_in['distance'] = pd.DataFrame(list(ls_dist)).astype(int)

	#print(ls_dist)
	return pd_pvalue_in


def write_output (list, out_file_name):
	with open(out_file_name, "w") as w:
		for i in range(len(list)):
			if i > 0:
				w.write("\n")
			for j in range(len(list[i])):
				if j > 0:
					w.write("\t")
				w.write(list[i][j])

### Main code ----------------------------------------------------------
if __name__ == '__main__':

	usage = """usage: %prog <Exp> <Signal> <Gene> <Range> <out_prefix>[options]
	Exp:    expression table (x-axis: gene or transcript IDs / y-axis: sample IDs)
	Signal: signal table (x-axis: peak IDs / y-axis: label, chr, start, end and sample IDs)
	Gene:   gene ID type (gene or transcript)
	Range:  range limit from tss which connect gene and peak (ex, 5000)
	out_prefix: prefix of output file name
	"""
	parser = OptionParser(usage=usage)
	#parser.add_option("-e", dest="exp_value", help="Cutoff of max expression value (default = 0.1)", default=0.1)
	#parser.add_option("-s", dest="signal_value", help="Cutoff of max signal value (default = 0.1)", default=0.1)
	parser.add_option("-z", dest="zero_filter", help="Cutoff ratio of degree containing non-zero in row lines (default = 0.2)", default=0.2)
	parser.add_option("-v", dest="variance", help="Cutoff ratio of high variance range 0 - 1 (default = 0.75)", default=0.75)
	parser.add_option("--ae", dest="add_RNA_value", help="Add specific expression value (default = 1)", default=1)
	parser.add_option("--as", dest="add_ATAC_value", help="Add specific signal value (default = 0)", default=0)


	parser.add_option("--null", dest="null_model_file", help="Null model file (default = False)", default=False)
	parser.add_option("-n", dest="num_genes", help="# of peaks used in null model (default = 500)", default=500)

	parser.add_option("-t", dest="threads", help="# of threads (default = 1)", default=1)
	
	(options, args) = parser.parse_args()
	if len(args) != 5:
		parser.print_help()
		sys.exit(1)
	# input variables
	in_exp_file = sys.argv[1]
	in_signal_file = sys.argv[2]
	gene_type = sys.argv[3]
	range_length = int(sys.argv[4])
	out_prefix = sys.argv[5]

	#exp_cut = float(options.exp_value)
	#signal_cut = float(options.signal_value)
	zero_line_cut = float(options.zero_filter)
	variance_cut = float(options.variance)
	add_RNA_score = float(options.add_RNA_value)
	add_ATAC_score =  float(options.add_ATAC_value)

	input_null_model = options.null_model_file
	num_random_genes = int(options.num_genes)

	num_threads = int(options.threads)
	# output_file_suffix
	suffix = f'l{range_length}_z{zero_line_cut}_v{variance_cut}_ae{add_RNA_score}_as{add_ATAC_score}_n{num_random_genes}'

	# global df and ls
	global pd_exp_st
	global pd_signal_st
	global pd_connect

	global ls_signal_idx
	global ls_connect
	
	global pd_null_model
	# check inputs number and type

	if gene_type != "gene" and gene_type != "transcript":
		print('Please write "gene" or "transcript" in Gene part\n')
		parser.print_help()
		sys.exit(1)

	## Start analysis ---------------------------------------------------
	start = time()
	# read input files
	pd_tss, pd_exp, pd_signal = read_input (gene_type, in_exp_file, in_signal_file)
	pd_signal_bed = pd_signal.iloc[:,[1,2,3,0]]

	# check overlapped samples (RNA vs ATAC) and genes (tss vs RNA)
	pd_tss_st, pd_exp_st, pd_signal_st = check_sample_genes_overlap (pd_tss, pd_exp, pd_signal)

	pd_exp_ft, pd_signal_ft = signal_filter_log2_variance (pd_tss_st, pd_exp_st, pd_signal_st, pd_signal_bed, zero_line_cut, add_RNA_score, add_ATAC_score, variance_cut, suffix)

	# make gene-peak network information
	pd_connect = gene_peak_network (range_length, suffix)

	# Calculate mean and std of random sampled gene-peak for null model
	mid_null = time()
	pool = ProcessPoolExecutor(max_workers = num_threads)
	ls_signal_idx = pd_signal_st.index.values.tolist()
	ls_connect = pd_connect.values.tolist()

	###################### input_null_model (X) ->> check null model file in directory
	if input_null_model == False:
		ls_null_pr_sp = prepare_null_model (pd_exp_st, pd_signal_bed, num_random_genes)

		write_output (ls_null_pr_sp, f'input_null_model_{suffix}')
		in_nm_name = f'input_null_model_{suffix}'
	else:
		in_nm_name = input_null_model

	# read null model information	
	pd_null_model = pd.read_csv(in_nm_name, sep = "\t", header = None, names = ['gene', 'pr_mean', 'pr_std', 'sp_mean', 'sp_std'], index_col = 'gene')
	mid_p = time()
	# calculate gene-peak correlation and p-value
	make_results (pd_exp_st, pd_signal_st)

	pd_qvalue = fdr ()
	# Distance 
	pd_qvalue = cal_distance(pd_tss_st, pd_signal_bed, pd_qvalue)
	print(pd_qvalue)

	# Write output file
	pd_qvalue.to_csv(f'{out_prefix}_{suffix}.txt', sep = '\t', index = False, )

	end = time()
	print('%.3f seconds : ALL' % (end - start))
	print('%.3f seconds : Null' % (mid_p - mid_null))
	print('%.3f seconds : Pval' % (end - mid_p))
	#timeNow ()
	print("-- END --")

	## Checking points
	# 1. variables location (global, local)
	# 2. decision rule of variables name
	# 3. abbreviation of variables
	# 4. multi-thread
	# 5. list <-> numpy, pandas
	# 6. pandas column / row recognition
