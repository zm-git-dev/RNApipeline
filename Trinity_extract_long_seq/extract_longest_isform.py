#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#####
#author:yueyao
#data:2017-07-06
#用于提取trinity组装之后的最长的转录本，并统计transcript和unigene的长度分布和作图
#python extract_longest_isform.py Trinity.fasta name
####

import re
import sys
import os
import numpy as np  
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt



def Stat_length(transcript_len):
	count_200,count_500,count_1000,count_2000=[0,0,0,0]
	for trans_len in transcript_len:
		trans_len=int(trans_len)
		if 200 <= trans_len <500 :
			count_200 =count_200 +1
		if trans_len < 1000 and trans_len >=500:
			count_500 = count_500 +1
		if trans_len <2000 and trans_len >=1000:
			count_1000 = count_1000 +1
		if trans_len >=2000:
			count_2000 = count_2000 +1
	total = len(transcript_len)
	stat_list = [count_200,count_500,count_1000,count_2000,total]
	return stat_list
	
def Calculate_N50(transcript_length,rate=0.5):
	Genome_totoal_length = sum(transcript_length)
	N50=0
	min_len = min(transcript_length)
	max_len = max(transcript_length)
	trans_num = len(transcript_length)
	mean_len = Genome_totoal_length / trans_num
	transcript_length.sort()
	transcript_length.reverse()
	for i_len in transcript_length:
		N50 = N50 + i_len
		if N50 >= Genome_totoal_length * rate:
			return [min_len,mean_len,max_len,i_len,Genome_totoal_length]

def Stat_Length_dis(trans_len):
	'''
	统计转录本的长度分布，用于作图，起始长度为200
	'''
	count=[]
	for i in range(2,21):
#		i is range from 2 to 20;
		count.append(0)
		for len in trans_len:
			len = int(len)
			if i == 20 and len >= i*100:
				count[i-2] = count[i-2] +1
			if 2<= i <= 19 and i*100 <= len < i*100+100:
				count[i-2] = count[i-2] +1
	return count

def draw_dis(count_plot,name,type):
	width=1
	plt.figure(figsize=(10,8))
	num = len(count_plot)+2
	num=int(num)
	idx = np.arange(2,num)
	id = idx * 100
	plt.bar(idx,count_plot,color='b')
	plt.xticks(range(2,num),id,rotation=40)
	plt.xlabel(type+" Length Interval")
	plt.ylabel("Number of "+type)
	plt.title(type+" Length Distribution")
	#plt.show()
	plt.savefig(name+'_'+type+'_length_distribution.png')
	plt.savefig(name+'_'+type+'_length_distribution.pdf')

	
with open (sys.argv[1],'r') as trinity_fa:
	transcript={}
	transcript_length=[]
	for line in trinity_fa.readlines():
		line = line.strip()
		m = re.match('^(>\w+)\s\w{3}=(\d+)\s',line)
		if m is not None:	
			gene_name = m.group(1)
			gene_len = m.group(2)
			gene_len = int(gene_len)
			transcript_length.append(gene_len)
			transcript[gene_name] = gene_len

name = sys.argv[2]

def output(name,transcript_length,type='Transcript'):
	#调用画图函数
	count_plot = Stat_Length_dis(transcript_length)
	draw_dis(count_plot,name,type)
	
	N50 = Calculate_N50(transcript_length,0.5)
	N70 = Calculate_N50(transcript_length,0.7)[3]
	N90 = Calculate_N50(transcript_length,0.9)[3]

	stat_trans_len = Stat_length(transcript_length)
	f1 = open (name+'_'+type+'_frequency_distribution.txt','w')
	f1.write("\t200-500bp\t500-1kb\t1k-2kb\t>2kbp\tTotal\n")
	f1.write(name+'\t'+str(stat_trans_len[0])+'\t'+str(stat_trans_len[1])+'\t'+str(stat_trans_len[2])+'\t'+str(stat_trans_len[3])+'\t'+str(stat_trans_len[4])+'\n')
	f1.close()
	
	f1 = open (name + '_'+type+'_length_distribution.txt','w')
	f1.write("\tmin length\tmean length\tmax length\tN50\tN70\tN90\tTotal\n")
	f1.write(name + '\t'+str(N50[0])+'\t'+str(N50[1])+'\t'+str(N50[2])+'\t'+str(N50[3])+'\t'+str(N70)+'\t'+str(N90)+'\t'+str(N50[4])+'\n')
	f1.close()

output(name,transcript_length)

def extract_longest_isform(transcript,name):
	gene_dit={}
	for gene in transcript.keys():
		gene_isform = gene.split('_')[-1]
		gene_id = gene.split('_')[0:4]
		gene_id = "_".join(gene_id)
		gene_len = transcript[gene]
		if gene_id in gene_dit.keys():
			if gene_dit[gene_id][1] >= gene_len:
				pass
			else:
				gene_dit[gene_id] = [gene_isform,gene_len,gene]
		else:
			gene_dit[gene_id] = [gene_isform,gene_len,gene]
#	unigene={}
	unigene_length=[]
	f1 = open (name + '_unigene_id.txt','w')
	for i_key in gene_dit.values():
			id = re.sub('>','',i_key[2])
			f1.write(id+'\t'+str(i_key[1])+'\n')
#			unigene[i_key[2]] = i_key[1]
			unigene_length.append(i_key[1])
	f1.close()
	output(name,unigene_length,type='Unigene')

extract_longest_isform(transcript,name)


#如果想要得到uigene的序列运行这一步即可
#os.system('perl ./common_bin/fishInWinter.pl -bf table -ff fasta '+name+'_unigene_id.txt '+sys.argv[1] +'> Unigene.fa')
