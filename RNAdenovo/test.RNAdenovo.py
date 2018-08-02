#!/usr/bin/python
#_*_coding:utf-8_*_

'''
author: yueyao 
data:2017-06-05

参考转录组流程编写的python版本，所设置的参数适用于植物，参考橡胶项目的参数设置，暂时还不能生成pdf报告，功能注释分析的脚本还需要后期继续完善
'''

import re
import string
import os
import sys
import subprocess
import ConfigParser
import commands

outdir = os.getcwd()


def Check_dir(path):
	if os.path.exists(path):
		os.system('rm -rf ' + path)
		os.makedirs(path)
	else:
		os.makedirs(path)

###import lib information
lib = sys.argv[1]

####lib deal###
Lib = open(lib,'r')
lib_dict = {}
Group_dict = {}
for line in Lib.readlines():
	line=line.strip("\n")
    	lib_arry = line.split("\t")
	Keygroup = lib_arry[0]
	KeyName = lib_arry[1]
	if Keygroup in Group_dict.keys():
		Group_dict[Keygroup].append(KeyName)
	else:
		Group_dict[Keygroup] = []
		Group_dict[Keygroup].append(KeyName)
	FastqA = lib_arry[2].split(',')[0]
	FastqA = os.path.abspath(lib_arry[2].split(',')[0])
	FastqB = lib_arry[2].split(',')[1]
	FastqB = os.path.abspath(lib_arry[2].split(',')[1])
	lib_dict[KeyName]=[FastqA,FastqB] 
Lib.close()

##read conf file ###
conf = sys.argv[2]
config = ConfigParser.ConfigParser()
config.read(conf)
bin_path = config.get('software_path','bin_path')

###Check_soft_path
def Check_Software_Path(filepath):
    if os.path.exists(filepath):
        pass
    else:
        print 'filepath  do not exists, please check you conf file carefully!'
        sys.exit(0)



###step00 filter rawdata###
###用于设置是否进行过滤操作
filter = config.get('software_parameter','filter')

cleandata_dict={}
if filter == 'yes':
	print "Filter work shell was creating!"
	#脚本路径
	fqcheck_distribute = config.get('script_path','fqcheck_distribute.pl')
	soapnuke_stat = config.get('script_path','soapnuke_stat.pl')
	findNtile = config.get('script_path','findNtile.pl')
	
	fqcheck33 = config.get('software_path','fqcheck33')
	soapnuke = config.get('software_path','soapnuke')
	soapnuke_para = config.get('software_parameter','soapnuke_para')

	Check_dir(outdir + '/01.Filter_SOAPnuke')
	Check_dir(outdir + '/01.Filter_SOAPnuke/result')
	filter_sh_path = outdir + '/01.Filter_SOAPnuke/'
	
	for lib_name in lib_dict.keys():
		filter_result = outdir + '/01.Filter_SOAPnuke/result/' + lib_name
		Check_dir(filter_result)
		cleandata_dict[lib_name]=[filter_result+'/'+os.path.split(lib_dict[lib_name][0])[1]+'.clean',filter_result+'/'+os.path.split(lib_dict[lib_name][1])[1]+'.clean']
		cleandata_name1 = os.path.split(lib_dict[lib_name][0])[1]+'.clean'
		cleandata_name2 = os.path.split(lib_dict[lib_name][1])[1]+'.clean'
		filter_sh = open(filter_sh_path+lib_name+'_filter.sh','w')
		filter_stat_sh = open(filter_sh_path+lib_name+'_filter_stat.sh','w')
		#create filter statistic script
		filter_stat_sh.write("echo ==========start at : `date` ========== && \\\n")
		filter_stat_sh.write("%s -r %s -c %s && \\\n" % (fqcheck33,cleandata_dict[lib_name][0],filter_result+'/'+cleandata_name1+'1.fqcheck'))
		filter_stat_sh.write("%s -r %s -c %s && \\\n" % (fqcheck33,cleandata_dict[lib_name][1],filter_result+'/'+cleandata_name2+'2.fqcheck'))
		filter_stat_sh.write("perl %s %s %s -o %s && \\\n" % (fqcheck_distribute,filter_result+'/'+cleandata_name1+'1.fqcheck',filter_result+'/'+cleandata_name2+'2.fqcheck',filter_result))
		filter_stat_sh.write("echo ==========end at : `date` ========== && \\\n")
		filter_stat_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
		filter_stat_sh.write("echo Still_waters_run_deep > " + filter_sh_path+lib_name+'_filter_stat.sh.sign\n')
		#create filter shell script 
		filter_sh.write("echo ==========start at : `date` ========== && \\\n")
		filter_sh.write("tile=`perl %s -fq1 %s -fq2 %s -seqType 0` && \\\n" % (findNtile,lib_dict[lib_name][0],lib_dict[lib_name][1]))
		filter_sh.write('%s filter %s -1 %s -2 %s -o %s -C %s -D %s $tile && \\\n' % (soapnuke,soapnuke_para,lib_dict[lib_name][0],lib_dict[lib_name][1],filter_result,cleandata_name1,cleandata_name2))
		filter_sh.write("echo ==========end at : `date` ========== && \\\n")
		filter_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
		filter_sh.write("echo Still_waters_run_deep > " + filter_sh_path+lib_name+'_filter.sh.sign\n')
		
		filter_sh.close()
		filter_stat_sh.close()
	print "Filter work shell complete!"
else:
	print "You skip Filter raw data, if you need Filter raw data, you could set filter = yes."
	for lib_name in lib_dict.keys():
		cleandata_dict[lib_name]=[lib_dict[lib_name][0],lib_dict[lib_name][1]]

###step01 trinity denovo###
print "Denovo work shell was creating!"
Check_dir(outdir + '/02.Denovo_Trinity')
Check_dir(outdir + '/02.Denovo_Trinity/result')

trinity_sh_path = outdir + '/02.Denovo_Trinity'
trinity_result = outdir + '/02.Denovo_Trinity/result'

trinity = config.get('software_path','trinity')
trinity_para = config.get('software_parameter','trinity_para')
tgicl_para = config.get('software_parameter','tgicl_para')
trinity_script_path = ''

for lib_group in Group_dict.keys():
#	print lib_group
	Check_dir(trinity_result + '/'+lib_group+'/Trinity')
	fq_repeat=[]
	leftfq=[]
	rightfq=[]
	for i in range(0,len(Group_dict[lib_group])):
		fq_repeat.append(str(Group_dict[lib_group][i]))
	for j in fq_repeat:
		leftfq.append(cleandata_dict[j][0])
		rightfq.append(cleandata_dict[j][1])
	left_fq = ','.join(leftfq)
	right_fq = ','.join(rightfq)

	trinity_sh = open(trinity_sh_path + '/' + lib_group + '_Denovo.sh','w')
	trinity_sh.write("echo ==========start at : `date` ========== && \\\n")
	trinity_sh.write("export PATH=/ifshk4/BC_PUB/biosoft/PIPE_RD/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/software:/ifshk4/BC_PUB/biosoft/PIPE_RD/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/software/bowtie:$PATH && \\\n")
	trinity_sh.write("export LD_LIBRARY_PATH=/ifshk4/BC_PUB/biosoft/PIPE_RD/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/software/RNA_lib:$LD_LIBRARY_PATH && \\\n")
	trinity_sh.write('%s %s --left %s --right %s -output %s && \\\n' % (trinity,trinity_para,left_fq,right_fq,trinity_result+'/'+lib_group+'/Trinity'))
	trinity_sh.write('if [ ! -f %s/Trinity/Trinity.fasta ];then echo \"%s/Trinity/Trinity.fasta not exists! exit...\" && exit;fi && \\\n' % (trinity_result+'/'+lib_group,trinity_result+'/'+lib_group))
	trinity_sh.write("perl "+bin_path + "Denovo/get_chosen_fa.pl -fa " + trinity_result + '/'+lib_group+'/Trinity/'+'Trinity.fasta -output '+trinity_result + '/'+lib_group+'/'+lib_group+'.Trinity.fa' +' -type TR -name '+lib_group+" -size && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/fa_quality.pl -len -Head -N -gc "+trinity_result + '/'+lib_group+'/'+lib_group+".Trinity.fa && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/barplot.pl "+trinity_result + '/'+lib_group+'/'+lib_group+".Trinity.fa.quality.xls "+lib_group+".Trinity && \\\n")
	trinity_sh.write("mkdir -p "+trinity_result + '/'+lib_group+"/Tgicl && \\\n")
	trinity_sh.write("cd "+trinity_result + '/'+lib_group+"/Tgicl && \\\n")
	trinity_sh.write("ln -s ../Trinity/Trinity.fasta "+'./'+lib_group+".for_cluster.fa && \\\n")
	trinity_sh.write("perl "+bin_path + "software/tgicl/tgicl.nozmsort.pl "+lib_group+".for_cluster.fa "+tgicl_para+" && \\\n")
	trinity_sh.write("cat asm_*/align > align && \\\n")
	trinity_sh.write("perl "+bin_path + "software/tgicl/phrap.id.list.pl align align && \\\n")
	trinity_sh.write("perl "+bin_path + "software/tgicl/get_single.pl align.cluster "+lib_group+".for_cluster.fa single && \\\n")
	trinity_sh.write("cat asm_*/contigs > asm_cluster.fa && cat single.prefect.fa >> single.fa && \\\n")
	trinity_sh.write("cat asm_*/contigs  single.fa > tgicl_cluster_and_single.fa  && \\\n")
	trinity_sh.write("cat align.cluster single.list > tgicl_cluster_and_single.fa.list && \\\n")
	trinity_sh.write(bin_path+"software/formatdb -p F -i tgicl_cluster_and_single.fa && \\\n")
	trinity_sh.write(bin_path+"software/blastall -p blastn -m 8 -e 1e-10 -F F -a 10 -d tgicl_cluster_and_single.fa -i tgicl_cluster_and_single.fa -o all_vs_all.blast.m8 && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/cluster_for_coverage.pl all_vs_all.blast.m8 "+trinity_result + '/'+lib_group+"/Tgicl && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/clusterOrSingle.pl "+trinity_result + '/'+lib_group+"/Tgicl "+lib_group +"200 && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/get_chosen_fa.pl -fa "+lib_group+".all.fa -output "+lib_group+"-Unigene.fa && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/get_chosen_fa.pl -fa cluster.all.fa -output "+lib_group+"-Cluster.fa && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/get_chosen_fa.pl -fa singleton.all.fa -output "+lib_group+"-Single.fa  && \\\n")
	trinity_sh.write("cd ../ && cp Tgicl/"+lib_group+"-Unigene.fa ./ && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/fa_quality.pl -len -Head -N  -gc "+lib_group+"-Unigene.fa && \\\n")
	trinity_sh.write("perl "+bin_path + "Denovo/barplot.pl "+lib_group+"-Unigene.fa.quality.xls "+lib_group+"-Unigene && \\\n")
	trinity_sh.write("\n")
	trinity_sh.write("echo ==========end at : `date` ========== && \\\n")
	trinity_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
	trinity_sh.write("echo Still_waters_run_deep > " + trinity_sh_path+'/'+lib_group+'_Denovo.sh.sign\n')
	trinity_sh.close()

###Merged for denovo result###
Sample_List = list(Group_dict.keys())
if (len(Sample_List) == 1):
	print "Single Sample was used!"
	denovo_fa = trinity_result + '/'+str(Sample_List)+'/'+str(Sample_List)+'-Unigene.fa'
###stat_for_denovo_result###
	assembly_stat_sh = open(trinity_sh_path+'/'+'assembly_stat.sh','w')
	assembly_stat_sh.write("echo ==========start at : `date` ========== && \\\n")
	assembly_stat_sh.write("perl "+bin_path + "Denovo/assembly_stat.pl -indir "+trinity_result +" -outdir "+trinity_result+" && \\\n")
	assembly_stat_sh.write("echo ==========end at : `date` ========== && \\\n")
	assembly_stat_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
	assembly_stat_sh.write("echo Still_waters_run_deep > " + trinity_sh_path+'/assembly_stat.sh.sign\n')
	assembly_stat_sh.close()
else:
	print "Repeat Sample have been used!"
	denovo_fa = trinity_result + '/Merge_for_Tgicl/All-Unigene.fa'
	all_tgicl_sh_path = outdir + '/02.Denovo_Trinity'
	all_tgicl_sh_result = outdir + '/02.Denovo_Trinity/result'
	Check_dir(outdir + '/02.Denovo_Trinity/result/Merge_for_Tgicl/Tgicl')
	total_fa = ' '
	for lib_group in Group_dict.keys():
		total_fa = total_fa + trinity_result+'/'+lib_group+'/'+lib_group+'-Unigene.fa '
	all_tgicl_sh = open(all_tgicl_sh_path + '/' + 'all_tgicl.sh','w')
	all_tgicl_sh.write("echo ==========start at : `date` ========== && \\\n")
	all_tgicl_sh.write('cat ' + total_fa + ' >all.for_claster.fa && \\\n')
	all_tgicl_sh.write("perl "+bin_path + "software/tgicl/tgicl.nozmsort.pl all.for_cluster.fa "+tgicl_para+" && \\\n")
	all_tgicl_sh.write("cat asm_*/align > align && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "software/tgicl/phrap.id.list.pl align align && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "software/tgicl/get_single.pl align.cluster all.for_cluster.fa single && \\\n")
	all_tgicl_sh.write("cat asm_*/contigs > asm_cluster.fa && cat single.prefect.fa >> single.fa && \\\n")
	all_tgicl_sh.write("cat asm_*/contigs  single.fa > tgicl_cluster_and_single.fa  && \\\n")
	all_tgicl_sh.write("cat align.cluster single.list > tgicl_cluster_and_single.fa.list && \\\n")
	all_tgicl_sh.write(bin_path+"software/formatdb -p F -i tgicl_cluster_and_single.fa && \\\n")
	all_tgicl_sh.write(bin_path+"software/blastall -p blastn -m 8 -e 1e-10 -F F -a 10 -d tgicl_cluster_and_single.fa -i tgicl_cluster_and_single.fa -o all_vs_all.blast.m8 && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "Denovo/cluster_for_coverage.pl all_vs_all.blast.m8 "+trinity_result + '/'+"Merge_for_Tgicl/Tgicl && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "Denovo/clusterOrSingle.pl "+trinity_result + '/'+"Merge_for_Tgicl/Tgicl All +200 && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "Denovo/get_chosen_fa.pl -fa All.all.fa -output All-Unigene.fa && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "Denovo/get_chosen_fa.pl -fa cluster.all.fa -output All-Cluster.fa && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "Denovo/get_chosen_fa.pl -fa singleton.all.fa -output All-Single.fa  && \\\n")
	all_tgicl_sh.write("cd ../ && cp Tgicl/All-Unigene.fa ./ && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "Denovo/fa_quality.pl -len -Head -N  -gc All-Unigene.fa && \\\n")
	all_tgicl_sh.write("perl "+bin_path + "Denovo/barplot.pl All-Unigene.fa.quality.xls All-Unigene && \\\n")
	all_tgicl_sh.write("\n")
	all_tgicl_sh.write("echo ==========end at : `date` ========== && \\\n")
	all_tgicl_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
	all_tgicl_sh.write("echo Still_waters_run_deep > " + all_tgicl_sh_path+'/all_tgicl.sh.sign\n')
	all_tgicl_sh.close()
####stat_for_denovo_result###
	assembly_stat_sh = open(trinity_sh_path+'/'+'assembly_stat.sh','w')
	assembly_stat_sh.write("echo ==========start at : `date` ========== && \\\n")
	assembly_stat_sh.write("perl "+bin_path + "Denovo/assembly_stat.pl -indir "+trinity_result +" -outdir "+trinity_result+" && \\\n")
	assembly_stat_sh.write("echo ==========end at : `date` ========== && \\\n")
	assembly_stat_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
	assembly_stat_sh.write("echo Still_waters_run_deep > " + trinity_sh_path+'/assembly_stat.sh.sign\n')
	assembly_stat_sh.close()
print "Denovo shell complete!"
###SSR detection###
print "SSR work shell is creating!"
Check_dir(outdir + '/SSR_MISA')
Check_dir(outdir + '/SSR_MISA/result')

ssr_sh = outdir + '/SSR_MISA'
ssr_result_path = outdir + '/SSR_MISA/result'

ssr = config.get('software_path','ssr')
ssr_parameter = config.get('software_parameter','ssr_parameter')

os.system("perl "+ssr +" -fa "+denovo_fa + " -opts "+ssr_parameter +" -shdir "+ssr_sh+" -outdir "+ssr_result_path)
print "SSR work shell complete!"
###TF###
print "TF work shell was creating!"
Check_dir(outdir + '/TFpredict_Hmmsearch')
Check_dir(outdir + '/TFpredict_Hmmsearch/result')

tf_sh_path = outdir + '/TFpredict_Hmmsearch'
tf_result_path = outdir + '/TFpredict_Hmmsearch/result'

fasta_path,fasta_prefix = os.path.split(denovo_fa)
fasta_prefix = str(fasta_prefix)
fasta_prefix = re.sub(r'\.fa','',fasta_prefix)
#print fasta_prefix
tf_sh = open(tf_sh_path+'/'+'TF_predict.sh','w')
tf_sh.write("echo ==========start at : `date` ========== && \\\n")
tf_sh.write(bin_path+"software/interproscan5/bin/nucleotide/getorf -minsize 150 -sequence "+denovo_fa+' -outseq '+ tf_result_path + '/' + fasta_prefix+'.orf  && \\\n')
tf_sh.write("perl "+bin_path+"/Annotation/select_orf.pl -seq " + denovo_fa + ' -input '+ tf_result_path + '/' + fasta_prefix +'.orf ' + '-output ' + tf_result_path + '/' + fasta_prefix + '.pep  && \\\n')
tf_sh.write("perl "+bin_path+"/TF/TFCodingGene_Predict.pl " + tf_result_path + '/' + fasta_prefix + '.pep ' + fasta_prefix +' ' + tf_result_path + " && \\\n")
tf_sh.write("echo ==========end at : `date` ========== && \\\n")
tf_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
tf_sh.write("echo Still_waters_run_deep > " + tf_sh_path+'/tf_sh.sign\n')
tf_sh.close()
print "TF work shell complete!"

###SNP###
Check_dir(outdir + '/SNP_GATK')
Check_dir(outdir + '/SNP_GATK/result/index-build')
Check_dir(outdir + '/SNP_GATK/result')
snp_sh_path = outdir + '/SNP_GATK'
snp_result_path = outdir + '/SNP_GATK/result'
IndexDir = outdir + '/SNP_GATK/result/index-build'

hisat_para = config.get('software_parameter','hisat_para')
gatk_para = config.get('software_parameter','gatk_para')


snp_index_build_sh = open (snp_sh_path + '/index_build.sh','w')
snp_index_build_sh.write("export LD_LIBRARY_PATH=$Bin/software/RNA_lib:\$LD_LIBRARY_PATH && \\\n")
snp_index_build_sh.write("cd "+ IndexDir+" && \\\n")
snp_index_build_sh.write("perl "+bin_path+"/SNP/sort_fa.pl "+ denovo_fa +  " refMrna.fa && \\\n")
snp_index_build_sh.write("# build hisat index && \\\n")
snp_index_build_sh.write(bin_path+"/software/hisat/hisat-build refMrna.fa refMrna && \\\n")
snp_index_build_sh.write("# build gatk index && \\\n")
snp_index_build_sh.write(bin_path+"software/samtools faidx refMrna.fa && \\\n")
snp_index_build_sh.write(bin_path+"/software/picard/CreateSequenceDictionary.jar R=refMrna.fa O=refMrna.dict ")


for lib_name in sorted(cleandata_dict.keys()):
	Check_dir(snp_sh_path + '/'+ lib_name)
	Check_dir(snp_result_path + '/'+ lib_name)
	Check_dir(snp_sh_path + '/'+ lib_name + '/java_tmp')
	JavaTmpDir = snp_sh_path + '/'+ lib_name + '/java_tmp'

	snp_sh = open (snp_sh_path + '/'+ lib_name+'/'+lib_name+'_hisat_align.sh','w')
	snp_sh.write("export LD_LIBRARY_PATH=$Bin/software/RNA_lib:\$LD_LIBRARY_PATH && \\\n")
	snp_sh.write("cd "+snp_result_path + '/'+ lib_name+" && \\\n")
	snp_sh.write(bin_path+"software/hisat/hisat "+hisat_para+" -x "+IndexDir+'/refMrna -1 '+cleandata_dict[lib_name][0]+' -2 '+cleandata_dict[lib_name][1]+" -S "+lib_name+".sam \n")
	snp_sh.close()

	gatk_sh = open (snp_sh_path + '/'+ lib_name+'/'+lib_name+'_gatk.sh','w')
	gatk_sh.write("cd "+snp_result_path+'/'+lib_name+" && \\\n")
	gatk_sh.write(bin_path+"/software/java -Djava.io.tmpdir="+JavaTmpDir+" -jar "+bin_path+"/software/picard/AddOrReplaceReadGroups.jar I="+lib_name+".sam O="+lib_name+".addRG.bam RGID="+lib_name+" RGLB="+' wenkumingzi '+" RGPL=illumina RGPU=machine RGSM="+lib_name+" VALIDATION_STRINGENCY=SILENT && \\\n")
	gatk_sh.write(bin_path+"/software/java -Djava.io.tmpdir="+JavaTmpDir+" -jar "+bin_path+"/software/picard/ReorderSam.jar I="+lib_name+".addRG.bam O="+lib_name+".addRG.Reorder.bam R="+IndexDir+"/refMrna.fa VALIDATION_STRINGENCY=SILENT && \\\n")
	gatk_sh.write(bin_path+"/software/java -Djava.io.tmpdir="+JavaTmpDir+" -jar "+bin_path+"/software/picard/SortSam.jar I="+lib_name+".addRG.Reorder.bam O="+lib_name+".addRG.Reorder.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT && \\\n")
	gatk_sh.write(bin_path+"/software/java -Djava.io.tmpdir="+JavaTmpDir+" -jar "+bin_path+"/software/picard/MarkDuplicates.jar REMOVE_DUPLICATES=false I="+lib_name+".addRG.Reorder.Sort.bam O="+lib_name+".addRG.Reorder.Sort.markDup.bam METRICS_FILE="+lib_name+".addRG.Reorder.Sort.markDup.metrics VALIDATION_STRINGENCY=SILENT && \\\n")
	gatk_sh.write(bin_path+"/software/samtools index "+lib_name+".addRG.Reorder.Sort.markDup.bam && \\\n")
	gatk_sh.write(bin_path+"/software/java1.8 -Djava.io.tmpdir="+JavaTmpDir+" -jar "+bin_path+"/software/GenomeAnalysisTK.jar -T SplitNCigarReads -R "+IndexDir+"/refMrna.fa -I "+lib_name+".addRG.Reorder.Sort.markDup.bam -o "+lib_name+".addRG.Reorder.Sort.markDup.splitN.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS && \\\n")
	gatk_sh.write(bin_path+"/software/samtools index "+lib_name+".addRG.Reorder.Sort.markDup.splitN.bam && \\\n")
	gatk_sh.write(bin_path+"/software/java1.8 -Xmx20G  -Djava.io.tmpdir="+JavaTmpDir+" -jar "+bin_path+"/software/GenomeAnalysisTK.jar -T HaplotypeCaller -R "+IndexDir+"/refMrna.fa -I "+lib_name+".addRG.Reorder.Sort.markDup.splitN.bam "+ gatk_para +"-o "+lib_name+".gatk.vcf && \\\n")
	gatk_sh.write(bin_path+"/software/java1.8 -Djava.io.tmpdir="+JavaTmpDir+" -jar "+bin_path+"/software/GenomeAnalysisTK.jar -T SelectVariants -R "+IndexDir+"/refMrna.fa -selectType SNP -V "+lib_name+".gatk.vcf -o "+lib_name+".gatk.select_snp.vcf && \\\n")
	gatk_sh.write(bin_path+"/software/java1.8 -Djava.io.tmpdir="+JavaTmpDir+" -jar "+bin_path+"/software/GenomeAnalysisTK.jar -T VariantFiltration -R "+IndexDir+"/refMrna.fa -V "+lib_name+".gatk.select_snp.vcf "+' canshu' +" -o "+lib_name+".raw.snp.vcf && \\\n")
	gatk_sh.write("awk '(/^#/ || \$7 == \"PASS\")' "+lib_name+".raw.snp.vcf >"+lib_name+".snp.vcf && \\\n")
	gatk_sh.close()

	depth_sh = open (snp_sh_path + '/'+ lib_name+'/depth_'+lib_name+'.sh','w')
	depth_sh.write(bin_path+"/software/samtools depth $processDir/$s/$s.addRG.Reorder.Sort.bam >$processDir/$s/$s.depth ")
	depth_sh.close()
	
	snp_basic_stat_sh=open(snp_sh_path + '/'+ lib_name+'/snp_basic_stat.sh','w')
	snp_basic_stat_sh.write("perl "+bin_path+"/SNP/snp_statistics.pl -sn\np -out && \\\n")
	snp_basic_stat_sh.close()
'''
	snp_population_stat_sh =open(,'w')
	snp_population_stat_sh.write("perl $Bin/SNP/snp_population.pl -snp " . join(",", map{$snp_info{$_}{'snp'}} @keys) . " -depth " . join(",", map{$snp_info{$_}{'depth'}} @keys) . " -outdir $processDir && \\\n")
	snp_population_stat_sh.close()
'''

###CDS###
Check_dir(outdir + '/CDSpredict_ESTscan')
Check_dir(outdir + '/CDSpredict_ESTscan/result')

cds = config.get('script_path','cds')

cds_sh_path = outdir + '/CDSpredict_ESTscan'
cds_result_path = outdir + '/CDSpredict_ESTscan/result'

Check_dir(outdir + '/CDSpredict_ESTscan/result/Blast')
Check_dir(outdir + '/CDSpredict_ESTscan/result/ESTscan/estscan')

os.system('perl '+cds +' -fa '+ denovo_fa + ' -annot ' + 'nr'+','+'swissprot'+','+'kegg'+','+'cog'+' -shdir '+ cds_sh_path +' -outdir ' + cds_result_path)

###GeneExpression###
Check_dir(outdir + '/05.GeneExp')
Check_dir(outdir + '/05.GeneExp/result')
Check_dir(outdir + '/05.GeneExp/result/bowtie-build')
Check_dir(outdir + '/05.GeneExp/result/rsem-build')
Check_dir(outdir + '/05.GeneExp/result/PCA')

geneexp_sh_path = outdir + '/05.GeneExp'
geneexp_result_path = outdir + '/05.GeneExp/result'
BowtieIndexDir = geneexp_result_path + '/bowtie-build'
RsemIndexDir = geneexp_result_path + '/rsem-build'
PCA_result_path = geneexp_result_path + '/PCA'

bowtie2_para = config.get('software_parameter','bowtie2_para')
rsem_para = config.get('software_parameter','rsem_para')


bulid_index_sh = open(geneexp_sh_path+'/bulid_index.sh','w')
bulid_index_sh.write("echo ==========start at : `date` ========== && \\\n")
bulid_index_sh.write("if [ ! -f " + BowtieIndexDir +'/refMrna.fa ];then ln -s '+ denovo_fa + ' ' +BowtieIndexDir +'/refMrna.fa;fi && \\\n')
bulid_index_sh.write("grep '^>' "+ denovo_fa +"| sed 's/>//g' | awk '{print \$1\"\\t\"\$1}' >" + geneexp_result_path +"/rsem-build/gene2tr.txt && \\\n")
bulid_index_sh.write("export LD_LIBRARY_PATH=$Bin/software/RNA_lib:\$LD_LIBRARY_PATH && \\\n")
bulid_index_sh.write("$Bin/software/bowtie2/bowtie2-build -f " + BowtieIndexDir +"/refMrna.fa "+ BowtieIndexDir +"/refMrna.fa && \\\n")
bulid_index_sh.write("$Bin/software/rsem/rsem-prepare-reference " + denovo_fa + ' ' + RsemIndexDir+ "/refMrna.fa --bowtie2 --bowtie2-path $Bin/software/bowtie2 --transcript-to-gene-map " + geneexp_result_path +"/rsem-build/gene2tr.txt && \\\n")
bulid_index_sh.write("echo ==========end at : `date` ========== && \\\n")
bulid_index_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
bulid_index_sh.write("echo Still_waters_run_deep > " + geneexp_sh_path + '/bulid_index.sh.sign\n')
bulid_index_sh.close()

for lib_name in sorted(cleandata_dict.keys()):

	Check_dir(geneexp_result_path + '/' + lib_name)
	fpkm_result = geneexp_result_path + '/fpkm'
	Check_dir(fpkm_result)
	GeneExp_sh = open(geneexp_sh_path + '/' + lib_name+'_GeneExp.sh','w')
	GeneExp_sh.write("echo ==========start at : `date` ========== && \\\n")
	GeneExp_sh.write("export LD_LIBRARY_PATH=$Bin/software/RNA_lib:\$LD_LIBRARY_PATH && \\\n")
	GeneExp_sh.write(bin_path+"/software/bowtie2/bowtie2 "+bowtie2_para+' -x '+BowtieIndexDir+'/refMrna.fa -1 '+cleandata_dict[lib_name][0]+' -2 '+cleandata_dict[lib_name][1] +'| $Bin/software/samtools view -S -b -o '+geneexp_result_path + '/' + lib_name+'/'+lib_name+'.bam - && \\\n')
	GeneExp_sh.write("perl "+bin_path+"/GeneExp/BowtieMapStat.pl -bam "+geneexp_result_path + '/' + lib_name+'/'+lib_name+'.bam -key '+geneexp_result_path + '/' + lib_name+'/'+lib_name+'.Bowtie2Gene -seqType PE -samtools $Bin/software/samtools -gene2tr '+geneexp_result_path+"rsem-build/gene2tr.txt && \\\n")
	GeneExp_sh.write(bin_path+"/software/rsem/rsem-calculate-expression --paired-end -p 8 --bam "+geneexp_result_path + '/' + lib_name+'/'+lib_name+'.bam '+RsemIndexDir+'/refMrna.fa '+geneexp_result_path + '/' + lib_name+'/'+lib_name+" && \\\n")
	GeneExp_sh.write("awk '{if(\$7!=0.00)print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$5\"\\t\"\$7}' "+geneexp_result_path + '/' + lib_name+'/'+lib_name+'.genes.results > '+geneexp_result_path + '/' + lib_name+'/'+lib_name+".gene.fpkm.xls && \\\n")
	GeneExp_sh.write("cp "+geneexp_result_path + '/' + lib_name+'/'+lib_name+'.gene.fpkm.xls '+geneexp_result_path + '/' + lib_name+'/'+lib_name+'.Bowtie2Gene.MapReadsStat.xls '+fpkm_result)
	GeneExp_sh.write("echo ==========end at : `date` ========== && \\\n")
	GeneExp_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
	GeneExp_sh.write("echo Still_waters_run_deep > " + geneexp_sh_path + '/' + lib_name+'_GeneExp.sh.sign\n')
	GeneExp_sh.close()

#PCA分析需要all_exp.lst	一般样品数在4个以上才会做PCA分析
	all_exp = open(geneexp_result_path+'/all.exp.lst','a')
	all_exp.write(lib_name+'\t'+geneexp_result_path+'/'+lib_name+'/'+lib_name+'.gene.fpkm.xls\n')
	all_exp.close()
	
'''
R1      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/R1/R1.genes.fpkm.xls
R2      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/R2/R2.genes.fpkm.xls
R3      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/R3/R3.genes.fpkm.xls
W1      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/W1/W1.genes.fpkm.xls
W2      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/W2/W2.genes.fpkm.xls
W3      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/W3/W3.genes.fpkm.xls
'''

Statistic_R_sh = open(geneexp_sh_path + '/Statistic_R.sh','w' )
Statistic_R_sh.write("echo ==========start at : `date` ========== && \\\n")
Statistic_R_sh.write("perl "+bin_path+"/GeneExp/AllGeneStat.pl "+ geneexp_result_path +' '+ geneexp_result_path +"/All.GeneExpression.FPKM.xls && \\\n")
Statistic_R_sh.write("perl "+bin_path+"/GeneExp/AllMapStat.pl "+ geneexp_result_path + ' ' + geneexp_result_path +"/MappingSummary.xls && \\\n")
Statistic_R_sh.write("perl "+bin_path+"/GeneExp/PCA_deal.pl -list "+geneexp_sh_path+'/all.exp.lst'+' -IDcolumn 1 -Exprcolumn 5 -outdir '+ PCA_result_path+ " && \\\n")
Statistic_R_sh.write("echo ==========end at : `date` ========== && \\\n")
Statistic_R_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
Statistic_R_sh.write("echo Still_waters_run_deep > " + geneexp_sh_path + '/Statistic_R.sh.sign\n')
Statistic_R_sh.close()

###GeneDiffExp###
#这里利用两个包进行差异表达计算，需要三个配置文件，一个是GroupList.txt文件，一个是SampleList.txt，还有CompareList.txt


Check_dir(outdir + '/06.GeneDiff_Allin')
Check_dir(outdir + '/06.GeneDiff_Allin/result')

dge_sh_path = outdir + '/06.GeneDiff_Allin' 
dge_result_path = outdir + '/06.GeneDiff_Allin/result'

#CompareList.txt
'''
R1,R2,R3        W1,W2,W3
'''
compare=config.get('compare_list','compare_list')
compare_list_sh = open(dge_result_path+'/CompareList.txt','w')
compare_list_sh.write(compare+'\n')
compare_list_sh.close()

#SampleList.txt
'''
R1      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/R1/R1.genes.fpkm.xls
R2      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/R2/R2.genes.fpkm.xls
R3      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/R3/R3.genes.fpkm.xls
W1      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/W1/W1.genes.fpkm.xls
W2      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/W2/W2.genes.fpkm.xls
W3      /ifshk5/PC_PA_US/PMO/F14FTSSCKF0222_PLAkqjR/xjhk7/19.yueyao/05.testRNAdenovo/FPKM/R_W_FPKM/W3/W3.genes.fpkm.xls
'''

for lib_name in sorted(cleandata_dict.keys()):
	sample_list = open(dge_result_path+'/SampleList.txt','a')
	sample_list.write(lib_name+'\t'+geneexp_result_path+'/'+lib_name+'/'+lib_name+'.gene.fpkm.xls\n')
	sample_list.close()


#GroupList.txt
'''
R       R1,R2,R3
W       W1,W2,W3
'''
for lib_name in Group_dict.keys():
	group_list = open(dge_result_path + '/GroupList.txt','a')
	group_list.write(lib_name+'\t'+','.join(Group_dict[lib_name])+'\n')
	group_list.close()

dgeseq_sh = open (dge_sh_path + '/deseq2.sh','w')
dgeseq_sh.write("echo ==========start at : `date` ========== && \\\n")
dgeseq_sh.write("export LD_LIBRARY_PATH=$Bin/software/RNA_lib:\$LD_LIBRARY_PATH && \\\n")
dgeseq_sh.write("perl "+bin_path+"/GeneDiffExp/DEseq2.pl -list "+ dge_sh_path+'/SampleList.txt -diff '+dge_sh_path+'/CompareList.txt -group '+dge_sh_path+'/GroupList.txt '+' -log2 1 -padj 0.05 -outdir '+ dge_result_path+" && \\\n")
dgeseq_sh.write("perl "+bin_path+"/GeneDiffExp/drawMA-plot.pl -indir " + dge_result_path +' -log2col 5 -exp1col 3 -exp2col 4 -udcol 8 -outdir '+ dge_result_path +" && \\\n")
dgeseq_sh.write("perl "+bin_path+"GeneDiffExp/drawVolcano-plot.pl -indir "+ dge_result_path + ' -log2col 5 -signcol 7 -udcol 8 -xlab \"log2(fold change)\" -ylab=\"-log10(Padj)\" -outdir '+dge_result_path +" && \\\n")
dgeseq_sh.write("echo ==========end at : `date` ========== && \\\n")
dgeseq_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
dgeseq_sh.write("echo Still_waters_run_deep >" + dge_sh_path + '/deseq2.sh.sign\n')
dgeseq_sh.close()


###Annotation_Blast###
Check_dir(outdir + '/Annotation_Blast')
Check_dir(outdir + '/Annotation_Blast/result/')
annotation_dbversion = config.get('software_parameter','annotation_dbversion')
annotation_db = config.get('software_parameter','annotation_db')
annotation_class = config.get('software_parameter','annotation_class')  # an(animal) or pl(plant) or fg(fungi)
annotation_sh_path = outdir + '/Annotation_Blast/'
annotation_result_path = outdir + '/Annotation_Blast/result/'
os.system('perl '+bin_path+"Annotation/Annotation_sh.pl -fa "+ denovo_fa+" -dbClass "+annotation_class + " -dbVersion " + annotation_dbversion + ' '+ annotation_db +" -shDir " + annotation_sh_path + " -outDir " + annotation_result_path )

###Venn###
Check_dir(outdir + '/Venn/')
Check_dir(outdir + '/Venn/result/')

venn_sh_path = outdir + '/Venn'
venn_result_path = outdir + '/Venn/result/'

venn_list = config.get('venn','venn')

venn_arr = venn_list.split('\t')

venn_sh = open (venn_sh_path+'/venn.sh','a')
venn_sh.write("echo ==========start at : `date` ========== && \\\n")
for venn_names in venn_arr:
	venn_file = venn_names.split(',')
	venn_file_list=[]
	for venn_i in venn_file:
		venn_i = geneexp_result_path+'/'+venn_i+'/'+venn_i+'.gene.fpkm.xls'
		venn_file_list.append(venn_i)

	venn_sh.write("perl "+bin_path+"/Annotation/venny.pl -infile " +','.join(venn_file_list) + " -name " + venn_names + " -header -color -outdir "+ venn_result_path+" -imgname "+'_'.join(venn_file)+".venn && \\\n")

venn_sh.write("echo ==========end at : `date` ========== && \\\n")
venn_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
venn_sh.write("echo Still_waters_run_deep >" + venn_sh_path + '/venn.sh.sign\n')
venn_sh.close()

###Cluster###
#需要差异表达的结果
'''
R-VS-TW.DEseq2_Method.GeneDiffExp.xls
前缀是CompareList.txt的文件
-list     input DEG file list,
format: a-vs-b.DEseq2_Method.GeneDiffExp.xls,
        b-vs-c.DEseq2_Method.GeneDiffExp.xls,
        c-vs-d.EBseq_Method.GeneDiffExp.xls,......
'''

Check_dir(outdir + '/07.Cluster_R')
Check_dir(outdir + '/07.Cluster_R/result')

cluster_r_sh_path = outdir + '/07.Cluster_R'
cluster_r_result_path = outdir + '/07.Cluster_R/result'

cluster = config.get('cluster_r','cluster_r')
cluster_arr = cluster.split(',')

cluster_list = ''
for cluster in cluster_arr:
	cluster_list = cluster_list + dge_result_path+'/'+cluster + '.DEseq2_Method.GeneDiffExp.xls,'
cluster_list = cluster_list.rstrip(',')

cluster_sh = open(cluster_r_sh_path + '/cluster_r.sh','w')
cluster_sh.write("echo ==========start at : `date` ========== && \\\n")
cluster_sh.write('perl '+bin_path+'/Cluster/Cluster_R.pl -list ' + cluster_list + " -log2col 5 -outdir " + cluster_r_result_path + ' && \\\n')
cluster_sh.write("echo ==========end at : `date` ========== && \\\n")
cluster_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
cluster_sh.write("echo Still_waters_run_deep >" + cluster_r_sh_path + '/cluster_r.sh.sign\n')
cluster_sh.close()


###GO Enrichment###
Check_dir(outdir + '/08.GO_enrichment')
Check_dir(outdir + '/08.GO_enrichment/result/')
Check_dir(outdir + '/08.GO_enrichment/tmp_file/')

go_enrichment_sh_path = outdir + '/08.GO_enrichment'
go_enrichment_result_path = outdir + '/08.GO_enrichment/result'
tmpfile_path = outdir + '/08.GO_enrichment/tmp_file/'
annotation_result_mrna = 'xxx'
spieces = 'xxx'

go_enrichment_sh = open (go_enrichment_sh_path + '/go_enrichment.sh','w')
go_enrichment_sh.write("echo ==========start at : `date` ========== && \\\n")
go_enrichment_sh.write("#export LD_LIBRARY_PATH=$Bin/software/RNA_lib:\$LD_LIBRARY_PATH && \\\n")
go_enrichment_sh.write("goclass=`cat $Bin/Annotation/db_version.txt | grep -v '^#' |awk '(\$2 == \"goclass\"){print \$3}'` && \\\n")
go_enrichment_sh.write("outdir="+go_enrichment_result_path+" && \\\n")
go_enrichment_sh.write("tmpdir="+tmpfile_path+" && \\\n")
go_enrichment_sh.write("prefix="+annotation_result_mrna+spieces+" && \\\n")
go_enrichment_sh.write("refGene="+denovo_fa+" && \\\n")
go_enrichment_sh.write("for i in `ls " + dge_result_path + "/*/*_Method.GeneDiffExpFilter.xls`\ndo\n\tkeyname=`echo \$i|sed 's/.*\\/\\(.*\\).GeneDiffExpFilter.xls/\\1/'`\n\tawk '{print \$1\"\\t\"\$5}' \$i >\$tmpdir/\$keyname.glist\n\tperl  $Bin/Annotation/drawGO.pl -list \$tmpdir/\$keyname.glist -goclass \$goclass -goprefix \$prefix -outprefix \$outdir/\$keyname\ndone && \\\n")
go_enrichment_sh.write("perl "+bin_path+"/Enrichment/go.pl -gldir \$tmpdir -sdir `dirname \$prefix` -species `basename \$prefix` -outdir \$outdir && \\\n")
go_enrichment_sh.write("grep '^>' \$refGene | sed 's/>//g' | awk '{print \$1}' >\$tmpdir/total_gene.id && \\\n")
go_enrichment_sh.write("perl "+bin_path+"/Enrichment/topGO.pl -gldir \$tmpdir -godir \$outdir -prefix \$prefix -list \$tmpdir/total_gene.id -outdir \$outdir && \\\n")
go_enrichment_sh.write("echo ==========end at : `date` ========== && \\\n")
go_enrichment_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
go_enrichment_sh.close()

###KEGG Enrichment###
Check_dir(outdir + '/09.kegg_enrichment')
Check_dir(outdir + '/09.kegg_enrichment/result/')
Check_dir(outdir + '/09.kegg_enrichment/tmp_file/')

kegg_enrichment_sh_path = outdir + '/09.kegg_enrichment'
kegg_enrichment_result_path = outdir + '/09.kegg_enrichment/result/'
kegg_tmp_file_path = outdir + '/09.kegg_enrichment/tmp_file/'

kegg_enrichment_sh = open(kegg_enrichment_sh_path + '/kegg_enrichment.sh','w')
kegg_enrichment_sh.write("echo ==========start at : `date` ========== && \\\n")
kegg_enrichment_sh.write("keggFa=`cat $Bin/Annotation/db_version.txt | grep -v '^#' |awk '(\$2 == \"$conf{'Annotation_dbClass'}\"){print \$3}'` && \\\n")
kegg_enrichment_sh.write("koMap=`cat $Bin/Annotation/db_version.txt | grep -v '^#' |awk '(\$2 == \"$conf{'Annotation_dbClass'}_komap\"){print \$3}'` && \\\n")
kegg_enrichment_sh.write("mapTitle=`cat $Bin/Annotation/db_version.txt | grep -v '^#' |awk '(\$2 == \"map_title\"){print \$3}'` && \\\n")
kegg_enrichment_sh.write("mapDir=`cat $Bin/Annotation/db_version.txt | grep -v '^#' |awk '(\$2 == \"map_dir\"){print \$3}'` && \\\n")
kegg_enrichment_sh.write("outdir="+kegg_enrichment_result_path+" && \\\n")
kegg_enrichment_sh.write("tmpdir="+kegg_tmp_file_path+" && \\\n")
kegg_enrichment_sh.write("bg="+annotation_result_mrna+spieces+".ko && \\\n")
kegg_enrichment_sh.write("for i in `ls "+ dge_result_path +"/*/*_Method.GeneDiffExpFilter.xls`\ndo\n\tkeyname=`echo \$i|sed 's/.*\\/\\(.*\\).GeneDiffExpFilter.xls/\\1/'`\n\tawk '{print \$1\"\\t\"\$5}' \$i >\$tmpdir/\$keyname.glist\n\tperl $Bin/Enrichment/getKO.pl -glist \$tmpdir/\$keyname.glist -bg \$bg -outdir \$tmpdir\n\tperl $Bin/Annotation/pathfind.pl -kegg \$keggFa -komap \$koMap -maptitle \$mapTitle -fg \$tmpdir/\$keyname.ko -bg \$bg -output \$outdir/\$keyname.path\n\t[ -d \$outdir/\$keyname\\_map ] && rm -rf \$outdir/\$keyname\\_map\n\tmkdir -p \$outdir/\$keyname\\_map\n\tperl $Bin/Enrichment/keggMap.pl -ko \$tmpdir/\$keyname.ko -diff \$tmpdir/\$keyname.glist -komap \$koMap -mapdir \$mapDir -outdir \$outdir/\$keyname\\_map\n\tperl $Bin/Annotation/drawKEGG.pl -path \$outdir/\$keyname.path -outprefix \$outdir/\$keyname -idCol 6 -level1Col 7 -level2Col 8 -geneCol 9\ndone && \\\n")
kegg_enrichment_sh.write("perl "+bin_path+"/Annotation/genPathHTML.pl -indir "+ kegg_enrichment_result_path+" && \\\n")
kegg_enrichment_sh.write("perl "+bin_path+"/Enrichment/pathway_enrichFigure.pl "+ kegg_enrichment_result_path+" && \\\n")
kegg_enrichment_sh.write("echo ==========end at : `date` ========== && \\\n")
kegg_enrichment_sh.write("echo Still_waters_run_deep 1>&2 && \\\n")
kegg_enrichment_sh.write("echo Still_waters_run_deep 1 > "+kegg_enrichment_sh_path + '/kegg_enrichment.sh.sign\n')
kegg_enrichment_sh.close()

###Report###
#print S3 "perl $Bin/ArfReport/genArfReport_cn.pl -indir $reportDir -conf $Config\n" if ($conf{'ChineseReport'} && $conf{'ChineseReport'} =~ /YES/i);
