#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os

""" 
@author:yueyao 
@file: SoftWare.py 
@time: 2018/05/18 
"""

class SoftWare():

    HISAT2="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/hisat2-2.0.4/"
    STRINGTIE="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/stringtie-1.2.3.Linux_x86_64/stringtie"
    BOWTIE2="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/bowtie2-2.2.5/"
    RSEM="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/rsem-1.2.12/"
    TRANSDECODER="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/TransDecoder-3.0.1/"
    DIAMOND="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/diamond-master/diamond"
    BLASTALL="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/blast-2.2.26/bin/blastall"
    BLASTN="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/ncbi-blast-2.4.0+/bin/blastn"
    BLASTP="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/ncbi-blast-2.4.0+/bin/blastp"
    BLASTX="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/ncbi-blast-2.4.0+/bin/blastx"
    FORMATDB="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/blast-2.2.26/bin/formatdb"
    BLAT=""
    FEATURECOUNT="/ldfssz1/ST_BIGDATA/USER/yueyao/software/subread-1.5.3-Linux-x86_64/bin/featureCounts"
    HTSEQ=""
    SOAPNUKE = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/SOAPnuke"
    FQCHECK="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/fqcheck"
    CUFFLINKS="/ldfssz1/ST_BIGDATA/USER/share_app/anaconda3/bin/"
    CPC="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/cpc-0.9-r2/"
    BWA="/ldfssz1/ST_BIGDATA/USER/yueyao/bin/miniconda2/bin/bwa"
    STAR="/ldfssz1/ST_BIGDATA/USER/yueyao/bin/miniconda2/bin/star"
    TBTOOLS=""
    TRINITY="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/trinityrnaseq-Trinity-v2.4.0/Trinity"
    OSAS=""
    RMATS="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/rMATS/RNASeq-MATS.py"
    TGICL=""
    CDHIT=""
    CORSET=""
    KOBAS=""
    KASS=""
    LEAFCUTTER=""
    SGSEQ=""
    SALMON=""
    SAILFISH=""
    KALLISTO=""
    SAMBAMBA=""
    HMMSCAN="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/hmmer/bin/hmmscan"
    SUBREAD="/ldfssz1/ST_BIGDATA/USER/yueyao/software/subread-1.5.3-Linux-x86_64/bin/"
    SUBREADINDEX = "/ldfssz1/ST_BIGDATA/USER/yueyao/software/subread-1.5.3-Linux-x86_64/bin/subreadindex"
    SAMTOOLS="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/samtools-0.1.19/samtools"
    CUTADAPTER=""
    TRIM_GALORE=""
    SOAPFUSE="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/SOAPfuse-v1.26/"
    CIRCOS="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/circos-0.66/"
    GATK="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar"
    PYTHON="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/Python-2.7.9/bin/python"
    CONVERT="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/convert"
    GETORF="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/getorf"
    SSR="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/SSR"
    JAVA="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/jre1.7.0/bin/java"
    JAVA18="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/jre1.8.0_45/bin/java"
    PICARD="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/picard-tools-1.54/"
    BEDTOOLS="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/bedtools-2.17.0/bin/bedtools"
    RSCRIPT="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/R-3.1.1/bin/Rscript"
    RSCRIPT325="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/R-3.2.5/bin/Rscript"
    BOWTIE=""
    ESTSCAN=""
    def __init__(self,hisat2=HISAT2,stringtie=STRINGTIE,bowtie=BOWTIE,bowtie2=BOWTIE2,rsem=RSEM, \
                 transdecoder=TRANSDECODER,estscan=ESTSCAN,diamond=DIAMOND,blatall=BLASTALL,blastn=BLASTN, \
                 blastp=BLASTP,blat=BLAT,featurecount=FEATURECOUNT,htseq=HTSEQ,soapnuke=SOAPNUKE,fqcheck=FQCHECK,\
                 soapfuse=SOAPFUSE,cufflinks=CUFFLINKS,cpc=CPC,bwa=BWA,star=STAR,tbtools=TBTOOLS,trinity=TRINITY, \
                 osas=OSAS,rmats=RMATS,tgicl=TGICL,cdhit=CDHIT,corset=CORSET,kobas=KOBAS,kass=KASS,leafcutter=LEAFCUTTER,\
                 sgseq=SGSEQ,salmon=SALMON,sailfish=SAILFISH,kallisto=KALLISTO,sambamba=SAMBAMBA,samtools=SAMTOOLS,\
                 subread=SUBREAD,cutadapter=CUTADAPTER,trim_galore=TRIM_GALORE,gatk=GATK,subreadbuildindex=SUBREADINDEX,python=PYTHON,
                 convert=CONVERT,getorf=GETORF,ssr=SSR,formatdb=FORMATDB,hmmscan=HMMSCAN,java=JAVA,java18=JAVA18,picard=PICARD,bedtools=BEDTOOLS,
                 rscript=RSCRIPT,circos=CIRCOS):
        self.CIRCOS=circos
        self.RSCRIPT=rscript
        self.BEDTOOLS=bedtools
        self.PICARD=picard
        self.JAVA18=java18
        self.JAVA=java
        self.HMMSCAN=hmmscan
        self.FORMATDB=formatdb
        self.HISAT2=hisat2
        self.STRINGTIE=stringtie
        self.BOWTIE=bowtie
        self.RSEM=rsem
        self.TRANSDECODER=transdecoder
        self.ESTSCAN=estscan
        self.BLASTALL=blatall
        self.BOWTIE2=bowtie2
        self.DIAMOND=diamond
        self.BLASTN=blastn
        self.BLASTP=blastp
        self.BLAT=blat
        self.FEATURECOUNT=featurecount
        self.HTSEQ=htseq
        self.SOAPNUKE=soapnuke
        self.FQCHECK=fqcheck
        self.SOAPFUSE=soapfuse
        self.CUFFLINKS=cufflinks
        self.CPC=cpc
        self.BWA=bwa
        self.STAR=star
        self.TBTOOLS=tbtools
        self.TRINITY=trinity
        self.OSAS=osas
        self.RMATS=rmats
        self.TGICL=tgicl
        self.CDHIT=cdhit
        self.CORSET=corset
        self.KOBAS=kobas
        self.KASS=kass
        self.LEAFCUTTER=leafcutter
        self.SGSEQ=sgseq
        self.SALMON=salmon
        self.SAILFISH=sailfish
        self.KALLISTO=kallisto
        self.SAMBAMBA=sambamba
        self.SAMTOOLS=samtools
        self.SUBREAD=subread
        self.CUTADAPTER=cutadapter
        self.TRIM_GALORE=trim_galore
        self.GATK=gatk
        self.SUBREADINDEX=subreadbuildindex
        self.CONVERT=convert
        self.PYTHON=python
        self.GETORF=getorf
        self.SSR=ssr