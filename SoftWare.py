#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os

""" 
@author:yueyao 
@file: SoftWare.py 
@time: 2018/05/18 
"""

class SoftWare():

    HISAT2="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/software/hisat2-2.0.4/"
    STRINGTIE="/ifs4/BC_PUB/biosoft/pipeline/Package/StringTie/stringtie-1.0.4/stringtie"
    BOWTIE=""
    BOWTIE2="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/software/bowtie2/"
    RSEM="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/software/rsem/"
    TRANSDECODER="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/software/TransDecoder-3.0.1/"
    ESTSCAN=""
    DIAMOND="/ifs4/BC_PUB/biosoft/pipeline/Package/diamond-master/diamond"
    BLASTALL="/ifs4/BC_PUB/biosoft/pipeline/Package/blast-2.2.26/bin/blastall"
    BLASTN="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/Annotation/../software/blastn"
    BLASTP="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/Annotation/../software/blastp"
    BLASTX="/ifs4/BC_PUB/biosoft/pipeline/Package/ncbi-blast-2.4.0+/bin/blastx"
    FORMATDB="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/SSR/../software/formatdb"
    BLAT=""
    FEATURECOUNT="/ldfssz1/ST_BIGDATA/USER/yueyao/software/subread-1.5.3-Linux-x86_64/bin/featureCounts"
    HTSEQ=""
    SOAPNUKE = "/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/software/SOAPnuke"
    FQCHECK="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/software/fqcheck"
    CUFFLINKS="/ldfssz1/ST_BIGDATA/USER/share_app/anaconda3/bin/"
    CPC="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/software/cpc/"
    BWA="/ldfssz1/ST_BIGDATA/USER/yueyao/bin/miniconda2/bin/bwa"
    STAR="/ldfssz1/ST_BIGDATA/USER/yueyao/bin/miniconda2/bin/star"
    TBTOOLS=""
    TRINITY="/ldfssz1/ST_BIGDATA/USER/yueyao/software/trinityrnaseq-Trinity-v2.4.0/Trinity"
    OSAS=""
    RMATS="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2016a/GeneDiffSplice/../software/rMATS/RNASeq-MATS.py"
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
    HMMSCAN="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/CDSpredict/../software/hmmscan"
    SUBREAD="/ldfssz1/ST_BIGDATA/USER/yueyao/software/subread-1.5.3-Linux-x86_64/bin/"
    SUBREADINDEX = "/ldfssz1/ST_BIGDATA/USER/yueyao/software/subread-1.5.3-Linux-x86_64/bin/subreadindex"
    SAMTOOLS="/ifs4/BC_PUB/biosoft/pipeline/Package/samtools-0.1.19/samtools"
    CUTADAPTER=""
    TRIM_GALORE=""
    SOAPFUSE="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2016a/GeneFusion/../software/SOAPfuse/"
    CIRCOS="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2016a/GeneFusion/../software/CircosGeneFusion/"
    GATK="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/software/GenomeAnalysisTK.jar"
    PYTHON="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2016a/PPI/../software/python"
    CONVERT="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2016a/PPI/../software/convert"
    GETORF="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/software/getorf"
    SSR="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/SSR/"
    JAVA="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/SNP/../software/java"
    JAVA18="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/SNP/../software/java1.8"
    PICARD="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/software/picard/"
    BEDTOOLS="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2016a/SnpIndel/../software/bedtools"
    RSCRIPT="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2016a/SNP/../software/Rscript"
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