#!/usr/bin/env python
# -*- coding:utf-8 -*-

import json
import logging
import os
import re
import sys


from DataBasePath import DataBasePath

from SoftWare import SoftWare

""" 
@author:yueyao 
@file: RNAseq.py 
@time: 2018/05/03 
"""

class common(object):

    def __init__(self):
        self.outdirMain = os.path.abspath('.')
        self.PEorSE="PE"
        self.ref = ""
        self.species = {}
        self.fqList = []
        self.fqLink = {}

    def makepair(self, inputfq):
        pairArray = []
        pair = {}
        order = {}
        p = 0
        for fq in inputfq:
            fqbase = os.path.basename(fq)            ## NA12878-L3_1.fq
            aa = re.search(r'_\d\.fq',fqbase)
            if aa:
                prefix = fqbase[0:aa.start()]        ## NA12878-L3
                readtype = int(fqbase[aa.start()+1]) ## 1
                try:
                    pair[prefix][readtype] = fq
                except:
                    pair[prefix] = {readtype:fq}
                if prefix in order:
                    pass
                else:
                    order[prefix] = p
                p += 1
        for pp in sorted(order.items(),key=lambda x : x[1]):
            pairArray.append(pair[pp[0]])
        return pairArray

    def prepare_lib(self):
        rawdict = {}
        for i,j in self.fqLink.items():
            rawdict[j[1]] = [j[2], j[3]]
        return rawdict

    def getsoftware(self):
        soft = SoftWare()
        return soft

class filter(common):

    def __init__(self):
        super(filter, self).__init__()
        self.parameter = "-l 15 -q 0.2 -n 0.05 -i -Q 1 -5 0  -c 0.2 " \
                         "-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA "
        self.soapnuke = self.getsoftware().SOAPNUKE
        self.fqcheck = self.getsoftware().FQCHECK
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAseq/Filter"
        self.program=[self.soapnuke,self.fqcheck]
        self.outdir = "RNAseq/Filter_SOAPnuke"


    def makeCommand(self, inputfq):
        SampleList= self.prepare_lib()
        CleanDict={}
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)

        cmd=[]
        output=[]
        filterShell=""
        filterstatshell=""
        if self.PEorSE == "PE":
            for SampleID,SampleFqList in SampleList.items():
                    FqA=SampleFqList[0]
                    FqB=SampleFqList[1]
                    cleanFqA=SampleID+'.clean.1.fq.gz'
                    cleanFqB=SampleID+'.clean.2.fq.gz'
                    rawFqA=SampleID+'.rawdata.1.fq.gz'
                    rawFqB=SampleID+'.rawdata.2.fq.gz'
                    filterShell += "{soapnuke} filter {soapnuke_para} -1 {fq1} -2 {fq2} -o {outDir} -C {clean1} -D {clean2} " \
                                  "-R {rawdata1} -W {rawdata2}\n".format(
                        soapnuke=self.soapnuke,
                        soapnuke_para=self.parameter,
                        fq1=FqA,
                        fq2=FqB,
                        outDir=self.outdir+'/'+SampleID+'/',
                        clean1=cleanFqA,
                        clean2=cleanFqB,
                        rawdata1=rawFqA,
                        rawdata2=rawFqB
                        )
                    filterstatshell += "cd {outdir};" \
                                       "{fqcheck} -r {outdir}/{clean1} -c {sampleid}_1.fqcheck;" \
                                       "{fqcheck} -r {outdir}/{clean2} -c {sampleid}_2.fqcheck;" \
                                       "{fqcheck} -r {outdir}/{rawdata1} -c {sampleid}_1.rawdata.fqcheck;" \
                                       "{fqcheck} -r {outdir}/{rawdata2} -c {sampleid}_2.rawdata.fqcheck;" \
                                       "perl {scriptbin}/fqcheck_distribute.pl {sampleid}_1.fqcheck {sampleid}_2.fqcheck -o {sampleid}.;" \
                                       "perl {scriptbin}/soapnuke_stat.pl {outdir} >{sampleid}.filter.stat.xls;" \
                                       "perl {scriptbin}/drawPizza.pl -infile {outdir}/{sampleid}.filter.stat.xls -outdir {outdir}\n".format(
                        outdir=self.outdir + '/' + SampleID + '/',
                        fqcheck=self.fqcheck,
                        clean1=cleanFqA,
                        clean2=cleanFqB,
                        rawdata1=rawFqA,
                        rawdata2=rawFqB,
                        scriptbin=self.scriptbin,
                        sampleid=SampleID
                    )
                    clean_fq1Path = "{outDir}/{clean1}".format(outDir=self.outdir+'/'+SampleID, clean1=cleanFqA)
                    clean_fq2Path = "{outDir}/{clean2}".format(outDir=self.outdir+'/'+SampleID, clean2=cleanFqB)
                    CleanDict[SampleID] = {"clean_fq1": clean_fq1Path, "clean_fq2": clean_fq2Path}
                    os.makedirs(self.outdir+'/'+SampleID+'/', mode=0o755, exist_ok=True)
            statshell="perl {scriptbin}/filter_stat.pl -indir {outdir} -output {outdir}/FilterSummary.xls".format(scriptbin=self.scriptbin,outdir=self.outdir)
            cmd.append(filterShell)
            cmd.append(filterstatshell)
            cmd.append(statshell)
            output.append(CleanDict)
        elif self.PEorSE == "SE":
            for SampleID, SampleFqList in SampleList.items():
                    Fq=SampleFqList[0]
                    cleanFq=SampleID+'.clean.fq.gz'
                    rawFq=SampleID+'.rawdata.fq.gz'
                    filterShell += "{soapnuke} filter {soapnuke_para} -1 {fq1} -o {outDir} -C {clean1} -R {rawdata1}\n".format(
                        soapnuke=self.soapnuke,
                        soapnuke_para=self.parameter,
                        fq1=Fq,
                        outDir=self.outdir+'/'+SampleID+'/',
                        clean1=cleanFq,
                        rawdata1=rawFq
                        )
                    clean_fqPath = "{outDir}/{clean1}".format(outDir=self.outdir+'/'+SampleID, clean1=cleanFq)
                    os.makedirs(clean_fqPath, mode=0o755, exist_ok=True)
                    CleanDict[SampleID] = {"clean_fq": clean_fqPath}
            cmd.append(filterShell)
            output.append(CleanDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)

        return cmd,output

    def makedefault(self, inputfq):
        SampleList= self.prepare_lib()
        CleanDict={}

        output=[]
        if self.PEorSE == "PE":
            for SampleID,SampleFqList in SampleList.items():
                    cleanFqA=SampleID+'.clean.1.fq.gz'
                    cleanFqB=SampleID+'.clean.2.fq.gz'
                    clean_fq1Path = "{outDir}/{clean1}".format(outDir=self.outdir+'/'+SampleID, clean1=cleanFqA)
                    clean_fq2Path = "{outDir}/{clean2}".format(outDir=self.outdir+'/'+SampleID, clean2=cleanFqB)
                    CleanDict[SampleID] = {"clean_fq1": clean_fq1Path, "clean_fq2": clean_fq2Path}
            output.append(CleanDict)
        elif self.PEorSE == "SE":
            for SampleID, SampleFqList in SampleList.items():
                    cleanFq=SampleID+'.clean.fq.gz'
                    clean_fqPath = "{outDir}/{clean1}".format(outDir=self.outdir+'/'+SampleID, clean1=cleanFq)
                    CleanDict[SampleID] = {"clean_fq": clean_fqPath}
            output.append(CleanDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)

        default={
            'input':inputfq,
            'parameter':self.parameter,
            'program':self.program,
            'resource':"1G,1CPU",
            'output':output
        }
        return default

class alignment(common):

    def __init__(self):
        super(alignment, self).__init__()
        soft = self.getsoftware()

        self.parameter=" --phred64 --sensitive --no-discordant --no-mixed -I 1 -X 1000 -p 8 "

        self.hisat2=soft.HISAT2+"/hisat2"
        self.samtools=soft.SAMTOOLS
        self.program=[self.hisat2,self.samtools]

        self.outdir = "RNAseq/GenomeMapping_HISAT"

    def makeCommand(self, inputfq):
        filter_para = filter()
        filter_para.fqLink = self.fqLink
        filter_para.species=self.species
        filter_para.outdir = self.outdir.replace("GenomeMapping_HISAT", "Filter_SOAPnuke")
        CleanDataDict=inputfq
        spe=filter_para.species
        database = DataBasePath()
        database.get_config(species=spe["RNAseq"][0])
        self.ref = database.HISAT_INDEX

        cmd=[]
        output=[]
        BamDict = {}
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        hisat2shell=""

        if self.PEorSE == "PE":
            for SampleID,Cleandata in CleanDataDict.items():
                cleanFqA=Cleandata["clean_fq1"]
                cleanFqB=Cleandata["clean_fq2"]
                hisat2shell+="{hisat2} {hisat2_para} -x {genome_index} -1 {fq1} -2 {fq2} " \
                             "2>{outdir}/{sampleid}.Map2GenomeStat.xls |  " \
                             "{samtools} view -b -S -o {outdir}/{sampleid}.bam -\n".format(
                    hisat2=self.hisat2,
                    hisat2_para=self.parameter,
                    genome_index=self.ref,
                    fq1=cleanFqA,fq2=cleanFqB,samtools=self.samtools,
                    outdir=self.outdir,sampleid=SampleID
                )

                BamPath="{outDir}/{sampleid}.bam" \
                        "".format(outDir=self.outdir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
            cmd.append(hisat2shell)
            output.append(BamDict)
        elif self.PEorSE == "SE":
            for SampleID, Cleandata in CleanDataDict.items():
                cleanFq=Cleandata["clean_fq"]
                hisat2shell+="{hisat2} {hisat2_para} -x {genome_index} -1 {fq1}  " \
                             "2>{outdir}/{sampleid}.Map2GenomeStat.xls |  " \
                             "{samtools} view -b -S -o {outdir}/{sampleid}.bam\n".format(
                    hisat2=self.hisat2,
                    hisat2_para=self.parameter,
                    genome_index=self.ref,
                    fq1=cleanFq,samtools=self.samtools,
                    outdir=self.outdir, sampleid=SampleID
                )
                BamPath = "{outDir}/{sampleid}.bam" \
                          "".format(outDir=self.outdir, sampleid=SampleID)
                BamDict[SampleID] = BamPath
            cmd.append(hisat2shell)
            output.append(BamDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)
        return cmd, output

    def makedefault(self, inputfq):

        filter_para = filter()
        filter_para.fqLink=self.fqLink
        filter_para.species=self.species
        filter_para.outdir=self.outdir.replace("GenomeMapping_HISAT","Filter_SOAPnuke")
        filter_output = filter_para.makedefault(inputfq)["output"]
        CleanDataDict=filter_output[0]

        database = DataBasePath()
        database.get_config(species=filter_para.species["RNAseq"][0])
        self.ref = database.GENOME

        output=[]
        BamDict = {}

        if self.PEorSE == "PE":
            for SampleID,Cleandata in CleanDataDict.items():
                BamPath="{outDir}/{sampleid}.bam" \
                        "".format(outDir=self.outdir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
            output.append(BamDict)
        elif self.PEorSE == "SE":
            for SampleID, Cleandata in CleanDataDict.items():
                cleanFq=Cleandata["clean_fq"]
                BamPath = "{outDir}/{sampleid}.bam" \
                          "".format(outDir=self.outdir, sampleid=SampleID)
                BamDict[SampleID] = BamPath
            output.append(BamDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)

        default={
            'input':CleanDataDict,
            'parameter':self.parameter,
            'program':self.program,
            'resource':"4G,8CPU",
            'output':output
        }
        return default

class geneexp(common):

    def __init__(self):
        super(geneexp,self).__init__()

        soft = self.getsoftware()

        self.cds = ""
        self.gene2tr = ""
        self.bowtie2 = soft.BOWTIE2
        self.samtools = soft.SAMTOOLS
        self.parameter="-q --phred64 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant  -p 8 -k 200"

        self.rsempreparereference=soft.RSEM+"/rsem-prepare-reference"
        self.rsem_calculate_expression=soft.RSEM+"/rsem-calculate-expression"
        self.program=[self.bowtie2,self.samtools,self.rsempreparereference,self.rsem_calculate_expression]

        self.outdir="RNAseq/GeneExp"

    def makeCommand(self,inputfq):
        filter_para = filter()
        filter_para.species=self.species
        filter_para.fqLink=self.fqLink
        filter_para.outdir = self.outdir.replace("GeneExp", "Filter_SOAPnuke")
        CleanDataDict=inputfq
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)

        spe=filter_para.species
        database = DataBasePath()
        database.get_config(species=spe["RNAseq"][0])

        self.cds = database.CDS
        self.gene2tr = database.Gene2Tr

        cmd=[]
        output=[]
        BamDict={}
        ExpDict={}

        buildindexshell="cd {outdir};{rsem_prepare_reference} {gene_fa} {outdir}/refMrna.fa --bowtie2 --bowtie2-path " \
                        "{bowtie2} --transcript-to-gene-map {gene2tr}".format(
            outdir=self.outdir,
            rsem_prepare_reference=self.rsempreparereference,
            gene_fa=self.cds,
            bowtie2=self.bowtie2,
            gene2tr=self.gene2tr
        )

        cmd.append(buildindexshell)
        rsemshell=""
        if self.PEorSE == "PE":
            for SampleID,Cleandata in CleanDataDict.items():
                os.makedirs(self.outdir + '/' + SampleID, mode=0o755, exist_ok=True)
                cleanFqA=Cleandata["clean_fq1"]
                cleanFqB=Cleandata["clean_fq2"]
                rsemshell+="{bowtie2} {bowtie2_para} -x {gene_index} -1 {fq1} -2 {fq2} 2>{outdir}/{sampleid}.Map2GeneStat.xls | " \
                          "{samtools} view -S -b -o {outdir}/{sampleid}.bam - ; " \
                          "{rsem_calculate_expression} --paired-end -p 8 --bam {bampath} {indexdir}/refMrna.fa {outdir}/{sampleid}; " \
                          "awk '{{if($7!=0.00) print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}' " \
                          "{outdir}/{sampleid}.genes.results  > {outdir}/{sampleid}.gene.fpkm.xls; " \
                          "awk '{{if($7!=0.00) print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}' " \
                          "{outdir}/{sampleid}.isoforms.results  > {outdir}/{sampleid}.transcript.fpkm.xls\n".format(
                    bowtie2=self.bowtie2+"/bowtie2",
                    bowtie2_para=self.parameter,
                    gene_index=self.outdir+"/refMrna.fa",
                    fq1=cleanFqA,fq2=cleanFqB,
                    samtools=self.samtools,
                    bampath=self.outdir+'/'+SampleID+"/"+SampleID+'.bam',
                    rsem_calculate_expression=self.rsem_calculate_expression,
                    indexdir=self.outdir,
                    outdir=self.outdir+'/'+SampleID,sampleid=SampleID
                )
                BamPath="{outDir}/{sampleid}.bam" \
                        "".format(outDir=self.outdir+'/'+SampleID,sampleid=SampleID)
                BamDict[SampleID]=BamPath
                ExpDict[SampleID]=self.outdir+'/'+SampleID+'/'+SampleID+'.gene.fpkm.xls'
            cmd.append(rsemshell)
            output.append(BamDict)
            output.append(ExpDict)
        elif self.PEorSE == "SE":
            for SampleID,Cleandata in CleanDataDict.items():
                os.makedirs(self.outdir + '/' + SampleID, mode=0o755, exist_ok=True)
                cleanFq=Cleandata["clean_fq1"]
                rsemshell="{bowtie2} {bowtie2_para} -x {gene_index} -1 {fq1} 2>{outdir}/{sampleid}.Map2GeneStat.xls | " \
                          "{samtools} view -S -b -o {outdir}/{sampleid}.bam - ; "" \
                          {rsem_calculate_expression} --paired-end -p 8 --bam {bampath} {indexdir}/refMrna.fa {outdir}; " \
                          "awk '{{if($7!=0.00) print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}' " \
                          "{outdir}/{sampleid}.genes.results  > {outdir}/{sampleid}.gene.fpkm.xls; " \
                          "awk '{{if($7!=0.00) print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}' " \
                          "{outdir}/{sampleid}.isoforms.results  > {outdir}/{sampleid}.gene.fpkm.xls\n".format(
                    bowtie2=self.bowtie2+"/bowtie2",
                    bowtie2_para=self.parameter,
                    gene_index=self.outdir+"/refMrna.fa",
                    fq1=cleanFq,
                    samtools=self.samtools,
                    indexdir=self.outdir,
                    bampath=self.outdir+'/'+SampleID+"/"+SampleID+'.bam',
                    rsem_calculate_expression=self.rsem_calculate_expression,
                    outdir=self.outdir+"/"+SampleID,sampleid=SampleID
                )
                BamPath="{outDir}/{sampleid}/{sampleid}.bam" \
                        "".format(outDir=self.outdir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
                ExpDict[SampleID]=self.outdir+'/'+SampleID+"/"+SampleID+'.gene.fpkm.xls'
            cmd.append(rsemshell)
            output.append(BamDict)
            output.append(ExpDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)
        return cmd, output

    def makedefault(self,inputfq):
        filter_para = filter()
        filter_para.species=self.species
        filter_para.fqLink=self.fqLink
        filter_para.outdir = self.outdir.replace("GeneExp", "Filter_SOAPnuke")
        filter_output = filter_para.makedefault(inputfq)["output"]
        CleanDataDict=filter_output[0]

        database = DataBasePath()
        database.get_config(species=filter_para.species["RNAseq"][0])

        BamDict={}
        ExpDict={}
        output=[]

        if self.PEorSE == "PE":
            for SampleID,Cleandata in CleanDataDict.items():
                BamPath="{outDir}/{sampleid}.bam" \
                        "".format(outDir=self.outdir+'/'+SampleID,sampleid=SampleID)
                BamDict[SampleID]=BamPath
                ExpDict[SampleID]=self.outdir+'/'+SampleID+'/'+SampleID+'.gene.fpkm.xls'
            output.append(BamDict)
            output.append(ExpDict)
        elif self.PEorSE == "SE":
            for SampleID,Cleandata in CleanDataDict.items():
                BamPath="{outDir}/{sampleid}.bam" \
                        "".format(outDir=self.outdir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
                ExpDict[SampleID]=self.outdir+'/'+SampleID+'.gene.fpkm.xls'
            output.append(BamDict)
            output.append(ExpDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)

        default={
            'input': CleanDataDict,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,8CPU",
            'output': output
        }
        return default

class genediffexp(common):

    def __init__(self):
        super(genediffexp,self).__init__()
        self.parameter={}
        self.program="DEGseq,DEseq2,EBseq,NOIseq,PossionDis"
        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAseq/GeneDiffExp"
        self.outdir="RNAseq/GeneDiffExp_Allin"

    def makeCommand(self, inputfq):
        gxp = geneexp()
        gxp.species=self.species
        gxp.fqLink=self.fqLink
        gxp.outdir = self.outdir.replace("GeneDiffExp_Allin", "GeneExp")
        ExpDict=inputfq
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        cmd=[]
        output=[]
        degshell=""

        deg_methods = self.program.split(',')

        explist = self.outdir + "/ExpList.txt"
        with open(explist, 'w') as exp:
            for Sample, GeneExp in ExpDict.items():
                exp.write(Sample + '\t' + GeneExp + '\n')

        for type in deg_methods:
            if type == "DEGseq":
                os.makedirs(self.outdir + '/DEGseq', mode=0o755, exist_ok=True)
                diffcompare = self.outdir + "/DEGseq/diffCompare"
                with open(diffcompare, 'w') as diffcomparelist:
                    compare_groups = self.parameter["DEGseq_VS"].split(',')
                    for single_group in compare_groups:
                        control = single_group.split('&')[0]
                        treat = single_group.split('&')[1]
                        diffcomparelist.write(control + "\t" + treat + "\n")
                        output.append(self.outdir + "/DEGseq/" + control + "-VS-" + treat + ".DEGseq_Method.GeneDiffExpFilter.xls")
                CompareList = self.outdir + "/DEGseq/CompareList.txt"
                with open(CompareList, 'w') as compare_list:
                    compare_groups = self.parameter["DEGseq_VS"].split(',')
                    for single_group in compare_groups:
                        control = single_group.split('&')[0]
                        treat = single_group.split('&')[1]
                        compare_list.write(
                            control + '\t' + ExpDict[control] + '\n' + treat + '\t' + ExpDict[treat] + '\n')

                degshell += "perl {GeneDiffExpBin}/DEGseq.pl -list {CompareList} -diff {diffcom} -GeneIDColumn 1 -GeneExpColumn 4 " \
                           " -GeneLenColumn 3 -method MARS -threshold 5 -pValue 1e-3 " \
                           "-zScore 4 {deg_para} -outdir {outdir}\n".format(
                    GeneDiffExpBin="/ldfssz1/ST_BIGDATA/USER/yueyao/bin",	#modify the DEGseq.pl by yueyao
                    CompareList=CompareList,
                    deg_para=self.parameter["DEGseq_Filter"],
                    outdir=self.outdir + "/DEGseq",
                    diffcom=diffcompare
                )

            if type == "DEseq2":
                os.makedirs(self.outdir + '/DEseq2', mode=0o755, exist_ok=True)
                grouplist = self.outdir + '/DEseq2/GroupList.txt'
                compare_list = self.outdir + "/DEseq2/CompareList.txt"
                tmpdict = {}
                with open(grouplist, 'w') as group:
                    groups = self.parameter["DEseq2_Group"].split(' ')
                    for single_group in groups:
                        groupsy = single_group.split(":")[0]
                        sampleid = single_group.split(":")[1]
                        tmpdict[groupsy] = sampleid
                        group.write(groupsy + "\t" + sampleid + "\n")

                with open(compare_list, 'w') as compare:
                    compare_groups = self.parameter["DEseq2_VS"].split(',')
                    for sgroup in compare_groups:
                        sampleA, sampleB = sgroup.split("&")
                        compare.write(tmpdict[sampleA] + "\t" + tmpdict[sampleB] + "\n")
                        output.append(
                            self.outdir + "/DEseq2/" + sampleA + "-VS-" + sampleB + ".DEseq2_Method.GeneDiffExpFilter.xls")
                degshell += "perl {GeneDiffExpBin}/DEseq2.pl -list {Samplelist} -diff {CompareList} -group {Grouplist} " \
                           "{deg_para} -outdir {outdir}\n".format(
                    GeneDiffExpBin=self.scriptbin,
                    Samplelist=explist,
                    CompareList=compare_list,
                    Grouplist=grouplist,
                    deg_para=self.parameter["DEseq2_Filter"],
                    outdir=self.outdir + '/DEseq2'
                )

            if type == "EBseq":
                os.makedirs(self.outdir + '/EBseq', mode=0o755, exist_ok=True)
                grouplist = self.outdir + '/EBseq/GroupList.txt'
                compare_list = self.outdir + "/EBseq/CompareList.txt"
                tmpdict = {}
                with open(grouplist, 'w') as group:
                    groups = self.parameter["EBseq_Group"].split(' ')
                    for single_group in groups:
                        groupsy = single_group.split(":")[0]
                        sampleid = single_group.split(":")[1]
                        tmpdict[groupsy] = sampleid
                        group.write(groupsy + "\t" + sampleid + "\n")

                with open(compare_list, 'w') as compare:
                    compare_groups = self.parameter["EBseq_VS"].split(',')
                    for sgroup in compare_groups:
                        sampleA, sampleB = sgroup.split("&")
                        compare.write(tmpdict[sampleA] + "\t" + tmpdict[sampleB] + "\n")
                        output.append(
                            self.outdir + "/EBseq/" + sampleA + "-VS-" + sampleB + ".EBseq_Method.GeneDiffExpFilter.xls")
                degshell += "perl {GeneDiffExpBin}/EBseq.pl -list {Samplelist} -diff {CompareList} -group {Grouplist}" \
                           " {deg_para} -outdir {outdir}\n".format(
                    GeneDiffExpBin=self.scriptbin,
                    Samplelist=explist,
                    CompareList=compare_list,
                    Grouplist=grouplist,
                    deg_para=self.parameter["EBseq_Filter"],
                    outdir=self.outdir + '/EBseq'
                )

            if type == "NOIseq":
                os.makedirs(self.outdir + '/NOIseq', mode=0o755, exist_ok=True)
                grouplist = self.outdir + '/NOIseq/GroupList.txt'
                compare_list = self.outdir + "/NOIseq/CompareList.txt"
                tmpdict = {}
                with open(grouplist, 'w') as group:
                    groups = self.parameter["NOIseq_Group"].split(' ')
                    for single_group in groups:
                        groupsy = single_group.split(":")[0]
                        sampleid = single_group.split(":")[1]
                        tmpdict[groupsy] = sampleid
                        group.write(groupsy + "\t" + sampleid + "\n")

                with open(compare_list, 'w') as compare:
                    compare_groups = self.parameter["NOIseq_VS"].split(',')
                    for sgroup in compare_groups:
                        sampleA, sampleB = sgroup.split("&")
                        compare.write(tmpdict[sampleA] + "\t" + tmpdict[sampleB] + "\n")
                        output.append(
                            self.outdir + "/NOIseq/" + sampleA + "-VS-" + sampleB + ".NOIseq_Method.GeneDiffExpFilter.xls")
                degshell += "perl {GeneDiffExpBin}/NOIseq.pl -list {Samplelist} -diff {CompareList} -group {Grouplist}" \
                           " {deg_para} -outdir {outdir}\n".format(
                    GeneDiffExpBin=self.scriptbin,
                    Samplelist=explist,
                    CompareList=compare_list,
                    Grouplist=grouplist,
                    deg_para=self.parameter["NOIseq_Filter"],
                    outdir=self.outdir + '/NOIseq'
                )

            if type == "PossionDis":
                os.makedirs(self.outdir + '/PossionDis', mode=0o755, exist_ok=True)
                CompareList = self.outdir + "/PossionDis/CompareList.txt"
                with open(CompareList, 'w') as compare_list:
                    compare_groups = self.parameter["PossionDis_VS"].split(',')
                    for single_group in compare_groups:
                        control = single_group.split('&')[0]
                        treat = single_group.split('&')[1]
                        compare_list.write(
                            control + '\t' + ExpDict[control] + '\n' + treat + '\t' + ExpDict[treat] + '\n')
                        output.append(
                            self.outdir + "/DEseq2/" + control + "-VS-" + treat + ".PossionDis_Method.GeneDiffExpFilter.xls")
                degshell += "{GeneDiffExpBin}/PossionDis.pl -list {CompareList} {deg_para} -outdir {outdir}\n".format(
                    GeneDiffExpBin=self.scriptbin,
                    CompareList=CompareList,
                    deg_para=self.parameter["PossionDis_Filter"],
                    outdir=self.outdir + '/PossionDis'
                )
        cmd.append(degshell)
        return cmd, output

    def makedefault(self,inputfq):

        gxp = geneexp()
        gxp.fqLink=self.fqLink
        gxp.species=self.species
        gxp.outdir=self.outdir.replace("GeneDiffExp_Allin","GeneExp")
        ExpDict = gxp.makedefault(inputfq)["output"][1]
        output=[]

        # self.fqLink #[A,A01,FQ1,FQ2,B]
        for i_key,line in self.fqLink.items():
            if line[0] == line[1]:
                self.program="DEGseq,PossionDis"
                self.parameter = {
                    "DEGseq_VS": "",
                    "DEGseq_Filter": "-foldChange 2 -qValue 0.001",
                    "PossionDis_VS": "",
                    "PossionDis_Filter": "-log2 1 -fdr 0.001",
                }
            else:
                self.program="DEseq2,NOIseq,EBseq"
                self.parameter={
                    "DEseq2_Group": "",
                    "DEseq2_VS": "",
                    "DEseq2_Filter": "-log2 1 -padj 0.05",
                    "EBseq_Group": "",
                    "EBseq_VS": "",
                    "EBseq_Filter": "-log2 1 -ppee 0.05",
                    "NOIseq_Group": "",
                    "NOIseq_VS": "",
                    "NOIseq_Filter": "-log2 1 -p 0.8",
                }
        deg_methods = self.program.split(',')

        for type in deg_methods:
            if type == "DEGseq":
                self.parameter["DEGseq_VS"] =""
                for i_key,line in self.fqLink.items():
                    if len(line) == 5:
                        control_tmp =line[0]
                        Treat=line[4]
                        Treat_list = Treat.split(',')
                        for treat_tmp in Treat_list:
                            self.parameter["DEGseq_VS"] += control_tmp + "&"+ treat_tmp+","
                            output.append(
                                self.outdir + "/DEGseq/" + control_tmp + "-VS-" + treat_tmp + ".DEGseq_Method.GeneDiffExpFilter.xls")
                self.parameter["DEGseq_VS"] = self.parameter["DEGseq_VS"].rstrip(",")
            if type == "PossionDis":
                self.parameter["PossionDis_VS"]=""
                for i_key,line in self.fqLink.items():
                    if len(line) == 5:
                        control_tmp =line[0]
                        Treat=line[4]
                        Treat_list = Treat.split(',')
                        for treat_tmp in Treat_list:
                            self.parameter["PossionDis_VS"] += control_tmp + "&"+ treat_tmp+","
                            output.append(
                                self.outdir + "/PossionDis/" + control_tmp + "-VS-" + treat_tmp + ".PossionDis_Method.GeneDiffExpFilter.xls")
                self.parameter["PossionDis_VS"]= self.parameter["PossionDis_VS"].rstrip(',')
            if type == "DEseq2":
                self.parameter["DEseq2_VS"]=""
                self.parameter["DEseq2_Group"]=""
                dic_gro = {}
                for i_key,line in self.fqLink.items():
                    groid=line[0]
                    samid=line[1]
                    if groid in dic_gro.keys():
                        dic_gro[groid].append(samid)
                    else:
                        dic_gro[groid]=[]
                        dic_gro[groid].append(samid)
                for i_key,line in self.fqLink.items():
                    if len(line) == 5:
                        control_tmp =line[0]
                        Treat=line[4]
                        Treat_list = Treat.split(',')
                        for treat_tmp in Treat_list:
                            self.parameter["DEseq2_VS"] += control_tmp + "&"+ treat_tmp+","
                            output.append(
                                self.outdir + "/DEseq2/" + control_tmp + "-VS-" + treat_tmp + ".DEseq2_Method.GeneDiffExpFilter.xls")
                            con_lis=",".join(dic_gro[control_tmp])
                            contm=control_tmp+":"+con_lis
                            tr_lis=",".join(dic_gro[treat_tmp])
                            trtm=treat_tmp+":"+tr_lis
                            self.parameter["DEseq2_Group"] = contm +" " +trtm +" "
                self.parameter["DEseq2_VS"] =self.parameter["DEseq2_VS"].rstrip(',')
                self.parameter["DEseq2_Group"] =self.parameter["DEseq2_Group"].rstrip(' ')
            if type == "EBseq":
                self.parameter["EBseq_VS"]=""
                self.parameter["EBseq_Group"]=""
                dic_gro = {}
                for i_key,line in self.fqLink.items():
                    groid=line[0]
                    samid=line[1]
                    if groid in dic_gro.keys():
                        dic_gro[groid].append(samid)
                    else:
                        dic_gro[groid]=[]
                        dic_gro[groid].append(samid)

                for i_key,line in self.fqLink.items():
                    if len(line) == 5:
                        control_tmp =line[0]
                        Treat=line[4]
                        Treat_list = Treat.split(',')
                        for treat_tmp in Treat_list:
                            self.parameter["EBseq_VS"] += control_tmp + "&"+ treat_tmp+","
                            output.append(
                                self.outdir + "/EBseq/" + control_tmp + "-VS-" + treat_tmp + ".EBseq_Method.GeneDiffExpFilter.xls")
                            con_lis=",".join(dic_gro[control_tmp])
                            contm=control_tmp+":"+con_lis
                            tr_lis=",".join(dic_gro[treat_tmp])
                            trtm=treat_tmp+":"+tr_lis
                            self.parameter["EBseq_Group"] = contm +" " +trtm +" "
                self.parameter["EBseq_VS"] =self.parameter["EBseq_VS"].rstrip(',')
                self.parameter["EBseq_Group"] = self.parameter["EBseq_Group"].rstrip(' ')
            if type == "NOIseq":
                self.parameter["NOIseq_VS"]=""
                self.parameter["NOIseq_Group"]=""
                dic_gro = {}
                for i_key, line in self.fqLink.items():

                    groid = line[0]
                    samid = line[1]
                    if groid in dic_gro.keys():
                        dic_gro[groid].append(samid)
                    else:
                        dic_gro[groid] = []
                        dic_gro[groid].append(samid)

                for i_key, line in self.fqLink.items():
                    if len(line) == 5:
                        control_tmp = line[0]
                        Treat = line[4]
                        Treat_list = Treat.split(',')
                        for treat_tmp in Treat_list:
                            self.parameter["NOIseq_VS"] += control_tmp + "&" + treat_tmp + ","
                            output.append(
                                self.outdir + "/NOIseq/" + control_tmp + "-VS-" + treat_tmp + ".NOIseq_Method.GeneDiffExpFilter.xls")
                            con_lis = ",".join(dic_gro[control_tmp])
                            contm = control_tmp + ":" + con_lis
                            tr_lis = ",".join(dic_gro[treat_tmp])
                            trtm = treat_tmp + ":" + tr_lis
                            self.parameter["NOIseq_Group"] = contm + " " + trtm + " "
                self.parameter["NOIseq_VS"] = self.parameter["NOIseq_VS"].rstrip(',')
                self.parameter["NOIseq_Group"] = self.parameter["NOIseq_Group"].rstrip(' ')

        default={
            'input': ExpDict,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "1G,1CPU",
            'output': output
        }
        return  default

class goenrichment(common):
    def __init__(self):
        super(goenrichment,self).__init__()

        soft = self.getsoftware()

        self.outdir="RNAseq/GO_Hypergeometric/GO"
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAseq/Enrichment"

        self.parameter={
            "Gene2Tr":"",
            "GO_Prefix":"",
            "GO_Class":"",
            "obo":"",
            "gene2go":"",
            "accession2go":""
        }
        self.program="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAseq/Enrichment/go.pl"

    def makeCommand(self,inputfq):
        deg=genediffexp()
        deg.fqLink=self.fqLink
        deg.species=self.species
        deg.outdir = self.outdir.replace("GO_Hypergeometric/GO", "GeneDiffExp_Allin")
        GeneDiffExpFilter = inputfq

        cmd=[]
        output=[]
        GODict={}
        go_shell="export PATH=/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/perl-V5/bin:$PATH; "

        go_tmpdir=self.outdir+"/tmp_file"
        os.makedirs(go_tmpdir, mode=0o755, exist_ok=True)

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            go_shell+="awk '{{if($5>0) print $1\"\\t\"$5\"\\tup\";else print $1\"\\t\"$5\"\\tdown\"}}'" \
                     " {diff_list} >{tmpdir}/{keyname}.glist; " \
                     "perl {scriptbin}/drawGO.pl -list {tmpdir}/{keyname}.glist -goclass {goclass} " \
                     "-goprefix {prefix} -outprefix {outdir}/{keyname};".format(
                    diff_list=diff_list,
                    keyname=diff_id,
                    tmpdir=go_tmpdir,
                    prefix=self.parameter["GO_Prefix"],
                    scriptbin=self.scriptbin,
                    outdir=self.outdir,
                    goclass=self.parameter["GO_Class"]
            )
            GODict[diff_id] = [self.outdir + "/" + diff_id + "_C.txt", self.outdir + "/" + diff_id + "_F.txt",
                               self.outdir + "/" + diff_id + "_P.txt"]
        go_shell+="perl {scriptbin}/go.pl -gldir {tmpdir} -sdir `dirname {prefix}` -species `basename {prefix}` -outdir {outdir};" \
                  "perl {scriptbin}/topGO.pl -gldir {tmpdir} -godir {outdir} -prefix {prefix}  -list {gene2tr} -outdir {outdir}".format(
            tmpdir=go_tmpdir,
            prefix=self.parameter["GO_Prefix"],
            scriptbin=self.scriptbin,
            outdir=self.outdir,
            gene2tr=self.parameter["Gene2Tr"]
        )

        cmd.append(go_shell)
        output.append(GODict)
        return cmd,output

    def makedefault(self,inputfq):
        deg=genediffexp()
        deg.fqLink=self.fqLink
        deg.species=self.species
        deg.outdir = self.outdir.replace("GO_Hypergeometric/GO", "GeneDiffExp_Allin")
        GeneDiffExpFilter = deg.makedefault(inputfq)["output"]

        database = DataBasePath()
        database.get_config(species=deg.species["RNAseq"][0])

        self.parameter={
            "Gene2Tr":database.Gene2Tr,
            "GO_Prefix":database.GO_PREFIX,
            "GO_Class":database.GOCLASS,
            "obo":database.OBO,
            "gene2go":database.GENE2GO,
            "accession2go":database.ACCESSION2GO
        }


        output=[]
        GODict={}

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            GODict[diff_id]=[self.outdir+"/"+diff_id+"_C.txt",self.outdir+"/"+diff_id+"_F.txt",self.outdir+"/"+diff_id+"_P.txt"]
        output.append(GODict)

        default={
            'input': GeneDiffExpFilter,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,1CPU",
            'output': output
        }
        return default

class pathwayenrichment(common):
    def __init__(self):
        super(pathwayenrichment,self).__init__()

        self.outdir="RNAseq/Pathway_Hypergeometric/KEGG"
        self.parameter={
            "KO":"",
            "map_title":"",
            "mapdir":"",
            "pl_komap":"",
            "kegg_fa":""
        }

        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAseq/Enrichment"
        self.program="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAseq/Enrichment/pathfind.pl"

    def makeCommand(self,inputfq):
        deg=genediffexp()
        deg.species=self.species
        deg.fqLink=self.fqLink
        deg.outdir = self.outdir.replace("Pathway_Hypergeometric/KEGG", "GeneDiffExp_Allin")
        GeneDiffExpFilter = inputfq

        cmd=[]
        output=[]
        KEGGDict={}
        kegg_shell=""
        kegg_tmpdir = self.outdir + "/tmp_file"
        os.makedirs(kegg_tmpdir,mode=0o755, exist_ok=True)

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            os.makedirs(self.outdir + '/'+diff_id+'_map',mode=0o755, exist_ok=True)
            kegg_shell+="awk '{{print $1\"\\t\"$5}}' {diff_list} >{tmpdir}/{keyname}.glist; " \
                       "/usr/bin/perl {scriptbin}/getKO.pl -glist {tmpdir}/{keyname}.glist -bg {bg_ko} -outdir {tmpdir}; " \
                       "/usr/bin/perl {scriptbin}/pathfind.pl -kegg {kegg_fa} -komap {ko_map} -maptitle {map_title} -fg {tmpdir}/{keyname}.ko " \
                       "-bg {bg_ko} -output {outdir}/{keyname}.path; " \
                       "awk '{{if($5>0) print $1\"\\t\"$5\"\\tup\";else print $1\"\\t\"$5\"\\tdown\"}} ' {diff_list}" \
                       " > {tmpdir}/{keyname}.glist.temp; " \
                       "/usr/bin/perl {scriptbin}/keggMap.pl -ko {tmpdir}/{keyname}.ko -diff {tmpdir}/{keyname}.glist -komap " \
                       "{ko_map} -mapdir {mapdir} -outdir {outdir}/{keyname}_map; " \
                       "perl {scriptbin}/drawKEGG.pl -path {outdir}/{keyname}.path -outprefix {outdir}/{keyname}  " \
                       "-idCol 6 -level1Col 7 -level2Col 8 -geneCol 9 -list {tmpdir}/{keyname}.glist.temp;".format(
                diff_list=diff_list,
                tmpdir=kegg_tmpdir,
                keyname=diff_id,
                bg_ko=self.parameter["KO"],
                scriptbin=self.scriptbin,
                outdir=self.outdir,
                kegg_fa=self.parameter["kegg_fa"],
                ko_map=self.parameter["pl_komap"],
                map_title=self.parameter["map_title"],
                mapdir=self.parameter["mapdir"]
            )
            KEGGDict[diff_id] = [self.outdir + "/" + diff_id + ".path"]
        kegg_shell +="perl {scriptbin}/genPathHTML.pl -indir {outdir};" \
                     "perl {scriptbin}/pathway_enrichFigure.pl {outdir}".format(
            scriptbin=self.scriptbin,
            outdir=self.outdir
        )
        cmd.append(kegg_shell)
        output.append(KEGGDict)
        return cmd,output

    def makedefault(self,inputfq):

        deg=genediffexp()
        deg.fqLink=self.fqLink
        deg.species=self.species
        deg.outdir = self.outdir.replace("Pathway_Hypergeometric/KEGG", "GeneDiffExp_Allin")
        GeneDiffExpFilter = deg.makedefault(inputfq)["output"]

        database = DataBasePath()
        database.get_config(species=deg.species["RNAseq"][0])

        self.parameter={
            "KO":database.KO,
            "map_title":database.KEGG_MAP_TITLE,
            "mapdir":database.KEGG_MAP_DIR,
            "pl_komap":database.KEGG_KOMAP,
            "kegg_fa":database.KEGG
        }

        output=[]
        KEGGDict={}

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            KEGGDict[diff_id] = [self.outdir + "/" + diff_id + ".path"]
        output.append(KEGGDict)


        default={
            'input': GeneDiffExpFilter,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,1CPU",
            'output': output
        }
        return default

class wgcna(common):
    def __init__(self):
        super(wgcna,self).__init__()
        self.parameter = {
            "GeneFracThreshold":0.9,
            "WeightThreshold":0.9,
            "NodesNumThreshold":1000
        }
        self.program="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAseq/GeneCoExpression/bin/WGCNA.pl"
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAseq/GeneCoExpression/bin"
        self.outdir = "RNAseq/GeneCoExpression_WGCNA"

    def makeCommand(self,inputfq):
        wgcna_shell=""

        cmd=[]
        output=[]
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        fpkm_lis = inputfq[0]
        with open(self.outdir+"/tmp.list",'w') as fw:
            for i_key,i_value in fpkm_lis.items():
                fw.write(i_key+"\t"+i_value+"\n")
        if len(fpkm_lis.keys()) >=6:
            wgcna_shell+="perl {scriptbin}/getExpMatrix.pl -list {fpkm_lis} -IDcolumn 1 -Exprcolumn 5 -geneFrac {genefrac} -outdir {outdir};" \
                        "perl {scriptbin}/WGCNA.pl -expMatrix {outdir}/geneExpMatrix.xls -weight {weight} -nodesNum {nodenum} -outdir {outdir};" \
                        "perl {scriptbin}/Coexpression_netwoerk0.pl -indir {outdir} -outdir {outdir};" \
                        "perl {scriptbin}/Coexpression_netwoerk.pl -indir {outdir} -outdir {outdir};".format(
                            scriptbin=self.scriptbin,
                            fpkm_lis=self.outdir+"/tmp.list",
                            genefrac=self.parameter["GeneFracThreshold"],
                            outdir=self.outdir,
                            weight=self.parameter["WeightThreshold"],
                            nodenum=self.parameter["NodesNumThreshold"]
                        )
            output.append(self.outdir + "/CytoscapeInput-edges.txt")
            output.append(self.outdir + "/Coexpression.network1.pdf")
            output.append(self.outdir + "/Modules_heatmap.pdf")
        else:
            wgcna_shell +="echo \"Your sample nume must above 6\" >{outdir}/wgcna.log".format(outdir=self.outdir)
            output.append(self.outdir+"/wgcna.log")
        cmd.append(wgcna_shell)


        return cmd,output

    def makedefault(self, inputfq):
        geneexp_o = geneexp()
        geneexp_o.species=self.species
        geneexp_o.fqLink=self.fqLink
        geneexp_o.outdir=self.outdir.replace("GeneCoExpression_WGCNA","GeneExp/")
        Expdict=geneexp_o.makedefault(inputfq)["output"][1]


        input=[]
        input.append(Expdict)
        output=[]
        if len(Expdict.keys()) >=6:
            output.append(self.outdir + "/CytoscapeInput-edges.txt")
            output.append(self.outdir + "/Coexpression.network1.pdf")
            output.append(self.outdir + "/Modules_heatmap.pdf")
        else:
            output.append(self.outdir + "/wgcna.log")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,1CPU",
            'output': output
        }
        return default

class preresult(common):
    def __init__(self):
        super(preresult, self).__init__()
        self.parameter = ""
        self.program=""
        self.outdir = "BGI_result"
    def makeCommand(self, inputfq):
        outd = self.outdir.replace("BGI_result", "RNAseq")
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/1.CleanData", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/2.MapStat", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/3.Structure", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/4.Quantify", mode=0o755, exist_ok=True)
        cpshell=""
        cmd=[]
        output=[]
        cpshell +="cp {filteroutdir}/*/*.filter.stat.xls {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.RawReadsClass.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.base.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.qual.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/FilterSummary.xls {outdir}/1.CleanData/;" \
                  "mkdir -p {outdir}/2.MapStat/GenomeMapping/;" \
                  "cp {alignmentoutdir}/GenomeMappingSummary.xls {outdir}/2.MapStat/GenomeMapping/;" \
                  "cp {alignmentoutdir}/*.Map2GenomeStat.xls {outdir}/2.MapStat/GenomeMapping/;" \
                  "cp {alignmentoutdir}/*/*AddRG.Reorder.Sort.bam {alignmentoutdir}/*/*AddRG.Reorder.Sort.bam.bai {outdir}/2.MapStat/GenomeMapping/" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/GeneExpression/;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/DEGList/;" \
                  "cp {geneexpoutdir}/*/*.fpkm.xls {outdir}/4.Quantify/GeneExpression/GeneExpression/;" \
                  "mkdir -p {outdir}/2.MapStat/GeneMapping/;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/GeneOntolotyEnrichment/;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/KeggPathwayEnrichmen/;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/GeneCoExpression_WGCNA/;" \
                  "cp {geneexpoutdir}/*/*.Map2GeneStat.xls {outdir}/2.MapStat/GeneMapping/;" \
                  "cp {genediffexpoutdir}/*/*GeneDiffExp*xls {genediffexpoutdir}/*/*.MA-plot.* {genediffexpoutdir}/*/*.Scatter-plot.* {genediffexpoutdir}/*/*.Volcano-plot.* {outdir}/4.Quantify/DifferentiallyExpressedGene/DEGList/;" \
                  "cp {goenrichmentoutdir}/* {outdir}/4.Quantify/DifferentiallyExpressedGene/GeneOntolotyEnrichment/;" \
                  "cp {pathwayenrichmentoutdir}/{{*.xls,*.htm,*.pdf,*.png,*.path,*map}} {outdir}/4.Quantify/DifferentiallyExpressedGene/KeggPathwayEnrichmen/;" \
                  "cp -r {wgcnaoutdir}/Cytoscape* {wgcnaoutdir}/Modules* {wgcnaoutdir}/Coexpression.network* {outdir}/4.Quantify/GeneExpression/GeneCoExpression_WGCNA/".format(
            filteroutdir=outd+"/Filter_SOAPnuke",
            alignmentoutdir=outd+"/GenomeMapping_HISAT",
            geneexpoutdir=outd+"/GeneExp",
            genediffexpoutdir=outd+"/GeneDiffExp_Allin",
            goenrichmentoutdir=outd+"/GO_Hypergeometric/GO",
            pathwayenrichmentoutdir=outd+"/Pathway_Hypergeometric/KEGG",
            wgcnaoutdir=outd+"/GeneCoExpression_WGCNA",
            outdir=self.outdir
        )
        cmd.append(cpshell)
        return cmd,output

    def makedefault(self, inputfq):
        input=[]
        output=[]
        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "0.5G,1CPU",
            'output': output
        }
        return default

class interface(common):
    def __init__(self):
        super(interface,self).__init__()
        self.step = [["filter"],["alignment"],["geneexp"],["genediffexp","wgcna"],["goenrichment","pathwayenrichment"]]
        # self.step = [["preresult"]]
        self.input = "%s/workflow.json" % (self.outdirMain)
        self.output = "%s/workflow.json" % (self.outdirMain)

    def dumpjson(self, outputfile=None):
        outputjson = self.output
        if outputfile is not None:
            outputjson = outputfile
        try:
            out = open(outputjson, mode='w')
            out.write("{\n")
            for stepL in self.step:
                for step in stepL:
                    astep = eval(step)
                    astepo = astep()
                    astepo.fqLink = self.fqLink
                    # cmd,output=astepo.makeCommand(self.fqList)
                    stepdict = json.dumps(astepo.makedefault(self.fqList))
                    out.write("\"%s\":%s,\n" % (step, stepdict))
            out.write("\"outdir\":\"%s\"\n" % (self.outdirMain))
            out.write("}\n")
            out.close()

        except IOError as e:
            raise e

    def loadjson(self, inputfile=None):
        inputjson = self.input
        if inputfile is not None:
            inputjson = inputfile
        try:
            infl = open(inputjson, mode='r')
            jsondict = json.load(infl)
        except ValueError as e:
            raise e
        return jsondict

if __name__=="__main__":
    pass
