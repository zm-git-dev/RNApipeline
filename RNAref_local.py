#!/usr/bin/env python
# -*- coding:utf-8 -*-

import json
import logging
import os
import re
import sys
import subprocess
import math
from multiprocessing import Pool

from DataBasePath import DataBasePath

from SoftWare import SoftWare

""" 
@author:yueyao 
@file: RNAref.py 
@time: 2018/05/15
"""


class common(object):

    def __init__(self):
        self.fqList = []
        self.fqLink = {}
        self.fqList = []
        self.fqLink = {}
        self.outdirMain = os.path.abspath('.')
        self.species = {}
        self.ref = ""
        self.PEorSE = "PE"

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
        self.parameter = "-l 15 -q 0.2 -n 0.05 -i -Q 1 -5 0 -c 0.1 " \
                         "-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA "
        self.soapnuke = self.getsoftware().SOAPNUKE
        self.fqcheck = self.getsoftware().FQCHECK
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/Filter"
        self.program=[self.soapnuke,self.fqcheck]
        self.outdir = "RNAref/Filter_SOAPnuke"

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

        self.parameter="--phred64 --no-discordant --no-mixed -I 1 -X 1000 -p 8 "

        self.samtools=soft.SAMTOOLS
        self.java=soft.JAVA
        self.picard=soft.PICARD
        self.hisat2=soft.HISAT2

        self.program=[self.samtools,self.java,self.picard,self.hisat2]

        self.outdir = "RNAref/GenomeMapping_HISAT"

    def makeCommand(self, inputfq):
        CleanDataDict=inputfq[0]

        spe=self.species
        database = DataBasePath()
        database.get_config(species=spe["RNAref"][0])
        hisat_index=database.HISAT_INDEX

        cmd=[]
        output=[]
        BamDict = {}
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        hisat2shell=""

        if self.PEorSE == "PE":
            for SampleID,Cleandata in CleanDataDict.items():
                sampledir=self.outdir+"/"+SampleID
                os.makedirs(sampledir, mode=0o755, exist_ok=True)
                os.makedirs(sampledir+"/java_tmp", mode=0o755, exist_ok=True)
                cleanFqA=Cleandata["clean_fq1"]
                cleanFqB=Cleandata["clean_fq2"]
                hisat2shell+="cd {sampledir}; {hisat2}/hisat2 {hisat2_para} -x {genome_index} -1 {fq1} -2 {fq2} " \
                             "2>{outdir}/{sampleid}.Map2GenomeStat.xls |  " \
                             "{samtools} view -b -S -o {sampledir}/{sampleid}.bam - ;" \
                             "{java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {picard}/AddOrReplaceReadGroups.jar I={sampleid}.bam " \
                             "O={sampleid}.AddRG.bam RGID={sampleid} RGLB={sampleid}_library RGPL=illumina RGPU=machine RGSM={sampleid} VALIDATION_STRINGENCY=SILENT;" \
                             "{java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {picard}/ReorderSam.jar I={sampleid}.AddRG.bam O={sampleid}.AddRG.Reorder.bam R={ref} VALIDATION_STRINGENCY=SILENT;" \
                             "{java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {picard}/SortSam.jar I={sampleid}.AddRG.Reorder.bam O={sampleid}.AddRG.Reorder.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT\n" \
                             .format(
                    sampledir=sampledir,
                    hisat2=self.hisat2,
                    java=self.java,
                    picard=self.picard,
                    ref=self.ref,
                    hisat2_para=self.parameter,
                    genome_index=hisat_index,
                    fq1=cleanFqA,fq2=cleanFqB,samtools=self.samtools,
                    outdir=self.outdir,sampleid=SampleID
                )

                BamPath="{outDir}/{sampleid}.AddRG.Reorder.Sort.bam" \
                        "".format(outDir=sampledir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
            cmd.append(hisat2shell)
            output.append(BamDict)
        elif self.PEorSE == "SE":
            for SampleID, Cleandata in CleanDataDict.items():
                sampledir = self.outdir + "/" + SampleID
                os.makedirs(sampledir, mode=0o755, exist_ok=True)
                os.makedirs(sampledir + "/java_tmp", mode=0o755, exist_ok=True)
                cleanFq = Cleandata["clean_fq"]
                hisat2shell += "cd {sampledir}; {hisat2}/hisat2 {hisat2_para} -x {genome_index} -1 {fq1} " \
                               "2>{outdir}/{sampleid}.Map2GenomeStat.xls |  " \
                               "{samtools} view -b -S -o {sampledir}/{sampleid}.bam;" \
                               "{java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {picard}/AddOrReplaceReadGroups.jar I={sampleid}.bam" \
                               "O={sampleid}.AddRG.bam RGID={sampleid} RGLB={sampleid}_library RGPL=illumina RGPU=machine RGSM={sampleid} VALIDATION_STRINGENCY=SILENT;" \
                               "{java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {picard}/ReorderSam.jar I={sampleid}.AddRG.bam O={sampleid}.AddRG.Reorder.bam R={ref} VALIDATION_STRINGENCY=SILENT;" \
                               "{java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {picard}/SortSam.jar I={sampleid}.AddRG.Reorder.bam O={sampleid}.AddRG.Reorder.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT;" \
                               "rm -rf java_tmp;rm {sampleid}.bam {sampleid}.AddRG.bam {sampleid}.AddRG.Reorder.bam\n".format(
                    sampledir=sampledir,
                    hisat2=self.hisat2,
                    java=self.java,
                    picard=self.picard,
                    ref=self.ref,
                    hisat2_para=self.parameter,
                    genome_index=hisat_index,
                    fq1=cleanFq, samtools=self.samtools,
                    outdir=self.outdir, sampleid=SampleID
                )

                BamPath = "{outDir}/{sampleid}.AddRG.Reorder.Sort.bam" \
                          "".format(outDir=sampledir, sampleid=SampleID)
                BamDict[SampleID] = BamPath
            cmd.append(hisat2shell)
            output.append(BamDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)

        return cmd, output

    def makedefault(self, inputfq):
        filter_para = filter()
        filter_para.species=self.species
        filter_para.fqLink=self.fqLink
        filter_para.outdir=self.outdir.replace("GenomeMapping_HISAT","Filter_SOAPnuke")
        filter_output = filter_para.makedefault(inputfq)["output"]
        CleanDataDict=filter_output[0]

        database = DataBasePath()
        database.get_config(species=filter_para.species["RNAref"][0])
        self.ref = database.GENOME

        input=[]
        output=[]
        BamDict={}

        input.append(CleanDataDict)

        if self.PEorSE == "PE":
            for SampleID,Cleandata in CleanDataDict.items():
                sampledir=self.outdir+"/"+SampleID
                BamPath="{outDir}/{sampleid}.AddRG.Reorder.Sort.bam" \
                        "".format(outDir=sampledir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
            output.append(BamDict)
        elif self.PEorSE == "SE":
            for SampleID, Cleandata in CleanDataDict.items():
                sampledir = self.outdir + "/" + SampleID
                BamPath = "{outDir}/{sampleid}.AddRG.Reorder.Sort.bam" \
                          "".format(outDir=sampledir, sampleid=SampleID)
                BamDict[SampleID] = BamPath
            output.append(BamDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)

        default={
            'input':input,
            'parameter':self.parameter,
            'program':self.program,
            'resource':"4G,8CPU",
            'output':output
        }
        return default

class novel_tr(common):

    def __init__(self):
        super(novel_tr,self).__init__()
        self.outdir = "RNAref/NovelTr_Cuffcompare"

        self.parameter = {
            "NovelTr_Method":"Cuffcompare",
            "Annotation_Options":"--evalue 1e-5",
            "Annotation_split": 2,
            "Annotation_dbClass":"pl",
        }

        self.cufflinks=self.getsoftware().CUFFLINKS
        self.stringtie=self.getsoftware().STRINGTIE
        self.diamond=self.getsoftware().DIAMOND
        self.cpc=self.getsoftware().CPC

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/NovelTr"
        self.program = [self.diamond,self.cpc,self.cufflinks,self.stringtie]

        self.gtf = ""
        self.pep = ""

    def makeCommand(self,inputfq):

        BamDict=inputfq

        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/Prediction/", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/Prediction/cuffmerge/", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir + "/Prediction/cpc/", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir + "/Annotation/fasta", mode=0o755, exist_ok=True)

        cmd=[]
        output=[]
        reconstruct_shell = ""
        num = self.parameter["Annotation_split"]

        spe=self.species
        database = DataBasePath()
        database.get_config(species=spe["RNAref"][0])

        self.ref=database.GENOME
        self.gtf=database.mRNA_GTF
        self.cds=database.CDS
        self.pep=database.Protein_Fasta

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)

        gtflist = open(self.outdir+"/Prediction/"+"gtf.list",'w')
        for SampleID,BamPath in BamDict.items():
            reconstruct_shell += "{stringtie} {bam_path} -G {gtf_file} -o {outdir}/{sampleid}.stringtie.gtf.tmp -f 0.3 -j 3 -c 5 -g 100 -p 5;" \
                                 "{scriptbin}/CorrectGTF.pl {ref} {outdir}/{sampleid}.stringtie.gtf.tmp {outdir}/{sampleid}.stringtie.gtf\n" \
                                 "".format(
                stringtie=self.stringtie,
                bam_path=BamPath,
                gtf_file=self.gtf,
                outdir=self.outdir+"/Prediction/",
                sampleid=SampleID,
                scriptbin=self.scriptbin,
                ref=self.ref
            )
            gtflist.write(self.outdir+"/Prediction/"+SampleID+".stringtie.gtf"+"\n")

        gtflist.close()

        cuffmerge_cpc_shell="{cufflinks}/cuffmerge -g {gtf_file} -s {ref} -o {outdir}/Prediction/cuffmerge -p 5 {gtf_list};" \
                            "perl {scriptbin}/getNovelTr_gtf.pl -input {outdir}/Prediction/cuffmerge/merged.gtf -output {outdir}/Prediction/novel_transcript.gtf;" \
                            "{cufflinks}/gffread -g {ref} -w {outdir}/Prediction/novel_transcript.fa {outdir}/Prediction/novel_transcript.gtf;" \
                            "sh {cpc} {outdir}/Prediction/novel_transcript.fa {outdir}/Prediction/cpc/novel_transcript.check {outdir}/Prediction/cpc {outdir}/Prediction/cpc/novel_transcript.evidence {pep_fa};" \
                            "perl {scriptbin}/getNovelTr_coding_FA_GTF.pl -in_fa {outdir}/Prediction/novel_transcript.fa -in_gtf {outdir}/Prediction/novel_transcript.gtf -in_check {outdir}/Prediction/cpc/novel_transcript.check -out_fa " \
                            " {outdir}/Prediction/novel_coding_transcript.fa -out_gtf {outdir}/Prediction/novel_coding_transcript.gtf -out_gene2tr {outdir}/Prediction/novel_coding_transcript.gene2tr -out_summary {outdir}/Prediction/NovelTranscriptSummary.xls;" \
                            "perl {scriptbin}/fastaDeal.pl --cutf {num} --outdir {outdir}/Annotation/fasta {outdir}/Prediction/novel_coding_transcript.fa;".format(
            cufflinks=self.cufflinks,
            scriptbin=self.scriptbin,
            gtf_file=self.gtf,
            outdir=self.outdir,
            num=num,
            gtf_list=self.outdir+"/Prediction/"+"gtf.list",
            ref=self.ref,
            cpc=self.cpc+"/bin/run_predict.db.sh",
            pep_fa = self.pep
        )

        os.makedirs(self.outdir + "/Prediction/cuffmerge/", mode=0o755, exist_ok=True)

        blast_kegg_go_cmd=""
        for i in range(1, num+1):
            blast_kegg_go_cmd += "{diamond} blastx {parameter} -d {kegg_db} -q {outdir}/Annotation/fasta/novel_coding_transcript.fa.{num} " \
                         "-o {outdir}/Annotation/fasta/novel_coding_transcript.fa.{num}.blast.kegg --threads 5 --outfmt  6 " \
                         "--seg no   --max-target-seqs 5 --more-sensitive -b 0.2 --salltitles\n".format(

                diamond=self.diamond,
                parameter=self.parameter["Annotation_Options"],
                kegg_db=dbclass.KEGG,
                outdir=self.outdir,
                num=i

            )

        for i in range(1, num+1):
            blast_kegg_go_cmd += "{diamond} blastx {parameter} -d {nr_db} -q {outdir}/Annotation/fasta/novel_coding_transcript.fa.{num} " \
                         "-o {outdir}/Annotation/fasta/novel_coding_transcript.fa.{num}.blast.nr --threads 5 --outfmt 5 " \
                         "--max-target-seqs 5 --more-sensitive -b 0.5 --salltitles;" \
                         "sed -i 's/>diamond .*</>BLASTX 2.2.28+</g' {outdir}/Annotation/fasta/novel_coding_transcript.fa.{num}.blast.nr;" \
                         "perl {scriptbin}/reformatBlast_m7.pl -input {outdir}/Annotation/fasta/novel_coding_transcript.fa.{num}.blast.nr -size 100 -output {outdir}/Annotation/fasta/novel_coding_transcript.fa.{num}.blast.nr2\n".format(
                diamond=self.diamond,
                parameter=self.parameter["Annotation_Options"],
                nr_db=dbclass.NR,
                outdir=self.outdir,
                num=i,
                scriptbin=self.scriptbin

            )

        blast_process_go_ko_cmd=""

        blast_process_go_ko_cmd += "cat {outdir}/Annotation/fasta/novel_coding_transcript.fa*.blast.nr > {outdir}/Annotation/fasta/all.blast.nr1;" \
                                   "perl {scriptbin}/blast_m7_parser.pl {outdir}/Annotation/fasta/all.blast.nr1 {outdir}/Annotation/fasta/all.blast.nr.tab;" \
                                   "/usr/bin/python {scriptbin}/BLAST3GO.py --NR_tab {outdir}/Annotation/fasta/all.blast.nr.tab --Accession_go {accession_go} --output {outdir}/Annotation/novel.annot;" \
                                   "perl {scriptbin}/annot2goa.pl {obo} {outdir}/Annotation/novel.annot {outdir}/Annotation/novel;" \
                                   "perl {scriptbin}/blast_m7_m8.pl -input \"{outdir}/Annotation/fasta/*.blast.nr\" -output {outdir}/Annotation/novel.nr.m8;" \
                                   "perl {scriptbin}/getNrDesc.pl -input {outdir}/Annotation/novel.nr.m8 -rank 1 -nr {nr_db} -output {outdir}/Annotation/novel.nr.desc\n" \
                                   "cat {outdir}/Annotation/fasta/novel_coding_transcript.fa*.blast.kegg >{outdir}/Annotation/novel.kegg;" \
                                   "perl {scriptbin}blast2ko.pl -input {outdir}/Prediction/novel_coding_transcript.fa -output {outdir}/Annotation/novel.ko -blastout {outdir}/Annotation/novel.kegg -kegg {kegg_db};".format(
            outdir=self.outdir,
            scriptbin=self.scriptbin,
            obo=dbclass.OBO,
            nr_db=dbclass.NR,
            kegg_db=dbclass.KEGG,
            accession_go=dbclass.ACCESSION2GO
        )

        cmd.append(reconstruct_shell)
        cmd.append(cuffmerge_cpc_shell)
        cmd.append(blast_kegg_go_cmd)
        cmd.append(blast_process_go_ko_cmd)

        output.append(self.outdir+"/Annotation/novel.C")
        output.append(self.outdir+"/Annotation/novel.nr.desc")
        output.append(self.outdir+"/Annotation/novel.ko")
        output.append(self.outdir + "/Prediction/novel_coding_transcript.gene2tr")
        output.append(self.outdir+"/Prediction/novel_coding_transcript.gtf")
        output.append(self.outdir + "/Prediction/novel_coding_transcript.fa")

        return  cmd,output

    def makedefault(self,inputfq):

        alignment_o = alignment()
        alignment_o.fqLink = self.fqLink
        alignment_o.species=self.species
        alignment_o.outdir=self.outdir.replace("NovelTr_Cuffcompare","GenomeMapping_HISAT")
        BamDict = alignment_o.makedefault(inputfq)["output"][0]

        output=[]

        output.append(self.outdir+"/Annotation/novel.C")
        output.append(self.outdir+"/Annotation/novel.nr.desc")
        output.append(self.outdir+"/Annotation/novel.ko")
        output.append(self.outdir + "/Prediction/novel_coding_transcript.gene2tr")
        output.append(self.outdir+"/Prediction/novel_coding_transcript.gtf")
        output.append(self.outdir + "/Prediction/novel_coding_transcript.fa")


        default={
            'input':BamDict,
            'parameter':self.parameter,
            'program':self.program,
            'resource':"4G,5CPU",
            'output':output
        }
        return default

class snpindel(common):
    def __init__(self):
        super(snpindel,self).__init__()

        self.outdir="RNAref/SnpIndel_GATK"

        self.gatk=self.getsoftware().GATK
        self.java=self.getsoftware().JAVA
        self.java18=self.getsoftware().JAVA18
        self.picard=self.getsoftware().PICARD
        self.samtools=self.getsoftware().SAMTOOLS
        self.bedtools=self.getsoftware().BEDTOOLS

        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/SnpIndel"

        self.parameter={}

        self.ref=""
        self.program=[self.gatk,self.java,self.picard,self.samtools,self.bedtools]

    def makeCommand(self,inputfq):

        BamDict=inputfq[0]
        novel_gtf=inputfq[1]

        database = DataBasePath()
        database.get_config(species=self.species["RNAref"][0])
        self.ref = database.GENOME

        cmd=[]
        output=[]
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)

        gatk_shell=""
        annot_sh=""
        stat_sh=""

        os.makedirs(self.outdir+"/java_tmp_indel", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/java_tmp_snp", mode=0o755, exist_ok=True)
        gatkdic={}
        for SampleID, BamPath in BamDict.items():
            os.makedirs(self.outdir+"/"+SampleID, mode=0o755, exist_ok=True)
            tmpdir = self.outdir + "/"+SampleID+"/java_tmp"
            os.makedirs(tmpdir, mode=0o755, exist_ok=True)

            gatk_shell += "cd {outdir}/{sampleid};" \
                          "{java} -Xmx10G -Djava.io.tmpdir={tmpdir} -jar {picard}/MarkDuplicates.jar REMOVE_DUPLICATES=false I={bam_file} O={sampleid}.markDup.bam METRICS_FILE={sampleid}.markDup.metrics VALIDATION_STRINGENCY=SILENT;" \
                          "{samtools} index {sampleid}.markDup.bam;" \
                          "{java} -Xmx10G -Djava.io.tmpdir={tmpdir} -jar {gatk} -T SplitNCigarReads -R {ref} -I {sampleid}.markDup.bam -o {sampleid}.markDup.splitN.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS --downsample_to_coverage 500;" \
                          "{samtools} index {sampleid}.markDup.splitN.bam;" \
                          "{java} -Xmx10G -Djava.io.tmpdir={tmpdir} -jar {gatk} -T HaplotypeCaller -R {ref} -I {sampleid}.markDup.splitN.bam -allowPotentiallyMisencodedQuals -stand_call_conf 20.0 -stand_emit_conf 20.0 -o {sampleid}.gatk.vcf;" \
                          "{java} -Xmx10G -Djava.io.tmpdir={tmpdir} -jar {gatk} -T SelectVariants -R {ref} -selectType SNP -V {sampleid}.gatk.vcf -o {sampleid}.gatk.select_snp.vcf;" \
                          "{java} -Xmx10G -Djava.io.tmpdir={tmpdir} -jar {gatk} -T VariantFiltration -R {ref} -V {sampleid}.gatk.select_snp.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o {sampleid}.raw.snp.vcf;" \
                          "awk '(/^#/ || $7 == \"PASS\")' {sampleid}.raw.snp.vcf >{sampleid}.snp.vcf;" \
                          "{java} -Xmx10G -Djava.io.tmpdir={tmpdir} -jar {gatk} -T SelectVariants -R {ref} -selectType INDEL -V {sampleid}.gatk.vcf -o {sampleid}.gatk.select_indel.vcf;" \
                          "{java} -Xmx10G -Djava.io.tmpdir={tmpdir} -jar {gatk} -T VariantFiltration -R {ref} -V {sampleid}.gatk.select_indel.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o {sampleid}.raw.indel.vcf;" \
                          "awk '(/^#/ || $7 == \"PASS\")' {sampleid}.raw.indel.vcf >{sampleid}.indel.vcf\n".format(
                picard=self.picard,
                outdir=self.outdir,
                sampleid=SampleID,
                java=self.java,
                samtools=self.samtools,
                tmpdir=tmpdir,
                ref=self.ref,
                bam_file=BamPath,
                gatk=self.gatk
            )
            gatkdic[SampleID]=[]
            gatkdic[SampleID].append(self.outdir + "/" + SampleID + "/" + SampleID + ".snp.vcf")
            gatkdic[SampleID].append(self.outdir + "/" + SampleID + "/" +SampleID+".indel.vcf")

            bed_dir = self.outdir + "/" + SampleID + "/bed"
            snp_annot=self.outdir + "/" + SampleID + "/snp_annot"
            indel_annot=self.outdir + "/" + SampleID + "/indel_annot"

            os.makedirs(bed_dir, mode=0o755, exist_ok=True)
            os.makedirs(snp_annot, mode=0o755, exist_ok=True)
            os.makedirs(indel_annot, mode=0o755, exist_ok=True)

            annot_sh += "cd {outdir}/{sampleid};" \
                        "cat {ref_gtf} {novel_gtf} >reference.gtf;" \
                        "perl {scriptbin}/gtf2bed.pl reference.gtf bed;" \
                        "perl {scriptbin}/VCF2Bed.pl {sampleid}.snp.vcf bed/{sampleid}.snp.bed;" \
                        "{bedtools} intersect -wb -a bed/{sampleid}.snp.bed -b bed/Exon.bed >snp_annot/{sampleid}_Exon.overlap;" \
                        "{bedtools} intersect -wb -a bed/{sampleid}.snp.bed -b bed/Intron.bed >snp_annot/{sampleid}_Intron.overlap;" \
                        "{bedtools} intersect -wb -a bed/{sampleid}.snp.bed -b bed/Up2k.bed >snp_annot/{sampleid}_Up2k.overlap;" \
                        "{bedtools} intersect -wb -a bed/{sampleid}.snp.bed -b bed/Down2k.bed >snp_annot/{sampleid}_Down2k.overlap;" \
                        "perl {scriptbin}/element.pl snp_annot snp_annot {sampleid} Exon,Intron,Up2k,Down2k {sampleid}.snp.vcf;" \
                        "perl {scriptbin}/VCF2Bed.pl {sampleid}.indel.vcf bed/{sampleid}.indel.bed;" \
                        "{bedtools} intersect -wb -a bed/{sampleid}.indel.bed -b bed/Exon.bed >indel_annot/{sampleid}_Exon.overlap;" \
                        "{bedtools} intersect -wb -a bed/{sampleid}.indel.bed -b bed/Intron.bed >indel_annot/{sampleid}_Intron.overlap;" \
                        "{bedtools} intersect -wb -a bed/{sampleid}.indel.bed -b bed/Up2k.bed >indel_annot/{sampleid}_Up2k.overlap;" \
                        "{bedtools} intersect -wb -a bed/{sampleid}.indel.bed -b bed/Down2k.bed >indel_annot/{sampleid}_Down2k.overlap;" \
                        "perl {scriptbin}/element.pl indel_annot indel_annot {sampleid} Exon,Intron,Up2k,Down2k {sampleid}.indel.vcf\n".format(
                outdir = self.outdir,
                sampleid=SampleID,
                ref_gtf=database.mRNA_GTF,
                novel_gtf=novel_gtf,
                bedtools=self.bedtools,
                scriptbin=self.scriptbin
            )
        snp_file=""
        v_snp_file=""
        v_indel_file=""
        snp_raw_file=""
        mark_bam_file=""
        for SampleID, BamPath in BamDict.items():
            snp_file +="{outdir}/{sampleid}/{sampleid}.snp.vcf,".format(outdir=self.outdir,sampleid=SampleID)

            v_snp_file += " --variant:{sampleid} {outdir}/{sampleid}/{sampleid}.snp.vcf ".format(
                sampleid=SampleID,outdir=self.outdir
            )
            v_indel_file += " --variant:{sampleid} {outdir}/{sampleid}/{sampleid}.indel.vcf ".format(
                sampleid=SampleID,outdir=self.outdir
            )
            snp_raw_file +="{outdir}/{sampleid}/{sampleid}.raw.snp.vcf,".format(
                outdir=self.outdir,
                sampleid=SampleID
            )
            mark_bam_file +="{outdir}/{sampleid}/{sampleid}.markDup.bam,".format(
                outdir=self.outdir,
                sampleid=SampleID
            )
        snp_file = snp_file.rstrip(",")
        snp_raw_file = snp_raw_file.rstrip(',')
        mark_bam_file = mark_bam_file.rstrip(',')

        stat_sh +="cd {outdir};" \
                  "perl {scriptbin}/snp_statistics.pl -snp {snp_file} -outdir {outdir};" \
                  "{java} -Xmx2G -Djava.io.tmpdir={outdir}/java_tmp_snp -jar {gatk} -T CombineVariants -R {ref} {v_snp_file} -o {outdir}/combine.snp.vcf -genotypeMergeOptions UNIQUIFY;" \
                  "{java} -Xmx2G -Djava.io.tmpdir={outdir}/java_tmp_indel -jar {gatk} -T CombineVariants -R {ref} {v_indel_file}-o {outdir}/combine.indel.vcf -genotypeMergeOptions UNIQUIFY;" \
                  "{scriptbin}/perl {scriptbin}/snp_population.pl -snp {raw_snp_file} -bam {bam_file} -outdir {outdir};".format(
            outdir=self.outdir,
            scriptbin=self.scriptbin,
            snp_file=snp_file,
            raw_snp_file=snp_raw_file,
            v_snp_file=v_snp_file,
            v_indel_file=v_indel_file,
            java=self.java,
            gatk=self.gatk,
            ref=self.ref,
            bam_file=mark_bam_file

        )
        output.append(gatkdic)
        cmd.append(gatk_shell)
        cmd.append(annot_sh)
        cmd.append(stat_sh)

        return  cmd,output

    def makedefault(self,inputfq):
        alignment_o=alignment()
        alignment_o.fqLink = self.fqLink
        alignment_o.species=self.species
        alignment_o.outdir=self.outdir.replace("SnpIndel_GATK","GenomeMapping_HISAT")
        BamDict = alignment_o.makedefault(inputfq)["output"][0]

        novel_o =novel_tr()
        novel_o.species = self.species
        novel_o.fqLink = self.fqLink
        novel_o.outdir=self.outdir.replace("SnpIndel_GATK","NovelTr_Cuffcompare")
        novel_gtf = novel_o.makedefault(inputfq)["output"][-2]

        output=[]
        input=[]
        input.append(BamDict)
        input.append(novel_gtf)
        gatkdic={}
        for SampleID, BamPath in BamDict.items():
            tmpdir = self.outdir + "/"+SampleID
            gatkdic[SampleID]=[]
            gatkdic[SampleID].append(tmpdir + "/" + SampleID + ".snp.vcf")
            gatkdic[SampleID].append(tmpdir+"/"+SampleID+".indel.vcf")
        output.append(gatkdic)
        default={
            'input':input,
            'parameter':self.parameter,
            'program':self.program,
            'resource':"15G,5CPU",
            'output':output
        }
        return default

class geneexp(common):
    def __init__(self):
        super(geneexp,self).__init__()
        self.outdir = "RNAref/GeneExp/"

        self.parameter="-q --phred64 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant  -p 8 -k 200"
        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/GeneExp"

        self.rsempreparereference=self.getsoftware().RSEM+"/rsem-prepare-reference"
        self.rsemcalculateexpression=self.getsoftware().RSEM+"/rsem-calculate-expression"
        self.samtools=self.getsoftware().SAMTOOLS
        self.bowtie2=self.getsoftware().BOWTIE2

        self.program=[self.rsempreparereference,self.samtools,self.bowtie2]

        self.cds=""
        self.gene2tr=""

    def makeCommand(self,inputfq):

        CleanDataDict = inputfq[0]
        novel_gene2tr= inputfq[1]
        novel_fa = inputfq[2]

        spe=self.species
        database = DataBasePath()
        database.get_config(species=spe["RNAref"][0])
        self.ref=database.GENOME
        self.cds=database.CDS
        self.gene2tr=database.Gene2Tr

        os.makedirs(self.outdir, mode=0o755, exist_ok=True)

        cmd=[]
        output=[]
        BamDict={}
        ExpDict={}
        stat_shell=""
        rsem_build=self.outdir+"/rsem-build"
        os.makedirs(rsem_build, mode=0o755, exist_ok=True)
        os.makedirs(self.outdir + '/GeneExpression', mode=0o755, exist_ok=True)

        buildindexshell="cat {cds_fa} {novel_fa} >{rsem_build}/refMrna.fa;" \
                        "cat {gene2tr} {novel_gene2tr} >{outdir}/gene2tr.txt;" \
                        "{rsem_prepare_reference} {rsem_build}/refMrna.fa {rsem_build}/refMrna.fa --bowtie2 --bowtie2-path " \
                        "{bowtie2} --transcript-to-gene-map {outdir}/gene2tr.txt;" \
                        "perl {scriptbin}/fastaDeal.pl -attr id:len {rsem_build}/refMrna.fa >{outdir}/TranscriptLength.txt;" \
                        "awk '{{print $1\"\\t1\\t\"$2}}' {outdir}/TranscriptLength.txt >{outdir}/TranscriptLength.bed".format(
            outdir=self.outdir,
            rsem_prepare_reference=self.rsempreparereference,
            rsem_build=rsem_build,
            cds_fa=self.cds,
            novel_fa=novel_fa,
            gene2tr=self.gene2tr,
            novel_gene2tr=novel_gene2tr,
            scriptbin=self.scriptbin,
            bowtie2=self.bowtie2,
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
                           "perl {scriptbin}/BowtieMapStat.pl -bam {bampath} -key {outdir}/{sampleid}.Bowtie2Gene -seqType PE -samtools {samtools} -gene2tr {gene2trdir}/gene2tr.txt;" \
                          "{rsem_calculate_expression} --paired-end -p 8 --bam {bampath} {rsem_build}/refMrna.fa {outdir}/{sampleid}; " \
                          "awk \'{{if($7!=0.00) print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}\' " \
                          "{outdir}/{sampleid}.genes.results  > {outdir}/{sampleid}.gene.fpkm.xls; " \
                          "awk \'{{if($7!=0.00) print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}\' " \
                          "{outdir}/{sampleid}.isoforms.results  > {outdir}/{sampleid}.transcript.fpkm.xls;" \
                           "cp {outdir}/{sampleid}.gene.fpkm.xls {outdir}/GeneExpression;" \
                           "cp {outdir}/{sampleid}.transcript.fpkm.xls {outdir}/GeneExpression\n".format(
                    bowtie2=self.bowtie2+"bowtie2",
                    scriptbin=self.scriptbin,
                    rsem_build=rsem_build,
                    bowtie2_para=self.parameter,
                    gene_index=rsem_build+"/refMrna.fa",
                    fq1=cleanFqA,fq2=cleanFqB,
                    samtools=self.samtools,
                    bampath=self.outdir+'/'+SampleID+'/'+SampleID+'.bam',
                    rsem_calculate_expression=self.rsemcalculateexpression,
                    outdir=self.outdir+'/'+SampleID,sampleid=SampleID,
                    gene2trdir=self.outdir
                )
                BamPath="{outDir}/{sampleid}/{sampleid}.bam" \
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
                          "{samtools} view -S -b -o {outdir}/{sampleid}.bam -; " \
                          "perl {scriptbin}/BowtieMapStat.pl -bam {bampath} -key {outdir}/{sampleid}.Bowtie2Gene -seqType PE -samtools {samtools} -gene2tr {gene2trdir}/gene2tr.txt;" \
                          "{rsem_calculate_expression} --paired-end -p 8 --bam {bampath} {rsem_build}/refMrna.fa {outdir}/{sampleid}; " \
                          "awk \'{{if($7!=0.00\)print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}\' " \
                          "{outdir}/{sampleid}.genes.results  > {outdir}/{sampleid}.gene.fpkm.xls; " \
                          "awk \'{{if($7!=0.00\)print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}\' " \
                          "{outdir}/{sampleid}.isoforms.results  > {outdir}/{sampleid}.transcript.fpkm.xls;" \
                          "cp {outdir}/{sampleid}.gene.fpkm.xls {ssoutdir}/GeneExpression;" \
                          "cp {outdir}/{sampleid}.transcript.fpkm.xls {ssoutdir}/GeneExpression\n".format(
                    bowtie2=self.bowtie2+"bowtie2",
                    scriptbin=self.scriptbin,
                    rsem_build=rsem_build,
                    bowtie2_para=self.parameter,
                    gene_index=rsem_build+"/refMrna.fa",
                    fq1=cleanFq,
                    samtools=self.samtools,
                    bampath=self.outdir+'/'+SampleID+'/'+SampleID+'.bam',
                    rsem_calculate_expression=self.rsemcalculateexpression,
                    outdir=self.outdir+'/'+SampleID,sampleid=SampleID,
                    ssoutdir=self.outdir,
                    gene2trdir=self.outdir
                )
                BamPath="{outDir}/{sampleid}/{sampleid}.bam" \
                        "".format(outDir=self.outdir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
                ExpDict[SampleID]=self.outdir+'/'+SampleID+'.gene.fpkm.xls'
            cmd.append(rsemshell)
            output.append(BamDict)
            output.append(ExpDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)
        os.makedirs(self.outdir+"/GeneExpression", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir + "/GeneExpression/CorrelationHeatmap", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir + "/GeneExpression/HclusterTree", mode=0o755, exist_ok=True)
        stat_shell+="perl {scriptbin}/AllGeneStat.pl {outdir} {outdir}/AllSamples.GeneExpression.FPKM.xls;" \
                    "perl {scriptbin}/AllTranscriptStat.pl {outdir} {outdir}/AllSamples.TranscriptExpression.FPKM.xls;" \
                    "perl {scriptbin}/ExpStat.pl -indir {outdir} -output {outdir}/GeneExpressionSummary.xls;" \
                    "perl {scriptbin}/barplot.GeneNumber.pl -summary {outdir}/GeneExpressionSummary.xls -outdir {outdir};" \
                    "perl {scriptbin}/stackedplot.SampleGene.pl -genelist {outdir}/AllSamples.GeneExpression.FPKM.xls -outdir {outdir};" \
                    "cp {outdir}/AllSamples.GeneExpression.FPKM.xls {outdir}/GeneExpression;" \
                    "cp {outdir}/AllSamples.TranscriptExpression.FPKM.xls {outdir}/GeneExpression;" \
                    "cp {outdir}/GeneExpressionSummary.xls  {outdir}/GeneExpression;" \
                    "cp {outdir}/*.png {outdir}/GeneExpression;" \
                    "cp {outdir}/*.pdf {outdir}/GeneExpression;" \
                    "perl {scriptbin}/drawCorrelationHeatmap.pl -indir {outdir}/GeneExpression -outdir {outdir}/GeneExpression/CorrelationHeatmap;" \
                    "perl {scriptbin}/drawHclustTree.pl -indir {outdir}/GeneExpression -outdir {outdir}/GeneExpression/HclusterTree;".format(
            scriptbin=self.scriptbin,
            outdir=self.outdir

        )

        cmd.append(stat_shell)
        output.append(self.outdir+"/GeneExpression/AllSamples.GeneExpression.FPKM.xls")
        output.append(self.outdir+"/GeneExpression/AllSamples.TranscriptExpression.FPKM.xls")
        return cmd, output

    def makedefault(self,inputfq):
        filter_para = filter()
        filter_para.outdir=self.outdir.replace("GeneExp","Filter_SOAPnuke")
        filter_para.species=self.species
        filter_para.fqLink=self.fqLink
        filter_output = filter_para.makedefault(inputfq)["output"]

        novel_o =novel_tr()
        novel_o.species = self.species
        novel_o.fqLink = self.fqLink
        novel_o.outdir=self.outdir.replace("GeneExp","NovelTr_Cuffcompare")
        novel_gene2tr= novel_o.makedefault(inputfq)["output"][-3]
        novel_fa = novel_o.makedefault(inputfq)["output"][-1]

        spe=self.species

        database = DataBasePath()
        database.get_config(species=spe["RNAref"][0])

        self.ref=database.GENOME
        self.cds=database.CDS
        self.gene2tr=database.Gene2Tr

        CleanDataDict=filter_output[0]

        BamDict={}
        ExpDict={}
        output=[]
        input=[]

        input.append(filter_output[0])
        input.append(novel_gene2tr)
        input.append(novel_fa)
        if self.PEorSE == "PE":
            for SampleID,Cleandata in CleanDataDict.items():
                # os.makedirs(self.outdir + '/' + SampleID, mode=0o755, exist_ok=True)
                BamPath="{outDir}/{sampleid}.bam" \
                        "".format(outDir=self.outdir+'/'+SampleID,sampleid=SampleID)
                BamDict[SampleID]=BamPath
                ExpDict[SampleID]=self.outdir+'/'+SampleID+'/'+SampleID+'.gene.fpkm.xls'
            output.append(BamDict)
            output.append(ExpDict)
        elif self.PEorSE == "SE":
            for SampleID,Cleandata in CleanDataDict.items():
                # os.makedirs(self.outdir + '/' + SampleID, mode=0o755, exist_ok=True)
                BamPath="{outDir}/{sampleid}.bam" \
                        "".format(outDir=self.outdir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
                ExpDict[SampleID]=self.outdir+'/'+SampleID+'.gene.fpkm.xls'
            output.append(BamDict)
            output.append(ExpDict)
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)
        output.append(self.outdir+"/GeneExpression/AllSamples.GeneExpression.FPKM.xls")
        output.append(self.outdir+"/GeneExpression/AllSamples.TranscriptExpression.FPKM.xls")
        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "4G,8CPU",
            'output': output
        }
        return default

class genediffexp(common):

    def __init__(self):
        super(genediffexp,self).__init__()
        self.parameter={}
        self.program="DEGseq,DEseq2,EBseq,NOIseq,PossionDis"
        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/GeneDiffExp"
        self.outdir="RNAref/GeneDiffExp_Allin"

    def makeCommand(self, inputfq):

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
                degseqconf=self.outdir+"DEGseq/DEGseq.conf"
                with open(degseqconf,'w') as degconf:
                    degconf.write("DEGseq_Filter = "+self.parameter["DEGseq_Filter"]+"\n")
                diffcompare = self.outdir + "DEGseq/diffCompare"
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

                degshell += "perl {GeneDiffExpBin}/DEGseq.pl  -list {CompareList} -diff {diffcom} -GeneIDColumn 1 -GeneExpColumn 4 " \
                           "-GeneLenColumn 3 -method MARS -threshold 5 -pValue 1e-3 " \
                           "-zScore 4 {deg_para} -outdir {outdir};" \
                            "perl {scriptbin}/drawMA-plot.pl -indir {outdir} -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawScatter-plot.pl -indir {outdir} -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawVolcano-plot.pl -Indir {outdir} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(FDR)\" -conf {conf} -outdir {outdir}\n" \
                            .format(
                    GeneDiffExpBin=self.scriptbin,
                    exp=CompareList,
                    conf=degseqconf,
                    scriptbin=self.scriptbin,
                    CompareList=CompareList,
                    deg_para=self.parameter["DEGseq_Filter"],
                    outdir=self.outdir + "DEGseq",
                    diffcom=diffcompare

                )

            if type == "DEseq2":
                os.makedirs(self.outdir + '/DEseq2', mode=0o755, exist_ok=True)
                grouplist = self.outdir + '/DEseq2/GroupList.txt'
                compare_list = self.outdir + "/DEseq2/CompareList.txt"
                tmpdict = {}
                deseq2conf=self.outdir+"/DEseq2/DEseq2.conf"
                with open(deseq2conf,'w') as degconf:
                    degconf.write("DEseq2_Filter = "+self.parameter["DEseq2_Filter"]+"\n")
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
                           "{deg_para} -outdir {outdir};" \
                            "perl {scriptbin}/drawMA-plot.pl -indir {outdir} -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawScatter-plot.pl -indir {outdir} -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawVolcano-plot.pl -Indir {outdir} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(FDR)\" -conf {conf} -outdir {outdir}\n" \
                    .format(
                    scriptbin=self.scriptbin,
                    conf=deseq2conf,
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
                degseqconf=self.outdir+"/EBseq/EBseq.conf"
                with open(degseqconf,'w') as degconf:
                    degconf.write("EBseq_Filter = "+self.parameter["EBseq_Filter"]+"\n")
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
                degshell += "{GeneDiffExpBin}/EBseq.pl -list {Samplelist} -diff {CompareList} -group {Grouplist}" \
                           " {deg_para} -outdir {outdir};" \
                            "perl {scriptbin}/drawMA-plot.pl -indir {outdir} -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawScatter-plot.pl -indir {outdir} -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawVolcano-plot.pl -Indir {outdir} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(FDR)\" -conf {conf} -outdir {outdir}\n" \
                    .format(
                    GeneDiffExpBin=self.scriptbin,
                    scriptbin=self.scriptbin,
                    conf=degseqconf,
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
                degseqconf=self.outdir+"/NOIseq/NOIseq.conf"
                with open(degseqconf,'w') as degconf:
                    degconf.write("NOIseq_Filter = "+self.parameter["NOIseq_Filter"]+"\n")
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
                degshell += "{GeneDiffExpBin}/NOIseq.pl -list {Samplelist} -diff {CompareList} -group {Grouplist} " \
                           "{deg_para} -outdir {outdir};" \
                            "perl {scriptbin}/drawMA-plot.pl -indir {outdir} -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawScatter-plot.pl -indir {outdir} -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawVolcano-plot.pl -Indir {outdir} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(FDR)\" -conf {conf} -outdir {outdir}\n" \
                    .format(
                    scriptbin=self.scriptbin,
                    conf=degseqconf,
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
                degseqconf=self.outdir+"PossionDis/PossionDis.conf"
                with open(degseqconf,'w') as degconf:
                    degconf.write("PossionDis_Filter = "+self.parameter["PossionDis_Filter"]+"\n")
                with open(CompareList, 'w') as compare_list:
                    compare_groups = self.parameter["PossionDis_VS"].split(',')
                    for single_group in compare_groups:
                        control = single_group.split('&')[0]
                        treat = single_group.split('&')[1]
                        compare_list.write(
                            control + '\t' + ExpDict[control] + '\n' + treat + '\t' + ExpDict[treat] + '\n')
                        output.append(
                            self.outdir + "/DEseq2/" + control + "-VS-" + treat + ".PossionDis_Method.GeneDiffExpFilter.xls")
                degshell += "{GeneDiffExpBin}/PossionDis.pl  -list {CompareList}  {deg_para}  -outdir {outdir};" \
                            "perl {scriptbin}/drawMA-plot.pl -indir {outdir} -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawScatter-plot.pl -indir {outdir} -exp1col 3 -exp2col 4 -udcol 7 -conf {conf} -outdir {outdir};" \
                            "perl {scriptbin}/drawVolcano-plot.pl -Indir {outdir} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(FDR)\" -conf {conf} -outdir {outdir}\n" \
                    .format(
                    GeneDiffExpBin=self.scriptbin,
                    scriptbin=self.scriptbin,
                    conf=degseqconf,
                    CompareList=CompareList,
                    deg_para=self.parameter["PossionDis_Filter"],
                    outdir=self.outdir + 'PossionDis'
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

class genediffsplice(common):
    def __init__(self):
        super(genediffsplice,self).__init__()
        self.parameter={
            "rMATS_Options":"-analysis U -t paired -a 8",
            "rMATS_Group":"",
            "rMATS_Filter":"-c 0.5",
            "rMATS_VS":""
        }
        self.program=self.getsoftware().RMATS
        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/GeneDiffSplice"
        self.outdir="RNAref/GeneDiffSplice"

    def makeCommand(self,inputfq):
        cmd=[]
        output=[]
        BamDict = inputfq[0]
        rmats_shell=""
        dsg_wego_shell=""
        novel_gtf = inputfq[1]
        novel_g2t = inputfq[2]
        novel_go=inputfq[3].split('.')[0]

        database = DataBasePath()
        database.get_config(species=self.species["RNAref"][0])
        spe_gtf = database.mRNA_GTF
        species_go = database.GO_PREFIX

        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/DifferentiallySplicingGene", mode=0o755, exist_ok=True)

        vs_lis = self.parameter["rMATS_VS"].split(',')
        group_lis=self.parameter["rMATS_Group"].split(' ')
        dic_gro={}
        for gr in group_lis:
            groid = gr.split(':')[0]
            grolis=gr.split(':')[1].split(',')
            dic_gro[groid]=grolis

        for vs in vs_lis:
            control_bam=""
            treat_bam=""
            control = vs.split('&')[0]
            treat = vs.split('&')[1]
            tmpdir = control+"-VS-"+treat
            if treat in dic_gro.keys():
                for repid in dic_gro[treat]:
                    control_bam +=BamDict[repid]+","
                    treat_bam += BamDict[repid]+","
            else:
                control_bam += BamDict[control]+","
                treat_bam += BamDict[treat]+","

            os.makedirs(self.outdir+"/"+tmpdir, mode=0o755, exist_ok=True)
            control_bam=control_bam.rstrip(',')
            treat_bam=treat_bam.rstrip(',')

            rmats_shell +="export PATH=/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/GeneDiffSplice/../software/Python/bin:$PATH;"
            rmats_shell +="export PATH=/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/GeneDiffSplice/../software:$PATH;"
            rmats_shell +="export  LD_LIBRARY_PATH=/share/app/gcc-4.9.3/lib64:$LD_LIBRARY_PATH;"
            rmats_shell +="cat {novel_gtf} {spe_gtf} >{outdir}/{tmpdir}/annot.gtf;" \
                          "python {rmats} -b1 {bam1} -b2 {bam2} -gtf {outdir}/{tmpdir}/annot.gtf {parameter} -len 90 {para} -o {outdir}/{tmpdir}"\
                .format(
                novel_gtf=novel_gtf,
                para=self.parameter["rMATS_Filter"],
                spe_gtf=spe_gtf,
                outdir=self.outdir,
                tmpdir=tmpdir,
                rmats=self.program,
                bam1=control_bam,
                bam2=treat_bam,
                parameter=self.parameter["rMATS_Options"]
            )

            dsg_wego_shell +="perl {scriptbin}/parse_rmats_result.pl -indir {outdir}/{tmpdir} -gtf {outdir}/{tmpdir}/annot.gtf  {para} -outdir {result_dir};" \
                             "perl {scriptbin}/merge_go.pl {novel_g2t} {novel_go},{species_go} {outdir}/{tmpdir}/annot/annot;" \
                             "perl {scriptbin}/drawGO.pl -indir {result_dir}/{tmpdir} -diff {tmpdir}.A3SS,{tmpdir}.A5SS,{tmpdir}.MXE,{tmpdir}.RI,{tmpdir}.SE -go" \
                             " {outdir}/{tmpdir}/annot/annot -goclass {goclass} -outdir {result_dir}/{tmpdir} -name {tmpdir};".format(
                scriptbin=self.scriptbin,
                outdir=self.outdir,
                result_dir=self.outdir+"DifferentiallySplicingGene",
                para=self.parameter["rMATS_Filter"],
                tmpdir=tmpdir,
                goclass=database.GOCLASS,
                novel_g2t=novel_g2t,
                species_go=species_go,
                novel_go=novel_go
            )

        asstat_shell="perl {scriptbin}/AS_statistics.pl -indir {result_dir} -outdir {result_dir}".format(
            scriptbin=self.scriptbin,
            result_dir=self.outdir + "DifferentiallySplicingGene"
        )
        cmd.append(rmats_shell)
        cmd.append(dsg_wego_shell)
        cmd.append(asstat_shell)
        output.append(self.outdir + "DifferentiallySplicingGene/AsSummary.xls")
        return cmd,output

    def makedefault(self,inputfq):
        alignment_o = alignment()
        alignment_o.fqLink=self.fqLink
        alignment_o.outdir=self.outdir.replace("GeneDiffSplice","GenomeMapping_HISAT")
        alignment_o.species=self.species
        BamDict = alignment_o.makedefault(inputfq)["output"][0]

        novel_o = novel_tr()
        novel_o.outdir=self.outdir.replace("GeneDiffSplice","NovelTr_Cuffcompare")
        novel_o.fqLink=self.fqLink
        novel_o.species=self.species
        novel_gtf = novel_o.makedefault(inputfq)["output"][-2]
        novel_g2t = novel_o.makedefault(inputfq)["output"][-3]
        novel_go=novel_o.makedefault(inputfq)["output"][0]

        spe=novel_o.species
        database = DataBasePath()
        database.get_config(species=spe["RNAref"][0])

        output=[]
        input=[]
        input.append(BamDict)
        input.append(novel_gtf)
        input.append(novel_g2t)
        input.append(novel_go)
        output.append(self.outdir + "DifferentiallySplicingGene/AsSummary.xls")
        dic_gro={}
        for i_key,line in self.fqLink.items():
            if line[0] == line[1]:
                if len (line) == 5:
                    control_tmp = line[0]
                    Treat = line[4]
                    Treat_list = Treat.split(',')
                    for treat_tmp in Treat_list:
                        self.parameter["rMATS_VS"] += control_tmp + "&" + treat_tmp + ","
                self.parameter["rMATS_VS"]=self.parameter["rMATS_VS"].rstrip(',')
            else:
                groid = line[0]
                samid = line[1]
                if groid in dic_gro.keys():
                    dic_gro[groid].append(samid)
                else:
                    dic_gro[groid] = []
                    dic_gro[groid].append(samid)

        for i_key, line in self.fqLink.items():
            if line[0] != line[1]:
                if len(line) == 5:
                    control_tmp = line[0]
                    Treat = line[4]
                    Treat_list = Treat.split(',')
                    for treat_tmp in Treat_list:
                        self.parameter["rMATS_VS"] += control_tmp + "&" + treat_tmp + ","
                        output.append(self.outdir+"/"+control_tmp+"-VS-"+treat_tmp)
                        con_lis = ",".join(dic_gro[control_tmp])
                        contm = control_tmp + ":" + con_lis
                        tr_lis = ",".join(dic_gro[treat_tmp])
                        trtm = treat_tmp + ":" + tr_lis
                        self.parameter["rMATS_Group"] = contm + " " + trtm + " "
                    self.parameter["rMATS_VS"] = self.parameter["rMATS_VS"].rstrip(',')
                    self.parameter["rMATS_Group"] = self.parameter["rMATS_Group"].rstrip(' ')


        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "1G,1CPU",
            'output': output
        }

        return  default

class tf(common):
    def __init__(self):
        super(tf, self).__init__()
        self.outdir = "RNAref/TFpredict"
        self.diamond=self.getsoftware().DIAMOND
        self.getorf=self.getsoftware().GETORF
        self.rscript=self.getsoftware().RSCRIPT
        self.convert=self.getsoftware().CONVERT
        self.program=[self.diamond,self.getorf]

        self.parameter = {
            "Annotation_dbClass": "pl"
        }

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/TF"
        self.cds = ""
        self.gene2tr = ""

    def makeCommand(self, inputfq):

        novel_fa=inputfq[0]
        novel_gene2tr=inputfq[1]
        all_trans_fpkm=inputfq[2]
        all=inputfq[3]

        genediffexp_o = genediffexp()
        genediffexp_o.species=self.species
        genediffexp_o.fqLink=self.fqLink
        genediffexp_o.outdir=self.outdir.replace("TFpredict","GeneDiffExp_Allin")
        diffdir = genediffexp_o.outdir

        os.makedirs(self.outdir, mode=0o775, exist_ok=True)

        self.parameter["Annotation_dbClass"] = self.species["RNAref"][1]

        spe=self.species
        database = DataBasePath(dbclass=self.parameter["Annotation_dbClass"])
        database.get_config(species=spe["RNAref"][0])

        self.cds=database.CDS
        self.gene2tr=database.Gene2Tr

        tf_sh = ""
        cmd = []
        output = []

        if self.parameter["Annotation_dbClass"] == "an":
            db = database.TF
            tf_sh += "cat {cds} {novel_fa} >{outdir}/reference.fa;" \
                     "cat {gene2tr} {novel_gene2tr} >{outdir}/reference.gene2tr;" \
                     "{getorf} -minsize 150 -sequence {outdir}/reference.fa -outseq {outdir}/reference.orf;" \
                     "perl {scriptbin}/select_orf.pl -seq {outdir}/reference.fa -input {outdir}/reference.orf -output {outdir}/reference.pep;" \
                     "{diamond} blastp --evalue 1e-5 --threads 4 --outfmt 6 -d {db} -q {outdir}/reference.pep -o {outdir}/outfmt6 --seg no --max-target-seqs 1 --more-sensitive -b 0.5 --salltitles;" \
                     ">{outdir}/result && for i in `ls /ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/TF/Animal-TF-tab/*xls`;do " \
                     "awk -F'\\t' '{{print $2\"\\t\"$(NF-1)\"\\t\"$NF}}' $i;done|sort|uniq|awk -F'\\t' 'NR==FNR{{a[$1]=$0;next}}{{split($2,array,/:/);if (array[1] in a) print $1\"\\t\"a[array[1]];else if(array[2] in a)print $1\"\\t\"a[array[2]]}}' - {outdir}/outfmt6 >>{outdir}/result;" \
                     "sed -r 's/_[0-9]+\\t/\\t/g' {outdir}/result |sort|uniq >{outdir}/Animal_TF_result.xls;" \
                     "perl {scriptbin}/get_tf_an.pl {outdir}/Animal_TF_result.xls {deg_dir} {outdir} {outdir}/reference.gene2tr {AllSamples_TranscriptExpression};" \
                     "perl {scriptbin}/TF_heatmap.R.pl {outdir}/TFheatmap.xls {outdir};" \
                     "{scriptbin}/../software/Rscript {scriptbin}/draw_bar.R {outdir}/All_DEG.TFGene.xls {outdir}/All_DEG.TFGene.pdf;" \
                     "{scriptbin}/../software/convert -density 300 -resize 30% {outdir}/All_DEG.TFGene.pdf {outdir}/All_DEG.TFGene.png;".format(
                cds=self.cds,
                novel_fa=novel_fa,
                outdir=self.outdir,
                gene2tr=self.gene2tr,
                novel_gene2tr=novel_gene2tr,
                scriptbin=self.scriptbin,
                AllSamples_TranscriptExpression=all_trans_fpkm,
                db=db,
                getorf=self.getorf,
                diamond=self.program[0],
                deg_dir=diffdir
            )
        elif self.parameter["Annotation_dbClass"] == "pl":
            tf_sh += "export LD_LIBRARY_PATH={scriptbin}/../software/RNA_lib:$LD_LIBRARY_PATH;" \
                     "cat {cds} {novel_fa} >{outdir}/reference.fa;" \
                     "cat {gene2tr} {novel_gene2tr} >{outdir}/reference.gene2tr;" \
                     "{getorf} -minsize 150 -sequence {outdir}/reference.fa -outseq {outdir}/reference.orf;" \
                     "perl {scriptbin}/select_orf.pl -seq {outdir}/reference.fa -input {outdir}/reference.orf -output {outdir}/reference.pep;" \
                     "perl {scriptbin}/TFCodingGene_Predict.pl {outdir}/reference.pep All {outdir};" \
                     "perl {scriptbin}/get_tf_pl.pl {outdir}/Plant_TF_result.xls {deg_dir} {outdir} {outdir}/reference.gene2tr {AllSamples_TranscriptExpression};" \
                     "perl {scriptbin}/TF_heatmap.R.pl {outdir}/TFheatmap.xls {outdir};" \
                     "{Rscript} {scriptbin}/draw_bar.R {outdir}/All_DEG.TFGene.xls {outdir}/All_DEG.TFGene.pdf;" \
                     "{convert} -density 300 -resize 30% {outdir}/All_DEG.TFGene.pdf {outdir}/All_DEG.TFGene.png;".format(
                cds=self.cds,
                novel_fa=novel_fa,
                outdir=self.outdir,
                gene2tr=self.gene2tr,
                novel_gene2tr=novel_gene2tr,
                scriptbin=self.scriptbin,
                AllSamples_TranscriptExpression=all_trans_fpkm,
                Rscript=self.rscript,
                convert=self.convert,
                getorf=self.getorf,
                diamond=self.program,
                deg_dir=diffdir
            )
        elif self.parameter["Annotation_dbClass"] == "fg":
            pass
        else:
            sys.exit(1)
        cmd.append(tf_sh)
        output.append(self.outdir + "/All_DEG.TFGene.xls")
        return cmd,output

    def makedefault(self, inputfq):
        novel_o = novel_tr()
        novel_o.species=self.species
        novel_o.fqLink=self.fqLink
        novel_o.outdir=self.outdir.replace("TFpredict","NovelTr_Cuffcompare")
        novel_fa = novel_o.makedefault(inputfq)["output"][-1]
        novel_gene2tr = novel_o.makedefault(inputfq)["output"][-3]

        gxp = geneexp()
        gxp.outdir=self.outdir.replace("TFpredict","GeneExp")
        gxp.fqLink=self.fqLink
        gxp.species=self.species
        all_trans_fpkm = gxp.makedefault(inputfq)["output"][-1]

        genediff_o = genediffexp()
        genediff_o.species=self.species
        genediff_o.outdir=self.outdir.replace("TFpredict","GeneDiffExp_Allin")
        genediff_o.fqLink=self.fqLink
        tmp = genediff_o.makedefault(inputfq)["output"][-1]

        input=[]
        output=[]
        output.append(self.outdir+"/All_DEG.TFGene.xls")
        self.parameter["Annotation_dbClass"]=self.species["RNAref"][1]

        input.append(novel_fa)
        input.append(novel_gene2tr)
        input.append(all_trans_fpkm)
        input.append(tmp)

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,1CPU",
            'output': output
        }

        return  default

class goenrichment(common):
    def __init__(self):
        super(goenrichment,self).__init__()
        self.outdir = "RNAref/GO_Hypergeometric/GO"
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/Enrichment"

        self.parameter = {
            "Annotation_dbClass": "pl"
        }

        self.program = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/Enrichment/go.pl"

    def makeCommand(self,inputfq):
        deg=genediffexp()
        deg.fqLink=self.fqLink
        deg.species=self.species
        deg.outdir = self.outdir.replace("GO_Hypergeometric/GO", "GeneDiffExp_Allin")
        GeneDiffExpFilter = inputfq[0]
        novel_g2t=inputfq[1]
        novel_go=inputfq[2].split(".")[0]

        self.parameter["Annotation_dbClass"]=self.species["RNAref"][1]

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)

        dbclass.get_config(species=self.species["RNAref"][0])


        cmd=[]
        output=[]
        GODict={}
        go_shell=""

        out_dir=self.outdir+"/GO"
        go_tmpdir=self.outdir+"/tmp_file"
        os.makedirs(go_tmpdir, mode=0o755, exist_ok=True)
        os.makedirs(out_dir,mode=0o755, exist_ok=True)
        go_shell+="export PATH=/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/perl-V5/bin:$PATH;" \
                      "perl {scriptbin}/merge_gene2tr.pl {novel_g2t},{gene2tr} {tmpdir}/gene2tr;" \
                      "perl {scriptbin}/merge_go.pl {tmpdir}/gene2tr {novel_go},{go_prefix} {tmpdir}/species;".format(
            scriptbin=self.scriptbin,
            novel_g2t=novel_g2t,
            gene2tr=dbclass.Gene2Tr,
            tmpdir=go_tmpdir,
            novel_go=novel_go,
            go_prefix=dbclass.GO_PREFIX
        )
        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            go_shell+="awk \'{{if($5>0) print $1\"\\t\"$5\"\\tup\";else print $1\"\\t\"$5\" down\"}}\'" \
                     " {diff_list} >{tmpdir}/{keyname}.glist; " \
                     "perl {scriptbin}/drawGO.pl -list {tmpdir}/{keyname}.glist -goclass {goclass} " \
                     "-goprefix {prefix} -outprefix {outdir}/{keyname}; ".format(
                    diff_list=diff_list,
                    keyname=diff_id,
                    tmpdir=go_tmpdir,
                    prefix=go_tmpdir+"/species",
                    scriptbin=self.scriptbin,
                    outdir=self.outdir,
                    goclass=dbclass.GOCLASS
            )
            GODict[diff_id]=[self.outdir+"/"+diff_id+"_C.txt",self.outdir+"/"+diff_id+"_F.txt",self.outdir+"/"+diff_id+"_P.txt"]
        go_shell+="perl {scriptbin}/go.pl -gldir {tmpdir} -sdir `dirname {prefix}` -species `basename {prefix}` -outdir {outdir};" \
                  "perl {scriptbin}/topGO.pl -gldir {tmpdir} -godir {outdir} -outdir {outdir} -prefix {prefix}  -list {gene2tr} -outdir {outdir}".format(
            tmpdir=go_tmpdir,
            prefix=go_tmpdir+"/species",
            scriptbin=self.scriptbin,
            outdir=self.outdir,
            goclass=dbclass.GOCLASS,
            gene2tr=go_tmpdir+"/gene2tr"
        )

        cmd.append(go_shell)
        output.append(GODict)
        return cmd,output

    def makedefault(self,inputfq):
        input=[]
        output=[]
        deg=genediffexp()
        deg.species=self.species
        deg.fqLink=self.fqLink
        deg.outdir=self.outdir.replace("GO_Hypergeometric/GO","GeneDiffExp_Allin")
        novel_o=novel_tr()
        novel_o.species=self.species
        novel_o.fqLink=self.fqLink
        novel_o.outdir=self.outdir.replace("GO_Hypergeometric/GO","NovelTr_Cuffcompare")

        GeneDiffExpFilter = deg.makedefault(inputfq)["output"]
        novel_g2t = novel_o.makedefault(inputfq)["output"][-3]
        novel_go_file = novel_o.makedefault(inputfq)["output"][0]

        self.parameter["Annotation_dbClass"] = self.species["RNAref"][1]

        input.append(GeneDiffExpFilter)
        input.append(novel_g2t)
        input.append(novel_go_file)

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)

        dbclass.get_config(species=self.species["RNAref"][0])

        GODict={}
        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            GODict[diff_id]=[self.outdir+"/"+diff_id+"_C.txt",self.outdir+"/"+diff_id+"_F.txt",self.outdir+"/"+diff_id+"_P.txt"]
        output.append(GODict)

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,1CPU",
            'output': output
        }
        return default

class pathwayenrichment(common):
    def __init__(self):
        super(pathwayenrichment,self).__init__()
        self.outdir = "RNAref/Pathway_Hypergeometric/Pathway"

        self.parameter = {
            "Annotation_dbClass":"pl"
        }

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/Enrichment"
        self.program = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/Enrichment/pathfind.pl"

    def makeCommand(self,inputfq):
        GeneDiffExpFilter = inputfq[0]
        novel_g2t = inputfq[1]
        novel_ko = inputfq[2]

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)

        dbclass.get_config(species=self.species["RNAref"][0])

        cmd=[]
        output=[]
        KEGGDict={}
        kegg_shell=""
        kegg_pre_shell=""
        keggdraw_shell=""
        kegg_tmpdir = self.outdir + "/tmp_file"
        os.makedirs(kegg_tmpdir, mode=0o755, exist_ok=True)
        kegg_pre_shell+="perl {scriptbin}/merge_gene2tr.pl {novel_g2t},{gene2tr} {tmpdir}/gene2tr;" \
                    "perl {scriptbin}/merge_ko.pl {tmpdir}/gene2tr {novel_ko},{species_ko} {tmpdir}/bg.ko;".format(
            scriptbin=self.scriptbin,
            novel_g2t=novel_g2t,
            gene2tr=dbclass.Gene2Tr,
            novel_ko=novel_ko,
            species_ko=dbclass.KO,
            tmpdir=kegg_tmpdir
        )
        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            os.makedirs(self.outdir + '/'+diff_id+'_map',mode=0o755, exist_ok=True)
            kegg_shell+="awk \'{{print $1\"\\t\"$5}}\' {diff_list} >{tmpdir}/{keyname}.glist; " \
                       "perl {scriptbin}/getKO.pl -glist {tmpdir}/{keyname}.glist -bg {bg_ko} -outdir {tmpdir}; " \
                       "perl {scriptbin}/pathfind.pl -kegg {kegg_fa} -komap {ko_map} -maptitle {map_title} -fg {tmpdir}/{keyname}.ko " \
                       "-bg {bg_ko} -output {outdir}/{keyname}.path; " \
                       "awk '{{if($5>0)print $1\"\\t\"$5\"\\tup\";else print $1\"\\t\"$5\"\\tdown\"}} \' {diff_list}" \
                       " > {tmpdir}/{keyname}.glist.temp; " \
                       "perl {scriptbin}/keggMap.pl -ko {tmpdir}/{keyname}.ko -diff {tmpdir}/{keyname}.glist -komap " \
                       "{ko_map} -mapdir {mapdir} -outdir {outdir}/{keyname}_map; " \
                       "perl {scriptbin}/drawKEGG.pl -path {outdir}/{keyname}.path -outprefix {outdir}/{keyname}  " \
                       "-idCol 6 -level1Col 7 -level2Col 8 -geneCol 9 -list {tmpdir}/{keyname}.glist.temp\n" \
                       "".format(
                diff_list=diff_list,
                tmpdir=kegg_tmpdir,
                keyname=diff_id,
                bg_ko=kegg_tmpdir+"/bg.ko",
                scriptbin=self.scriptbin,
                outdir=self.outdir,
                kegg_fa=dbclass.KEGG,
                ko_map=dbclass.KEGG_KOMAP,
                map_title=dbclass.KEGG_MAP_TITLE,
                mapdir=dbclass.KEGG_MAP_DIR
            )
            KEGGDict[diff_id] = [self.outdir + "/" + diff_id + ".path"]

        keggdraw_shell +="perl {scriptbin}/genPathHTML.pl -indir {outdir};" \
                     "perl {scriptbin}/pathway_enrichFigure.pl {outdir}".format(
            scriptbin=self.scriptbin,
            outdir=self.outdir
        )

        cmd.append(kegg_pre_shell)
        cmd.append(kegg_shell)
        cmd.append(keggdraw_shell)
        output.append(KEGGDict)
        return cmd,output

    def makedefault(self,inputfq):

        deg=genediffexp()
        deg.species=self.species
        deg.fqLink=self.fqLink
        deg.outdir=self.outdir.replace("Pathway_Hypergeometric/Pathway","GeneDiffExp_Allin")

        novel_o=novel_tr()
        novel_o.species=self.species
        novel_o.fqLink=self.fqLink
        novel_o.outdir=self.outdir.replace("Pathway_Hypergeometric/Pathway","NovelTr_Cuffcompare")

        GeneDiffExpFilter = deg.makedefault(inputfq)["output"]
        novel_g2t = novel_o.makedefault(inputfq)["output"][-3]
        novel_ko = novel_o.makedefault(inputfq)["output"][-4]

        input=[]
        input.append(GeneDiffExpFilter)
        input.append(novel_g2t)
        input.append(novel_ko)

        self.parameter["Annotation_dbClass"]=self.species["RNAref"][1]

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)

        dbclass.get_config(species=self.species["RNAref"][0])

        output=[]
        KEGGDict={}
        out_dir=self.outdir

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            KEGGDict[diff_id] = [out_dir + "/" + diff_id + ".path"]
        output.append(KEGGDict)


        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,1CPU",
            'output': output
        }
        return default

class ppi(common):
    def __init__(self):
        super(ppi,self).__init__()
        self.outdir = "RNAref/PPI_Interaction"

        self.parameter=" --evalue 1e-5 --outfmt 6 --threads 5  --more-sensitive -b 0.5 --salltitles "

        self.diamond=self.getsoftware().DIAMOND
        self.python=self.getsoftware().PYTHON
        self.convert=self.getsoftware().CONVERT

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/PPI"

        self.program = [self.diamond,self.python,self.convert]

    def makeCommand(self,inputfq):

        database = DataBasePath()
        database.get_config(species=self.species["RNAref"][0])

        cmd = []
        output = []
        cds_fa=inputfq[0]
        gene2tr=inputfq[1]
        novel_fa = inputfq[2]
        novel_g2t = inputfq[3]
        GeneDiffExpFilter = inputfq[4]
        stringdb = inputfq[5]

        net_prepare_shell=""
        networkshell=""
        os.makedirs(self.outdir,mode=0o755,exist_ok=True)
        net_prepare_shell +="cat {cds_fa} {novel_fa} >{outdir}/all.fa;" \
                            "cat {gene2tr} {novel_g2t} >{outdir}/gene2tr;" \
                            "{diamond} blastx {parameter} -d {stringdb} -q {outdir}/all.fa -o {outdir}/ppi.outfmt6 --seg no --max-target-seqs 1;" \
                            "perl {scriptbin}/dealfmt6.pl {outdir}/gene2tr {outdir}/ppi.outfmt6 {outdir}/ppi.outfmt6.deal;".format(
            cds_fa=cds_fa,
            diamond=self.diamond,
            novel_fa=novel_fa,
            outdir=self.outdir,
            gene2tr=gene2tr,
            novel_g2t=novel_g2t,
            parameter=self.parameter,
            stringdb=stringdb,
            scriptbin=self.scriptbin
        )
        cmd.append(net_prepare_shell)

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            networkshell +="{python} {scriptbin}/RNA_DENOVO_PPI.py --m6 {outdir}/ppi.outfmt6.deal --db {stringdbre} --protein_list " \
                           "{diff_list} --diff_list {diff_list} --result {outdir}/{diffid}.network.relation.xls;" \
                           "cat {outdir}/{diffid}.network.relation.xls|sed -r 's/\\t/\\ttt\\t/'|awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$6}}' >{outdir}/{diffid}.Cytoscape.xls;" \
                           "awk 'NR>1{{print $0}}' {outdir}/{diffid}.network.relation.xls|sort -nrk 5|head -100 >{outdir}/{diffid}.network.relation.xls.100;" \
                           "{python} {scriptbin}/fd_relaion.py {outdir}/{diffid}.network.relation.xls.100 {diff_list} >{outdir}/{diffid}.network.relation.xls.100.fd;" \
                           "{scriptbin}/Rscript {scriptbin}/PPI.R {outdir}/{diffid}.network.relation.xls.100 {outdir}/{diffid}.network.relation.xls.100.fd {outdir}/{diffid}.network.pdf;" \
                           "{convert} -density 300 -resize 30% {outdir}/{diffid}.network.pdf {outdir}/{diffid}.network.png\n".format(
                python=self.python,
                scriptbin=self.scriptbin,
                outdir=self.outdir,
                stringdbre=database.STRING_RELATION,
                diff_list=diff_list,
                diffid=diff_id,
                convert=self.convert
            )
            output.append(self.outdir+"/"+diff_id+".network.relation.xls")
            output.append(self.outdir + "/" + diff_id + ".Cytoscape.xls")
            output.append(self.outdir + "/" + diff_id + ".network.pdf")
            output.append(self.outdir + "/" + diff_id + ".network.png")
        cmd.append(networkshell)
        return  cmd,output

    def makedefault(self,inputfq):

        database = DataBasePath()
        stringdb = database.STRING
        database.get_config(species=self.species["RNAref"][0])

        input=[]
        output = []
        cds_fa=database.CDS
        gene2tr=database.Gene2Tr
        novel_o = novel_tr()
        novel_o.species=self.species
        novel_o.fqLink=self.fqLink
        novel_o.outdir=self.outdir.replace("PPI_Interaction","NovelTr_Cuffcompare")
        novel_fa = novel_o.makedefault(inputfq)["output"][-1]
        novel_g2t = novel_o.makedefault(inputfq)["output"][-3]
        genediffexp_o = genediffexp()
        genediffexp_o.species=self.species
        genediffexp_o.fqLink=self.fqLink
        genediffexp_o.outdir=self.outdir.replace("PPI_Interaction","GeneDiffExp_Allin")
        GeneDiffExpFilter = genediffexp_o.makedefault(inputfq)["output"]

        input.append(cds_fa)
        input.append(gene2tr)
        input.append(novel_fa)
        input.append(novel_g2t)
        input.append(GeneDiffExpFilter)
        input.append(stringdb)

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            output.append(self.outdir+"/"+diff_id+".network.relation.xls")
            output.append(self.outdir + "/" + diff_id + ".Cytoscape.xls")
            output.append(self.outdir + "/" + diff_id + ".network.pdf")
            output.append(self.outdir + "/" + diff_id + ".network.png")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "4G,5CPU",
            'output': output
        }
        return default

class genefusion(common):
    def __init__(self):
        super(genefusion,self).__init__()
        self.outdir = "RNAref/GeneFusion_SOAPfuse"

        self.parameter=""

        self.soapfuse=self.getsoftware().SOAPFUSE
        self.circos=self.getsoftware().CIRCOS

        self.program=[self.soapfuse,self.circos]

    def makeCommand(self, inputfq):
        os.makedirs(self.outdir,mode=0o755,exist_ok=True)
        soapfuse_shell=""
        cmd=[]
        output=[]
        spe=""
        for spe in self.species.keys():
            pass
        if spe == "human" or spe== "Homo_sapiens" or spe == "mouse" or spe == "rat" or spe == "OryzaSativa" or spe == "arabidopsis":
            genefusedic={}
            for SampleID, Cleandata in inputfq[0]:
                os.makedirs(self.outdir + SampleID + "/Library_" + SampleID,mode=0o755,exist_ok=True)
                soapfuse_shell += "ln -s {cleandata_fq1} {sample_lib_dir}/{sampleid}_1.fq.gz;" \
                                  "ln -s {cleandata_fq2} {sample_lib_dir}/{sampleid}_2.fq.gz;" \
                                  "echo -e \"{sampleid}\\tLibrary_{sampleid}\\t{sampleid}\\t90\" >{outdir}/{sampleid}/InputFq_{sampleid}.list;" \
                                  "perl {soapfuse}/SOAPfuse-RUN.pl -fd {outdir} -l {outdir}/{sampleid}/InputFq_{sampleid}.list -c {soapfuse}/config/config.txt -fs 1 -tp 9 -fm -o {outdir}/{sampleid};" \
                                  "perl {Circos}/Perl2CircosGeneFusion.pl --input {outdir}/final_fusion_genes/{sampleid}/{sampleid}.final.Fusion.specific.for.genes --species {species} --output {outdir}/{sampleid};" \
                                  "".format(
                    cleandata_fq1=Cleandata["clean_fq1"],
                    cleandata_fq2=Cleandata["clean_fq2"],
                    sample_lib_dir=self.outdir + SampleID + "/Library_" + SampleID,
                    sampleid=SampleID,
                    soapfuse=self.soapfuse,
                    Circos=self.circos,
                    outdir=self.outdir,
                    species=spe
                )
                genefusedic[SampleID] = []
                genefusedic[SampleID].append(
                    self.outdir + "/" + SampleID + "final_fusion_genes/" + SampleID + "/" + SampleID + ".final.Fusion.specific.for.genes")
                genefusedic[SampleID].append(
                    self.outdir + "/" + SampleID + "final_fusion_genes/" + SampleID + "/" + SampleID + ".final.Fusion.specific.trans")
                genefusedic[SampleID].append(
                    self.outdir + "/" + SampleID + "final_fusion_genes/" + SampleID + "/" + SampleID + ".GeneFusion.png")
            output.append(genefusedic)
        else:
            soapfuse_shell += "echo \"This is just for model organism such as human,Homo_sapiens,mouse,arabidopsis analysis\" >{0}/genefusion.log".format(self.outdir)
            output.append(self.outdir+"/genefusion.log")
        cmd.append(soapfuse_shell)
        return  cmd,output

    def makedefault(self, inputfq):
        filter_o = filter()
        filter_o.species=self.species
        filter_o.outdir=self.outdir.replace("GeneFusion_SOAPfuse","Filter_SOAPnuke")
        filter_o.fqLink=self.fqLink
        CleanDataDict = filter_o.makedefault(inputfq)["output"][0]

        input=[]
        input.append(CleanDataDict)
        output=[]
        spe=""
        for spe in self.species.keys():
            pass
        if spe == "human" or spe == "Homo_sapiens" or spe == "mouse" or spe == "rat" or spe == "OryzaSativa" or spe == "arabidopsis":
            genefusedic={}
            for SampleID, Cleandata in CleanDataDict.items():
                genefusedic[SampleID] =[]
                genefusedic[SampleID].append(self.outdir+"/"+SampleID+"final_fusion_genes/"+SampleID+"/"+SampleID+".final.Fusion.specific.for.genes")
                genefusedic[SampleID].append(self.outdir + "/" + SampleID + "final_fusion_genes/" + SampleID + "/" + SampleID + ".final.Fusion.specific.trans")
                genefusedic[SampleID].append(self.outdir + "/" + SampleID + "final_fusion_genes/" + SampleID+"/"+SampleID+".GeneFusion.png")
            output.append(genefusedic)
        else:
            output.append(self.outdir + "/genefusion.log")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,1CPU",
            'output': output
        }
        return default

class prg(common):
    def __init__(self):
        super(prg, self).__init__()
        self.outdir = "RNAref/PRG_Blastx"
        self.diamond=self.getsoftware().DIAMOND
        self.program=[self.diamond]
        self.ref=""
        self.cds=""
        self.gene2tr=""
        self.parameter = {
            "Annotation_dbClass": "pl"
        }

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/PRG/bin"

    def makeCommand(self, inputfq):

        novel_tr_fa = inputfq[0]
        novel_gene2tr = inputfq[1]

        deglis = inputfq[2]

        os.makedirs(self.outdir, mode=0o775, exist_ok=True)

        self.parameter["Annotation_dbClass"] = self.species["RNAref"][1]

        spe=self.species
        database = DataBasePath(dbclass=self.parameter["Annotation_dbClass"])
        database.get_config()
        database.get_config(species=spe["RNAref"][0])

        self.ref=database.GENOME
        self.cds=database.CDS
        self.gene2tr=database.Gene2Tr

        prg_sh = ""
        cmd = []
        output = []

        if self.parameter["Annotation_dbClass"] == "an":
            prg_sh +="echo \"Your sample is not plant. You don't need to do this step\" >{0}/prg.log;".format(self.outdir)
            output.append(self.outdir+"/prg.log")
        elif self.parameter["Annotation_dbClass"] == "pl":
            os.makedirs(self.outdir+"/tmp", mode=0o775, exist_ok=True)
            for deg in deglis:
                prg_sh+="cp {deg} {outdir}/tmp/;".format(deg=deg,outdir=self.outdir)
            prg_sh += "cat {cds} {novel_fa} >{outdir}/reference.fa;" \
                      "cat {cds_gene2tr} {novel_gene2tr} >{outdir}/reference.gene2tr;" \
                      "perl {scriptbin}/get_all_deg.pl {outdir}/reference.fa {deg_dir} {outdir}/degtrseq.fa {outdir}/reference.gene2tr;" \
                      "{blastx} blastx --query {outdir}/degtrseq.fa --db {prgdb} --evalue 1e-5 --threads 6 --outfmt 6 --max-target-seqs 1 --more-sensitive --out {outdir}/degtr2PRG.out;" \
                      "perl {scriptbin}/get_PRG.pl -trseq {outdir}/degtrseq.fa -blresult {outdir}/degtr2PRG.out -gene2tr {outdir}/reference.gene2tr -degdir {deg_dir} -anno {anno} -cov 50 -iden 40 -output {outdir}"\
                .format(
                blastx=self.diamond,
                cds=self.cds,
                novel_fa=novel_tr_fa,
                anno=database.PRG_AN,
                cds_gene2tr=self.gene2tr,
                novel_gene2tr=novel_gene2tr,
                deg_dir=self.outdir+"/tmp",
                outdir=self.outdir,
                prgdb=database.PRG,
                scriptbin=self.scriptbin
            )
            output.append(self.outdir + "/reference.fa")
        elif self.parameter["Annotation_dbClass"] == "fg":
            prg_sh += "echo \"You don't need to do this step\" >{0}/prg.log".format(self.outdir)
            output.append(self.outdir + "/prg.log")
        else:
            sys.exit(1)
        cmd.append(prg_sh)

        return cmd,output

    def makedefault(self, inputfq):
        novel_o = novel_tr()
        novel_o.species=self.species
        novel_o.fqLink=self.fqLink
        novel_o.outdir=self.outdir.replace("PRG_Blastx","NovelTr_Cuffcompare")
        novel_tr_fa = novel_o.makedefault(inputfq)["output"][-1]
        novel_gene2tr = novel_o.makedefault(inputfq)["output"][-3]

        genediffexp_o = genediffexp()
        genediffexp_o.species=self.species
        genediffexp_o.fqLink=self.fqLink
        genediffexp_o.outdir=self.outdir.replace("PRG_Blastx","GeneDiffExp_Allin")
        deg_lis = genediffexp_o.makedefault(inputfq)["output"]
        input=[]
        output=[]

        self.parameter["Annotation_dbClass"]=self.species["RNAref"][1]

        input.append(novel_tr_fa)
        input.append(novel_gene2tr)
        input.append(deg_lis)
        if self.parameter["Annotation_dbClass"] == "an":
            output.append(self.outdir + "/prg.log")
        elif self.parameter["Annotation_dbClass"] == "pl":
            output.append(self.outdir + "/reference.fa")
        elif self.parameter["Annotation_dbClass"] == "fg":
            output.append(self.outdir + "/prg.log")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "1G,4CPU",
            'output': output
        }

        return  default

class phi(common):
    def __init__(self):
        super(phi, self).__init__()
        self.outdir = "RNAref/PHI_Blast"
        self.diamond=self.getsoftware().DIAMOND
        self.program=[self.diamond]
        self.ref=""
        self.cds=""
        self.gene2tr=""
        self.parameter = {
            "Annotation_dbClass": "fg"
        }

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/PHI/bin"

    def makeCommand(self, inputfq):

        novel_tr_fa = inputfq[0]
        novel_gene2tr = inputfq[1]

        deglis = inputfq[2]

        os.makedirs(self.outdir, mode=0o775, exist_ok=True)

        self.parameter["Annotation_dbClass"] = self.species["RNAref"][1]

        spe=self.species
        database = DataBasePath(dbclass=self.parameter["Annotation_dbClass"])
        database.get_config()
        database.get_config(species=spe["RNAref"][0])

        self.ref=database.GENOME
        self.cds=database.CDS
        self.gene2tr=database.Gene2Tr

        phi_sh = ""
        cmd = []
        output = []

        if self.parameter["Annotation_dbClass"] == "an":
            phi_sh +="echo \"You don't need to do this step, this is for fg.\" >{0}/phi.log".format(self.outdir)
            output.append(self.outdir+"/phi.log")
        elif self.parameter["Annotation_dbClass"] == "fg":
            os.makedirs(self.outdir+"/tmp", mode=0o775, exist_ok=True)
            for deg in deglis:
                phi_sh+="cp {deg} {outdir}/tmp/".format(deg=deg,outdir=self.outdir)
            phi_sh += "cat {cds} {novel_fa} >{outdir}/reference.fa;" \
                      "cat {cds_gene2tr} {novel_gene2tr} >{outdir}/reference.gene2tr;" \
                      "perl {scriptbin}/get_all_deg.pl {outdir}/reference.fa {deg_dir} {outdir}/degtrseq.fa {outdir}/reference.gene2tr;" \
                      "{blastx} blastx -query {outdir}/degtrseq.fa -db {phidb} --evalue 1e-5 --threads 4 --outfmt 6 --max-target-seqs 1 --more-sensitive -out {outdir}/degtr2PHI.out;" \
                      "perl {scriptbin}/get_PHI.pl -trseq {outdir}/degtrseq.fa -blresult {outdir}/degtr2PHI.out -gene2tr {outdir}/reference.gene2tr -degdir {deg_dir} -cov 50 -iden 40 -output {outdir}"\
                .format(
                blastx=self.diamond,
                cds=self.cds,
                novel_fa=novel_tr_fa,
                cds_gene2tr=self.gene2tr,
                novel_gene2tr=novel_gene2tr,
                deg_dir=self.outdir+"/tmp",
                outdir=self.outdir,
                phidb=database.PHI,
                scriptbin=self.scriptbin,
            )
            output.append(self.outdir + "/reference.fa")
        elif self.parameter["Annotation_dbClass"] == "pl":
            phi_sh +="echo \"You don't need to do this step, this is for fg\" >{0}/phi.log".format(self.outdir)
            output.append(self.outdir+"/phi.log")
        else:
            sys.exit(1)
        cmd.append(phi_sh)
        return cmd,output

    def makedefault(self, inputfq):
        novel_o = novel_tr()
        novel_o.species=self.species
        novel_o.fqLink=self.fqLink
        novel_o.outdir=self.outdir.replace("PHI_Blast","NovelTr_Cuffcompare")
        novel_tr_fa = novel_o.makedefault(inputfq)["output"][-1]
        novel_gene2tr = novel_o.makedefault(inputfq)["output"][-3]

        genediffexp_o = genediffexp()
        genediffexp_o.species=self.species
        genediffexp_o.fqLink=self.fqLink
        genediffexp_o.outdir=self.outdir.replace("PHI_Blast","GeneDiffExp_Allin")
        deg_lis = genediffexp_o.makedefault(inputfq)["output"]
        input=[]
        output=[]

        self.parameter["Annotation_dbClass"]=self.species["RNAref"][1]

        input.append(novel_tr_fa)
        input.append(novel_gene2tr)
        input.append(deg_lis)
        if self.parameter["Annotation_dbClass"] == "an":
            output.append(self.outdir + "/phi.log")
        elif self.parameter["Annotation_dbClass"] == "fg":
            output.append(self.outdir+"/reference.fa")
        elif self.parameter["Annotation_dbClass"] == "pl":
            output.append(self.outdir + "/phi.log")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "1G,4CPU",
            'output': output
        }

        return  default

class circos(common):
    def __init__(self):
        super(circos,self).__init__()
        self.outdir = "RNAref/CircosFig"

        self.parameter = {
            "Figure_Type": 2,
            "Window_Size" : 1000000,
            "Chromosome_Unit" : 1000000,
            "Chr_Label_Size" : "15p",
            "Chr_Tick_Size" : "8p"
        }

        self.ref=""
        self.gtf=""

        self.soapfuse=self.getsoftware().SOAPFUSE
        self.circos=self.getsoftware().CIRCOS

        self.program=[self.soapfuse,self.circos]
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/Circos/bin"

    def makeCommand(self, inputfq):
        os.makedirs(self.outdir,mode=0o755,exist_ok=True)
        cmd=[]
        output=[]

        geneexp_o = inputfq[0]
        genefuse_o = inputfq[1]
        snpindel_o = inputfq[2]
        circos_shell = ""

        spe = self.species["RNAref"][0]
        type = self.species["RNAref"][1]

        database = DataBasePath(dbclass=type)
        database.get_config(species=spe)
        if spe == "human" and self.parameter["Figure_Type"] == 1:
            for SampleID in geneexp_o.keys():
                os.makedirs(self.outdir + "/" + SampleID, mode=0o755, exist_ok=True)
                circos_shell += "perl {scriptbin}/seq_basenum.pl {ref} {outdir}/{sampleid}/chr_basenum.txt;" \
                                "perl {scriptbin}/circos_karyotype.pl {outdir}/{sampleid}/chr_basenum.txt {outdir}/{sampleid}/karyotype.txt;" \
                                "grep -v '^#' {snp_vcf} | awk '{{print $1\"\\t\"$2\"\\t\"$4\"\\t\"$5}}' > {outdir}/{sampleid}/{sampleid}.snp.txt;" \
                                "grep -v '^#' {indel_vcf} | awk '{{print $1\"\\t\"$2\"\\t\"$4\"\\t\"$5}}' >{outdir}/{sampleid}/{sampleid}.indel.txt;" \
                                "perl {scriptbin}/snp_indel_num.pl {outdir}/{sampleid}/{sampleid}.snp.txt {outdir}/{sampleid}/chr_basenum.txt {outdir}/{sampleid} {window};" \
                                "cat {outdir}/{sampleid}/{sampleid}.snp.*.xls >{outdir}/{sampleid}/{sampleid}.snp.circos.txt;" \
                                "rm {outdir}/{sampleid}/{sampleid}.snp.*.xls;" \
                                "perl {scriptbin}/snp_indel_num.pl {outdir}/{sampleid}/{sampleid}.indel.txt {outdir}/{sampleid}/chr_basenum.txt {outdir}/{sampleid} {window};" \
                                "cat {outdir}/{sampleid}/{sampleid}.indel.*.xls >{outdir}/{sampleid}/{sampleid}.indel.circos.txt;" \
                                "rm {outdir}/{sampleid}/{sampleid}.indel.*.xls;" \
                                "perl {scriptbin}/gtf_geneloci.pl {gtf} {outdir}/{sampleid}/geneloci.txt;" \
                                "perl {scriptbin}/geneloci_exp.pl {fpkm} {outdir}/{sampleid}/geneloci.txt {outdir}/{sampleid}/{sampleid}.exp.txt;" \
                                "perl {scriptbin}/exp_value.pl {outdir}/{sampleid}/{sampleid}.exp.txt {outdir}/{sampleid}/chr_basenum.txt {outdir}/{sampleid} {window};" \
                                "cat {outdir}/{sampleid}/{sampleid}.exp.*.xls >{outdir}/{sampleid}/{sampleid}.exp.sum.txt;" \
                                "rm {outdir}/{sampleid}/{sampleid}.exp.*.xls;" \
                                "perl {scriptbin}/exp_logvalue.pl {outdir}/{sampleid}/{sampleid}.exp.sum.txt {outdir}/{sampleid}/{sampleid}.exp.circos.txt;" \
                                "grep -v 'up_gene' {genefuse} | awk '{{print $2\"\\t\"$4\"\\t\"$4\"\\t\"$1\"\\n\"$7\"\\t\"$9\"\\t\"$9\"\\t\"$6}}' >{outdir}/{sampleid}/{sampleid}.labels.circos.txt;" \
                                "grep -v 'up_gene' {genefuse} | awk '{{print $2\"\\t\"$4\"\\t\"$4\"\\t\"$7\"\\t\"$9\"\\t\"$9}}' > {outdir}/{sampleid}/{sampleid}.fusions.circos.txt;" \
                                "perl {scriptbin}/draw_circos.pl --kary {outdir}/{sampleid}/karyotype.txt --unit {unit} --chrlab {chrlab} --tick {tick} --snp {outdir}/{sampleid}/{sampleid}.snp.circos.txt --indel {outdir}/{sampleid}/{sampleid}.indel.circos.txt --exp {outdir}/{sampleid}/{sampleid}.exp.txt --label {outdir}/{sampleid}/{sampleid}.labels.circos.txt --fusion {outdir}/{sampleid}/{sampleid}.fusions.circos.txt --output {outdir}/{sampleid}\n".format(
                    scriptbin=self.scriptbin,
                    ref=database.GENOME,
                    outdir=self.outdir,
                    sampleid=SampleID,
                    snp_vcf=snpindel_o[SampleID][0],
                    indel_vcf=snpindel_o[SampleID][1],
                    gtf=database.mRNA_GTF,
                    fpkm=geneexp_o[SampleID],
                    genefuse=genefuse_o[SampleID][0],
                    unit=self.parameter["Chromosome_Unit"],
                    window=self.parameter["Window_Size"],
                    chrlab=self.parameter["Chr_Label_Size"],
                    tick=self.parameter["Chr_Tick_Size"]
                )
                output.append(self.outdir + "/" + SampleID + "/" + SampleID + ".circos.svg")
                output.append(self.outdir + "/" + SampleID + "/" + SampleID + ".circos.png")
        elif spe == "human" and self.parameter["Figure_Type"] == 2:
            for SampleID in geneexp_o.keys():
                os.makedirs(self.outdir + "/" + SampleID, mode=0o755, exist_ok=True)
                circos_shell += "perl {scriptbin}/seq_basenum.pl {ref} {outdir}/{sampleid}/chr_basenum.txt;" \
                                "perl {scriptbin}/circos_karyotype.pl {outdir}/{sampleid}/chr_basenum.txt {outdir}/{sampleid}/karyotype.txt;" \
                                "grep -v '^#' {snp_vcf} | awk '{{print $1\"\\t\"$2\"\\t\"$4\"\\t\"$5}}' > {outdir}/{sampleid}/{sampleid}.snp.txt;" \
                                "grep -v '^#' {indel_vcf} | awk '{{print $1\"\\t\"$2\"\\t\"$4\"\\t\"$5}}' >{outdir}/{sampleid}/{sampleid}.indel.txt;" \
                                "perl {scriptbin}/snp_indel_num.pl {outdir}/{sampleid}/{sampleid}.snp.txt {outdir}/{sampleid}/chr_basenum.txt {outdir}/{sampleid} {window};" \
                                "cat {outdir}/{sampleid}/{sampleid}.snp.*.xls >{outdir}/{sampleid}/{sampleid}.snp.circos.txt;" \
                                "rm {outdir}/{sampleid}/{sampleid}.snp.*.xls;" \
                                "perl {scriptbin}/snp_indel_num.pl {outdir}/{sampleid}/{sampleid}.indel.txt {outdir}/{sampleid}/chr_basenum.txt {outdir}/{sampleid} {window};" \
                                "cat {outdir}/{sampleid}/{sampleid}.indel.*.xls >{outdir}/{sampleid}/{sampleid}.indel.circos.txt;" \
                                "rm {outdir}/{sampleid}/{sampleid}.indel.*.xls;" \
                                "perl {scriptbin}/gtf_geneloci.pl {gtf} {outdir}/{sampleid}/geneloci.txt;" \
                                "perl {scriptbin}/geneloci_exp.pl {fpkm} {outdir}/{sampleid}/geneloci.txt {outdir}/{sampleid}/{sampleid}.exp.txt;" \
                                "perl {scriptbin}/exp_value.pl {outdir}/{sampleid}/{sampleid}.exp.txt {outdir}/{sampleid}/chr_basenum.txt {outdir}/{sampleid} {window};" \
                                "cat {outdir}/{sampleid}/{sampleid}.exp.*.xls >{outdir}/{sampleid}/{sampleid}.exp.sum.txt;" \
                                "rm {outdir}/{sampleid}/{sampleid}.exp.*.xls;" \
                                "perl {scriptbin}/exp_logvalue.pl {outdir}/{sampleid}/{sampleid}.exp.sum.txt {outdir}/{sampleid}/{sampleid}.exp.circos.txt;" \
                                "perl {scriptbin}/draw_circos.pl --kary {outdir}/{sampleid}/karyotype.txt --unit {unit} --chrlab {chrlab} --tick {tick} --snp {outdir}/{sampleid}/{sampleid}.snp.circos.txt --indel {outdir}/{sampleid}/{sampleid}.indel.circos.txt --exp {outdir}/{sampleid}/{sampleid}.exp.txt --label {outdir}/{sampleid}/{sampleid}.labels.circos.txt --output {outdir}/{sampleid}\n".format(
                    scriptbin=self.scriptbin,
                    ref=database.GENOME,
                    outdir=self.outdir,
                    sampleid=SampleID,
                    snp_vcf=snpindel_o[SampleID][0],
                    indel_vcf=snpindel_o[SampleID][1],
                    gtf=database.mRNA_GTF,
                    fpkm=geneexp_o[SampleID],
                    unit=self.parameter["Chromosome_Unit"],
                    window=self.parameter["Window_Size"],
                    chrlab=self.parameter["Chr_Label_Size"],
                    tick=self.parameter["Chr_Tick_Size"]
                )
                output.append(self.outdir + "/" + SampleID + "/" + SampleID + ".circos.svg")
                output.append(self.outdir + "/" + SampleID + "/" + SampleID + ".circos.png")
        else:
            circos_shell+="echo \"This analysis only for human, your species need to set as human and Figure_Type only set as 1 or 2.\" >{0}/circos.log".format(self.outdir)
            output.append(self.outdir+"/circos.log")
        cmd.append(circos_shell)
        return  cmd,output

    def makedefault(self, inputfq):
        genefusion_o = genefusion()
        genefusion_o.outdir=self.outdir.replace("CircosFig","GeneFusion_SOAPfuse")
        genefusion_o.fqLink=self.fqLink
        genefusion_o.species=self.species
        genefusion_output=genefusion_o.makedefault(inputfq)["output"][0]

        snpindel_o = snpindel()
        snpindel_o.outdir=self.outdir.replace("CircosFig","SnpIndel_GATK")
        snpindel_o.species=self.species
        snpindel_o.fqLink=self.fqLink
        snpindel_output=snpindel_o.makedefault(inputfq)["output"][0]

        geneexp_o = geneexp()
        geneexp_o.outdir=self.outdir.replace("CircosFig","GeneExp")
        geneexp_o.species=self.species
        geneexp_o.fqLink=self.fqLink
        geneexp_output=geneexp_o.makedefault(inputfq)["output"][1]

        input=[]
        output=[]
        spe = self.species["RNAref"][0]
        if spe == "human" and self.parameter["Figure_Type"] == 1:
            for sampleid in geneexp_output.keys():
                output.append(self.outdir + "/" + sampleid + "/" + sampleid + ".circos.svg")
                output.append(self.outdir + "/" + sampleid + "/" + sampleid + ".circos.png")
        elif spe == "human" and self.parameter["Figure_Type"] == 2:
            for sampleid in geneexp_output.keys():
                output.append(self.outdir + "/" + sampleid + "/" + sampleid + ".circos.svg")
                output.append(self.outdir + "/" + sampleid + "/" + sampleid + ".circos.png")
        else:
            output.append(self.outdir + "/circos.log")

        input.append(geneexp_output)
        input.append(genefusion_output)
        input.append(snpindel_output)

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "1G,4CPU",
            'output': output
        }

        return  default

class cluster(common):
    def __init__(self):
        super(cluster, self).__init__()
        self.outdir = "RNAref/Clustering_Mfuzz"

        self.program="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/TimeCluster/TimeClustering.pl"
        self.parameter = {
            "Mfuzz_Options":"-c 12 -m 1.25",
            "Mfuzz_plan":""
        }

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/TimeCluster"

    def makeCommand(self, inputfq):
        geneexpdic=inputfq[0]
        cmd=[]
        output=[]
        cluster_shell=""
        plan_lis = self.parameter["Mfuzz_plan"].split(" ")
        pattern=re.compile(r'^AAAA')

        if len(plan_lis) >1:
            for i in range(len(plan_lis)):
                get_a=re.match(pattern,plan_lis[i])
                if get_a:
                    os.makedirs(self.outdir + "/AAAA", mode=0o775, exist_ok=True)
                    plan_id = plan_lis[i].split("(")[0].split(":")[1]
                    para = plan_lis[i].split("(")[1].replace(")", "").replace(":", " ").replace(",", " ")
                    for sampleid in geneexpdic.keys():
                        cluster_shell += "cp {fpkm} {outdir}/AAAA;".format(fpkm=geneexpdic[sampleid],
                                                                           outdir=self.outdir)
                    cluster_shell += "perl {program} -indir {outdir}/AAAA  -suffix .gene.fpkm.xls -plan {planid} {para} -allname AAAA -outdir {outdir}/AAAA".format(
                        program=self.program,
                        outdir=self.outdir,
                        planid=plan_id,
                        para=para
                    )
                    output.append(self.outdir + "/AAAA.cluster.xls")
                    output.append(self.outdir + "/AAAA.membership.xls")
                else:
                    os.makedirs(self.outdir + "/" + plan_lis[i], mode=0o775, exist_ok=True)
                    tmp = plan_lis[i].split(",")
                    for j in tmp:
                        cluster_shell+="cp {fpkm} {outdir}/{plandir};".format(fpkm=geneexpdic[j],outdir=self.outdir,plandir=plan_lis[i])
                    cluster_shell+="perl {program} -indir {outdir}/{plandir} -suffix .gene.fpkm.xls -plan {planid} {para}  -allname {planid} -outdir {outdir}/{plandir}\n".format(
                        program=self.program,
                        outdir=self.outdir,
                        plandir=plan_lis[i],
                        planid=plan_lis[i],
                        para=self.parameter["Mfuzz_Options"]
                    )
                    output.append(self.outdir+"/"+plan_lis[i]+".cluster.xls")
                    output.append(self.outdir + "/" + plan_lis[i] + ".membership.xls")
        else:
            os.makedirs(self.outdir+"/AAAA", mode=0o775, exist_ok=True)
            plan_id = self.parameter["Mfuzz_plan"].split("(")[0].split(":")[1]
            para=self.parameter["Mfuzz_plan"].split("(")[1].replace(")","").replace(":"," ").replace(","," ")
            for sampleid in geneexpdic.keys():
                cluster_shell+="cp {fpkm} {outdir}/AAAA;".format(fpkm=geneexpdic[sampleid],outdir=self.outdir)
            cluster_shell+="perl {program} -indir {outdir}/AAAA  -suffix .gene.fpkm.xls -plan {planid} {para} -allname AAAA -outdir {outdir}/AAAA".format(
                program=self.program,
                outdir=self.outdir,
                planid=plan_id,
                para=para
            )
            output.append(self.outdir+"/AAAA/AAAA.cluster.xls")
            output.append(self.outdir+"/AAAA/AAAA.membership.xls")
        cmd.append(cluster_shell)
        return cmd,output
    def makedefault(self, inputfq):
        geneexp_o =geneexp()
        geneexp_o.outdir=self.outdir.replace("Clustering_Mfuzz","GeneExp")
        geneexp_o.fqLink=self.fqLink
        geneexp_o.species=self.species
        geneexpdic = geneexp_o.makedefault(inputfq)["output"][1]
        input=[]
        output=[]
        dic_gro={}
        for i_key, line in self.fqLink.items():
            groid = line[0]
            samid = line[1]
            if groid in dic_gro.keys():
                dic_gro[groid].append(samid)
            else:
                dic_gro[groid] = []
                dic_gro[groid].append(samid)
        for groupid,samplelis in dic_gro.items():
            if len(samplelis) <3:
                pass
            else:
                tmp = ",".join(samplelis)
                self.parameter["Mfuzz_plan"] +=tmp+" "
        tmplis=[]
        for sampleid in geneexpdic.keys():
            tmplis.append(sampleid)
        self.parameter["Mfuzz_plan"] +="AAAA:"+",".join(tmplis)+"(-c:16,-m:1.35)"

        input.append(geneexpdic)
        output.append(self.outdir+"/AAAA/AAAA.cluster.xls")
        output.append(self.outdir+"/AAAA/AAAA.membership.xls")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "1G,1CPU",
            'output': output
        }

        return  default

class preresult(common):
    def __init__(self):
        super(preresult, self).__init__()
        self.parameter = ""
        self.program=""
        self.outdir = "BGI_result"
    def makeCommand(self, inputfq):
        outd = self.outdir.replace("/BGI_result", "")
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/1.CleanData", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/2.MapStat", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/3.Structure", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/4.Quantify", mode=0o755, exist_ok=True)
        cpshell=""
        cmd=[]
        output=[]
        cpshell +="mkdir -p {outdir}/2.MapStat/GeneMapping;" \
                  "mkdir -p {outdir}/2.MapStat/GenomeMapping;" \
                  "mkdir -p {outdir}/3.Structure/Circos;" \
                  "mkdir -p {outdir}/3.Structure/DifferentiallySplicingGene;" \
                  "mkdir -p {outdir}/3.Structure/GeneFusion;" \
                  "mkdir -p {outdir}/3.Structure/NovelTranscript;" \
                  "mkdir -p {outdir}/3.Structure/SnpIndel;" \
                  "mkdir -p {outdir}/3.Structure/GeneFusion;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/DEGList;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/DEGVenn;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/GeneOntolotyEnrichment;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/HierarchicalCluster;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/KeggPathwayEnrichment;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/PHI;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/PRG;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/PPI;" \
                  "mkdir -p {outdir}/4.Quantify/DifferentiallyExpressedGene/TFprediction;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/CorrelationHeatmap;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/GeneExpression;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/Clustering_Mfuzz;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/HclusterTree;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/PCA;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/ReadsCoverage;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/ReadsRandom;" \
                  "mkdir -p {outdir}/4.Quantify/GeneExpression/VennDiagram;" \
                  "cp {filteroutdir}/*/*.filter.stat.xls {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.RawReadsClass.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.base.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.qual.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/FilterSummary.xls {outdir}/1.CleanData/;" \
                  "cp {alignmentoutdir}/GenomeMappingSummary.xls {outdir}/2.MapStat/GenomeMapping/;" \
                  "cp {alignmentoutdir}/*/*AddRG.Reorder.Sort.bam {alignmentoutdir}/*/*AddRG.Reorder.Sort.bam.bai {outdir}/2.MapStat/GenomeMapping/;" \
                  "cp {noveloutdir}/Prediction/{{NovelTranscriptSummary.xls,novel_coding_transcript.fa}} {outdir}/3.Structure/NovelTranscript/;" \
                  "cp {snpindeloutdir}/*/{{*.snp.vcf,*.indel.vcf,*.annot.xls,*.annot.stat.pdf,*.annot.stat.png }} {outdir}/3.Structure/SnpIndel/;" \
                  "cp {snpindeloutdir}/{{combine.snp.vcf,combine.indel.vcf,snp_population.xlsx}} {outdir}/3.Structure/SnpIndel/;" \
                  "cp {geneexpoutdir}/*/{{*.Bowtie2Gene.MapReadsStat.xls,*.gene.fpkm.xls,*transcript.fpkm.xls}} {outdir}/2.MapStat/GeneMapping/;" \
                  "cp {geneexpoutdir}/*/{{*.gene.fpkm.xls,*transcript.fpkm.xls}} {outdir}/4.Quantify/GeneExpression/GeneExpression/;" \
                  "cp {geneexpoutdir}/*/{{*ReadsCoverage.pdf *ReadsCoverage.png}} {outdir}/4.Quantify/GeneExpression/ReadsCoverage/;" \
                  "cp {geneexpoutdir}/*/{{*ReadsRandom.pdf *ReadsRandom.png}} {outdir}/4.Quantify/GeneExpression/ReadsRandom/;" \
                  "cp {geneexpoutdir}/{{*.pdf,*.png,AllSamples.GeneExpression.FPKM.xls,AllSamples.TranscriptExpression.FPKM.xls,GeneExpressionSummary.xls,*pdf,*png}} {outdir}/4.Quantify/GeneExpression/GeneExpression/;" \
                  "cp {geneexpoutdir}/GeneExpression/CorrelationHeatmap/* {outdir}/4.Quantify/GeneExpression/CorrelationHeatmap/;" \
                  "cp {geneexpoutdir}/GeneExpression/HclusterTree/* {outdir}/4.Quantify/GeneExpression/HclusterTree/;" \
                  "cp {genediffexpoutdir}/*/*GeneDiffExp*xls {genediffexpoutdir}/*/*.MA-plot.* {genediffexpoutdir}/*/*.Scatter-plot.* {genediffexpoutdir}/*/*.Volcano-plot.* {outdir}/4.Quantify/DifferentiallyExpressedGene/DEGList;" \
                  "cp {genediffspliceoutdir}/DifferentiallySplicingGene {outdir}/3.Structure/DifferentiallySplicingGene;" \
                  "cp {tfoutdir}/{{*pdf,*png,*xls}} {outdir}/4.Quantify/DifferentiallyExpressedGene/TFprediction;" \
                  "cp {goenrichmentoutdir}/* {outdir}/4.Quantify/DifferentiallyExpressedGene/GeneOntolotyEnrichment;" \
                  "cp {pathwayenrichmentoutdir}/{{*.xls,*.htm,*.pdf,*.png,*.path,*map}} {outdir}/4.Quantify/DifferentiallyExpressedGene/KeggPathwayEnrichment;" \
                  "cp {ppioutdir}/{{*.network.pdf,*.network.png,*.network.relation.txt}} {outdir}/4.Quantify/DifferentiallyExpressedGene/PPI;" \
                  "cp {genefusionoutdir}/*/{{*.final.*,*.GeneFusion.png}} {outdir}/3.Structure/GeneFusion;" \
                  "cp {phioutdir}/Unigene2PHI.result.xls {outdir}/4.Quantify/DifferentiallyExpressedGene/PHI;" \
                  "cp {prgoutdir}/Unigene2PRG.result.xls {outdir}/4.Quantify/DifferentiallyExpressedGene/PRG;" \
                  "cp {circosoutdir}/*/{{*.circos.png ,*.circos.svg}} {outdir}/3.Structure/Circos;" \
                  "cp {clusteroutdir}/*/*{{*.cluster.xls,*.membership.xls,*.mfuzz.plot.pdf,*.mfuzz.plot.png}} {outdir}/4.Quantify/GeneExpression/Clustering_Mfuzz;".format(
            filteroutdir=outd+"/Filter_SOAPnuke",
            alignmentoutdir=outd+"/GenomeMapping_HISAT",
            noveloutdir = outd+"/NovelTr_Cuffcompare",
            snpindeloutdir=outd+"/SnpIndel_GATK",
            geneexpoutdir=outd+"/GeneExp",
            genediffexpoutdir=outd+"/GeneDiffExp_Allin",
            genediffspliceoutdir=outd+"/GeneDiffSplice",
            goenrichmentoutdir=outd+"/GO_Hypergeometric/GO",
            pathwayenrichmentoutdir=outd+"/Pathway_Hypergeometric/KEGG",
            ppioutdir=outd+"/PPI_Interaction",
            tfoutdir=outd+"/TFpredict",
            phioutdir=outd+"/PHI_Blastx",
            prgoutdir=outd+"/PRG_Blastx",
            circosoutdir=outd+"/CircosFig",
            genefusionoutdir=outd+"/GeneFusion_SOAPfuse",
            clusteroutdir=outd+"/Clustering_Mfuzz",
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
        self.input = "%s/workflow.json" % (self.outdirMain)
        self.output = "%s/workflow.json" % (self.outdirMain)

    def runlocal(self):
        outdir = self.outdirMain
        try:
            for stepL in self.step:
                for step in stepL:
                    astep = eval(step)
                    astepo = astep()
                    astepo.fqLink = self.fqLink
                    astepo.species= self.species
                    astepo.fqList=self.fqList
                    astepo.outdir = outdir +"/"+astepo.outdir
                    default = astepo.makedefault(self.fqList)
                    tmpcmd,tmpout = astepo.makeCommand(default['input'])
                    print (step +" start")
                    for i in range(len(tmpcmd)):
                        runcmdlist = tmpcmd[i].split("\n")
                        processNum=len(runcmdlist)
                        realp = math.ceil(processNum/2)
                        with Pool(realp) as p:
                            for j in range(len(runcmdlist)):
                                p.apply_async(run_cmd,args=(runcmdlist[j],))
                            print("Waiting for all subprocess done...")
                            p.close()
                            p.join()
                            print ("All subprocess done")
                    print (step + " finish")
        except IOError as e:
            raise e

    def makeshell(self, outputfile=None):
        outdir = self.outdirMain
        outputjson="%s/workflow.json" % (outdir)
        os.makedirs(outdir+"/shell",exist_ok=True,mode=0o755)
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
                    astepo.species= self.species
                    astepo.fqList=self.fqList
                    astepo.outdir = outdir +"/"+astepo.outdir
                    default = astepo.makedefault(self.fqList)
                    tmpcmd,tmpout = astepo.makeCommand(default['input'])
                    for i in range(len(tmpcmd)):
                        outshell = open(outdir + "/shell/" + step + str(i) +".sh", mode='w')
                        outshell.write(tmpcmd[i])
                        outshell.close()
                    stepdict = json.dumps(astepo.makedefault(self.fqList))
                    out.write("\"%s\":%s,\n" % (step, stepdict))
                    os.system("sed -i \'s/;/\\n/g\' "+ outdir+"/shell/" +step+"*.sh")
            out.write("\"outdir\":\"%s\"\n" % (self.outdirMain))
            out.write("}\n")
            out.close()


        except IOError as e:
            raise e

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

def checkSpecies(species):
    if species is None:
        return {"none":"none"}
    else:
        parta=species.split(',')
        soft="none"
        specie="none"
        sdict={}
        for content in parta:
            if re.search(r'\:',content):
                contents=content.split(":")
                soft=contents[0]
                specie=contents[1]
            else:
                specie=content
            try:
                sdict[soft].append(specie)
            except:
                sdict[soft]=[specie]
        return sdict

def loadFqList(fqList):
    if fqList is None:
        return ["test_1.fq.gz", "test_2.fq.gz"], {'a': ["1", "2", "3", "4"]}
    else:
        try:
            lines = open(fqList, mode='r').readlines()
        except IOError as e:
            raise e
        stat = {}
        fq = []
        for line in lines:
            linep = line.split()
            fq1 = linep[2]
            fq1base = os.path.basename(fq1)
            fq1prefix = re.sub(r'_\d\..*$', r'', fq1base)
            stat[fq1prefix] = linep
            fq += [linep[2], linep[3]]
        return stat, fq

if __name__=="__main__":
    if len(sys.argv) == 1 :
        print("You need to set some parameter")
        sys.exit()
    pwd = os.path.abspath('.')
    parser = argparse.ArgumentParser(description="pipeline annotator help")
    parser.add_argument('--mode', dest='runMode', type=str,help='how to action: run/makejson. run means executing the workflow. makejson means make a json with default parameter with defualt name')
    parser.add_argument('--fqlist',dest='fqList',type=str,help="the list file that could contain five column: sampleID libraryID fq1path fq2path [specieA,specieB].[] means optional. species will used in RNAseq commparison test.")
    parser.add_argument('--genomefa',dest='genomeFa',type=str,help="the genome fa used in workflow.\n HG19:/hwfssz1/BIGDATA_COMPUTING/GaeaProject/reference/hg19/hg19.fasta\n HG38:/hwfssz1/BIGDATA_COMPUTING/GaeaProject/reference/hg38/hg38.fa")
    parser.add_argument('--species',dest='species',type=str,help="set species name where needed. format: augustus:A,genewise:B,C,D,fgene:E,F. A will used in augustus,BCD will used in genwise and so on.")
    parser.add_argument('--outdir',dest='outdir',type=str,default=pwd,help='the output directory,default current directory')
    localeArg=parser.parse_args()
    absoutdir = os.path.abspath(localeArg.outdir)

    os.makedirs(absoutdir, mode=0o755, exist_ok=True)
    allSpecies = checkSpecies(localeArg.species)

    a = interface()
    a.species = allSpecies
    a.ref=localeArg.genomeFa
    fqlist=localeArg.fqList
    a.fqLink,a.fqList=loadFqList(fqlist)
    a.outdirMain=absoutdir

    if localeArg.runMode is None:
        print("need set --mode")
        sys.exit()
    if localeArg.runMode == 'makeshell':
        if localeArg.fqList is None and localeArg.genomeFa is None:
            print("need fqlist or genome to makeshell")
            sys.exit()
        else:
            a.makeshell()
    elif localeArg.runMode == 'run':
        a.makeshell()
        a.runlocal()
    else:
        print("Unknow mode ==> %s" %(localeArg.runMode))
