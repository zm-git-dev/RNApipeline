#!/usr/bin/env python
# -*- coding:utf-8 -*-

import json
import logging
import os
import re
import sys


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
        self.species={}
        self.fqLink = {"GTpi":["GTpi","GTpi","/ldfssz1/ST_BIGDATA/PMO/WORKFLOW/subworkflow_test_raw_data/rnaseq_tests/GTpi_1.fq.gz","/ldfssz1/ST_BIGDATA/PMO/WORKFLOW/subworkflow_test_raw_data/rnaseq_tests/GTpi_2.fq.gz","GTye"]}
        self.outdirMain = os.path.abspath('.')
        self.ref = "/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/GT.genome.fa"
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

class filter(common):

    def __init__(self):
        super(filter, self).__init__()
        self.parameter = "-l 15 -q 0.2 -n 0.05 -i -Q 1 -5 0 -c 0.1 " \
                         "-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA "
        self.soapnuke = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/SOAPnuke"
        self.fqcheck = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/fqcheck"
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/Filter/"
        self.program=[self.soapnuke,self.fqcheck]
        self.outdir = "Transcript_Assembly/Filter_SOAPnuke"


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

        self.parameter="--phred64 --no-discordant --no-mixed -I 1 -X 1000 -p 8 "

        self.samtools="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/samtools"
        self.hisat2="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/hisat2-2.0.4"
        self.java="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/java"
        self.picard="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/picard"
        self.program=[self.samtools,self.hisat2,self.java,self.picard]

        self.outdir = "Transcript_Assembly/GenomeMapping_HISAT/"

    def makeCommand(self, inputfq):
        filter_para = filter()
        filter_para.species=self.species
        filter_para.fqLink = self.fqLink
        filter_para.outdir=self.outdir.replace("GenomeMapping_HISAT","Filter_SOAPnuke")
        CleanDataDict=inputfq[0]

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
                hisat2shell+="cd {sampledir}; {hisat2}/hisat2 {hisat2_para} -x {ref} -1 {fq1} -2 {fq2} " \
                             "2>{outdir}/{sampleid}.Map2GenomeStat.xls |  " \
                             "{samtools} view -b -S -o {sampledir}/{sampleid}.bam - ;" \
                             "{java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {picard}/SortSam.jar I={sampleid}.bam O={sampleid}.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT\n".format(
                    sampledir=sampledir,
                    hisat2=self.hisat2,
                    java=self.java,
                    picard=self.picard,
                    ref=self.ref,
                    hisat2_para=self.parameter,
                    fq1=cleanFqA,fq2=cleanFqB,samtools=self.samtools,
                    outdir=self.outdir,sampleid=SampleID
                )

                BamPath="{outDir}/{sampleid}.Sort.bam" \
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
                hisat2shell += "cd {sampledir}; {hisat2}/hisat2 {hisat2_para} -x {ref} -1 {fq1} " \
                               "2>{outdir}/{sampleid}.Map2GenomeStat.xls |  " \
                               "{samtools} view -b -S -o {sampledir}/{sampleid}.bam - " \
                               "{java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {picard}/SortSam.jar I={sampleid}.bam O={sampleid}.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT\n".format(
                    sampledir=sampledir,
                    hisat2=self.hisat2,
                    ref=self.ref,
                    hisat2_para=self.parameter,
                    fq1=cleanFq, samtools=self.samtools,
                    outdir=self.outdir, sampleid=SampleID
                )

                BamPath = "{outDir}/{sampleid}.Sort.bam" \
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

        input=[]
        output=[]
        BamDict={}

        input.append(CleanDataDict)

        if self.PEorSE == "PE":
            for SampleID,Cleandata in CleanDataDict.items():
                sampledir=self.outdir+"/"+SampleID
                BamPath="{outDir}/{sampleid}.Sort.bam" \
                        "".format(outDir=sampledir,sampleid=SampleID)
                BamDict[SampleID]=BamPath
            output.append(BamDict)
        elif self.PEorSE == "SE":
            for SampleID, Cleandata in CleanDataDict.items():
                sampledir = self.outdir + "/" + SampleID
                BamPath = "{outDir}/{sampleid}.Sort.bam" \
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
        self.outdir = "Transcript_Assembly/NovelTr_Stringtie"

        self.parameter = "-f 0.3 -j 3 -c 5 -g 100 -p 8"

        self.stringtie="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/stringtie"

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/NovelTr"

        self.program = [self.stringtie]

    def makeCommand(self,inputfq):

        alignment_o = alignment()
        alignment_o.fqLink = self.fqLink
        alignment_o.species=self.species
        BamDict=inputfq

        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/Prediction/", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/Prediction/merge/", mode=0o755, exist_ok=True)

        cmd=[]
        output=[]
        reconstruct_shell = ""

        gtflist = open(self.outdir+"/Prediction/"+"gtf.list",'w')
        for SampleID,BamPath in BamDict.items():
            reconstruct_shell += "{stringtie} {bam_path}  -o {outdir}/{sampleid}.stringtie.gtf.tmp {parameter};" \
                                 "perl {scriptbin}/CorrectGTF.pl {ref} {outdir}/{sampleid}.stringtie.gtf.tmp {outdir}/{sampleid}.stringtie.gtf\n" \
                                 "".format(
                stringtie=self.stringtie,
                parameter=self.parameter,
                bam_path=BamPath,
                outdir=self.outdir+"/Prediction/",
                sampleid=SampleID,
                scriptbin=self.scriptbin,
                ref=self.ref
            )
            gtflist.write(self.outdir+"/Prediction/"+SampleID+".stringtie.gtf"+"\n")
        gtflist.close()

        cuffmerge_cpc_shell="{stringtie} --merge -o {outdir}/Prediction/merge/merged.gtf -p 8 {gtf_list}".format(
            stringtie=self.stringtie,
            scriptbin=self.scriptbin,
            outdir=self.outdir,
            gtf_list=self.outdir+"/Prediction/"+"gtf.list",
            ref=self.ref,
        )

        os.makedirs(self.outdir + "/Prediction/merge/", mode=0o755, exist_ok=True)


        cmd.append(reconstruct_shell)
        cmd.append(cuffmerge_cpc_shell)

        output.append(self.outdir+"/Prediction/merge/merged.gtf")

        return  cmd,output

    def makedefault(self,inputfq):

        alignment_o = alignment()
        alignment_o.fqLink = self.fqLink
        alignment_o.species=self.species
        alignment_o.outdir=self.outdir.replace("NovelTr_Stringtie","GenomeMapping_HISAT")
        BamDict = alignment_o.makedefault(inputfq)["output"][0]

        output=[]
        output.append(self.outdir+"/Prediction/merge/merged.gtf")


        default={
            'input':BamDict,
            'parameter':self.parameter,
            'program':self.program,
            'resource':"4G,8CPU",
            'output':output
        }
        return default

class interface(common):

    def __init__(self):
        common.__init__(self)
        self.step = [["filter"],["alignment"],["novel_tr"]]
        self.input = "%s/workflow.json" % (self.outdirMain)
        self.output = "%s/workflow.json" % (self.outdirMain)

    def makeshell(self, outputfile=None):
        outputjson = self.output
        outdir = self.outdirMain
        os.makedirs(outdir+"/shell",exist_ok=True,mode=0o755)
        if outputfile is not None:
            outputjson = outputfile
        try:
            out = open(outputjson, mode='w')
            out.write("{\n")

            for stepL in self.step:
                for step in stepL:
                    outshell = open(outdir+"/shell/"+step+".sh", mode='w')
                    astep = eval(step)
                    astepo = astep()
                    astepo.fqLink = self.fqLink
                    astepo.outdir = outdir +"/"+astepo.outdir
                    default = astepo.makedefault(self.fqList)
                    tmpcmd,tmpout = astepo.makeCommand(default['input'])
                    for i in tmpcmd:
                        outshell.write(i+"\n")
                    stepdict = json.dumps(astepo.makedefault(self.fqList))
                    out.write("\"%s\":%s,\n" % (step, stepdict))
                    outshell.close()
                    os.system("sed -i \'s/;/\\n/g\' "+ outdir+"/shell/" +step+".sh")
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

    a = interface()
    a.makeshell()
