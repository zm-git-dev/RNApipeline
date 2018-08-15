# !/usr/bin/env python
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
@file: RNAdenovo.py 
@time: 2018/05/10 
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
        self.parameter = "-l 15 -q 0.2 -n 0.05 -i -Q 1 -5 0  -c 2 " \
                         "-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA "
        self.soapnuke = self.getsoftware().SOAPNUKE
        self.fqcheck = self.getsoftware().FQCHECK
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/Filter"
        self.program=[self.soapnuke,self.fqcheck]
        self.outdir = "RNAdenovo/Filter_SOAPnuke"


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

class trinity_assemble(common):

    def __init__(self):
        super(trinity_assemble,self).__init__()
        soft = self.getsoftware()
        self.program=soft.TRINITY
        self.parameter="--seqType fq --max_memory 10G --min_contig_length 250  --CPU 8 --min_kmer_cov 3 --min_glue 3 --bfly_opts '-V 5 --edge-thr=0.1 --stderr'"
        self.outdir="RNAdenovo/Denovo_Trinity"
        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/Denovo"

    def makeCommand(self,inputfq):
        CleanDataDict=inputfq[0]
        os.makedirs(self.outdir,mode=0o755, exist_ok=True)
        cmd=""
        output=[]
        fq1_list=[]
        fq2_list=[]
        for SampleID, Cleandata in CleanDataDict.items():
            cleanFqA = Cleandata["clean_fq1"]
            cleanFqB = Cleandata["clean_fq2"]
            fq1_list.append(cleanFqA)
            fq2_list.append(cleanFqB)
        leftfq=",".join(fq1_list)
        rightfq=",".join(fq2_list)

        cmd+="cd {outdir};" \
             "{trinity} {parameter} --left {fq1_list} --right {fq2_list} --output {outdir}/Trinity; " \
             "ln -s {outdir}/Trinity/Trinity.fasta Unigene.fa; " \
             "/ldfssz1/ST_BIGDATA/USER/yueyao/bin/miniconda2/bin/python {scriptbin}/extract_longest_isform.py Unigene.fa All ;" \
             "perl {scriptbin}/fishInWinter.pl -bf table -ff fasta All_Unigene_id.txt Unigene.fa >Temp.Unigene.fa; " \
             "perl -lane 'if(/(>.*c\d+_g\d+)(_i\d+)\slen/){{print $1}}else{{print}}' Temp.Unigene.fa >All-Unigene.fa; " \
             "perl {scriptbin}/get_Trinity_gene_to_trans_map.pl Unigene.fa >All-Unigene.gene2mark;" \
             "rm Temp.Unigene.fa;".format(

            trinity=self.program,
            parameter=self.parameter,
            fq1_list = leftfq,
            fq2_list = rightfq,
            outdir = self.outdir,
            scriptbin=self.scriptbin
        )
        output.append(self.outdir+"/All-Unigene.fa")
        output.append(self.outdir+"/Unigene.fa")
        output.append(self.outdir+"/All-Unigene.gene2mark")
        return [cmd],output

    def makedefault(self,inputfq):
        filter_o = filter()
        filter_o.species=self.species
        filter_o.fqLink=self.fqLink
        filter_o.outdir = self.outdir.replace("Denovo_Trinity", "Filter_SOAPnuke")
        CleanDataDict=filter_o.makedefault(inputfq)["output"]
        output=[]
        output.append(self.outdir+"/All-Unigene.fa")
        output.append(self.outdir+"/Unigene.fa")
        output.append(self.outdir+"/All-Unigene.gene2mark")

        default = {
            'input': CleanDataDict,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "10G,8CPU",
            'output': output
        }
        return default

class annotation(common):
    def __init__(self):
        super(annotation,self).__init__()
        self.outdir="RNAdenovo/Annotation_Blast"
        self.parameter={
            "Annotation_Method" :"Diamond",
            "Annotation_Options":"--evalue 1e-5",
            "Annotation_split":2,
            "Annotation_db":"nr,nt,swissprot,go,kegg,kog",
            "Annotation_dbClass":"pl"
        }

        self.diamond=self.getsoftware().DIAMOND
        self.blastp=self.getsoftware().BLASTP
        self.blastn=self.getsoftware().BLASTN
        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/Annotation"
        self.program = [self.diamond, self.blastn, self.blastp]

    def makeCommand(self,inputfq):
        cmd=[]
        output=[]
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        Unigene=inputfq

        split_fasta_path=self.outdir+'/split_fasta'
        os.makedirs(split_fasta_path,mode=0o755, exist_ok=True)
        split_cmd="{scriptbin}/fastaDeal.pl --cutf {num} --outdir {split_path} {unigene}".format(
            scriptbin=self.scriptbin,
            num=self.parameter["Annotation_split"],
            split_path=split_fasta_path,
            unigene=Unigene
        )
        cmd.append(split_cmd)

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)

        db_list=self.parameter["Annotation_db"].split(",")
        num=self.parameter["Annotation_split"]+1
        blast_cmd=""
        blast_process_cmd=""

        stat_cmd =""
        stat_cmd += "grep '^>' {Unigene} |sed 's/>//' |awk '{{print $1}}' > {outdir}/all_gene.id;" \
                    "perl {scriptbin}/all_stat.pl -list {outdir}/all_gene.id ".format(
            Unigene=Unigene,
            scriptbin=self.scriptbin,
            outdir=self.outdir
        )


        for db_val in db_list:
            if db_val == "nt":
                nt_path = self.outdir + '/nt'
                os.makedirs(nt_path, mode=0o755, exist_ok=True)
                blast_cmd += "{blastn} -dust no {parameter} -num_threads 5 -db {nt_db} -query {unigene} " \
                         "-out {outdir}/All-Unigene.fa.blast.nt -outfmt 5 \n".format(
                    blastn=self.blastn,
                    parameter=self.parameter["Annotation_Options"].replace("-","",1),
                    nt_db=dbclass.NT,
                    unigene=Unigene,
                    outdir=nt_path
                )

                blast_process_cmd+="perl {scriptbin}/blast_m7_parser.pl {outdir}/nt/All-Unigene.fa.blast.nt {outdir}/All-Unigene.fa.blast.nt.temp;" \
                                   "head -n1 {outdir}/All-Unigene.fa.blast.nt.temp >{outdir}/All-Unigene.fa.blast.nt.xls;" \
                                   "/usr/bin/python {scriptbin}/filter_top5nt.py {outdir}/All-Unigene.fa.blast.nt.temp >>{outdir}/All-Unigene.fa.blast.nt.xls;" \
                                   "rm {outdir}/All-Unigene.fa.blast.nt.temp\n".format(
                    outdir=self.outdir,
                    scriptbin=self.scriptbin
                )
                output.append(self.outdir+"/All-Unigene.fa.blast.nt.xls")
                stat_cmd += " -nt "+self.outdir+"/All-Unigene.fa.blast.nt.xls "
            if db_val == "nr":
                nr_path = self.outdir + '/nr'
                os.makedirs(nr_path, mode=0o755, exist_ok=True)
                for i in range(1, num):
                    blast_cmd+="{diamond} blastx {parameter} -d {nr_db} -q {split_fa} " \
                             "-o {nr_outdir}/All-Unigene.fa.{num}.blast.nr --threads 5 --outfmt  5 " \
                             "--seg no   --max-target-seqs 5 --more-sensitive -b 0.2 --salltitles;" \
                             "cat {nr_outdir}/All-Unigene.fa.{num}.blast.nr >>{outdir}/All-Unigene.fa.blast.nr\n".format(

                        diamond=self.diamond,
                        parameter=self.parameter["Annotation_Options"],
                        nr_db=dbclass.NR,
                        split_fa=split_fasta_path+"/All-Unigene.fa."+str(i),
                        nr_outdir=nr_path,
                        outdir=self.outdir,
                        num=i

                    )

                blast_process_cmd+="perl {scriptbin}/blast_m7_parser.pl {outdir}/All-Unigene.fa.blast.nr {outdir}/All-Unigene.fa.blast.nr.xls;" \
                                   "perl {scriptbin}/blast_nr_class.pl -nr {outdir}/All-Unigene.fa.blast.nr.xls -outdir {outdir}\n".format(
                    outdir=self.outdir,
                    scriptbin=self.scriptbin
                )
                output.append(self.outdir + "/All-Unigene.fa.blast.nr.xls")
                stat_cmd += " -nr " + self.outdir + "/All-Unigene.fa.blast.nr.xls "
            if db_val == "swissprot":
                swissprot_path=self.outdir+'/swissport'
                os.makedirs(swissprot_path,mode=0o755,exist_ok=True)
                for i in range(1, num):
                    blast_cmd+="{diamond} blastx {parameter} -d {nr_db} -q {split_fa} " \
                             "-o {outdir}/All-Unigene.fa.{num}.blast.swissprot --threads 5 --outfmt  6 " \
                             "--seg no   --max-target-seqs 5 --more-sensitive -b 0.2 --salltitles\n".format(

                        diamond=self.diamond,
                        parameter=self.parameter["Annotation_Options"],
                        nr_db=dbclass.SWISSPROT,
                        split_fa=split_fasta_path+"/All-Unigene.fa."+str(i),
                        outdir=swissprot_path,
                        num=i
                    )

                blast_process_cmd+="cat {swissprot_path}/*blast.swissprot > {outdir}/All-Unigene.fa.blast.swissprot;" \
                                   "perl {scriptbin}/get_annot_info.pl -tophit 5 -topmatch 1 " \
                                   "-id {swissprot_id} -input " \
                                   "{outdir}/All-Unigene.fa.blast.swissprot -out {outdir}/All-Unigene.fa.blast.swissprot.xls\n".format(
                    swissprot_path=swissprot_path,
                    scriptbin=self.scriptbin,
                    swissprot_id=dbclass.SWISSPROT_GENEID,
                    outdir=self.outdir
                )
                output.append(self.outdir + "/All-Unigene.fa.blast.swissprot.xls")
                stat_cmd += " -swissprot " + self.outdir + "/All-Unigene.fa.blast.swissprot.xls "
            if db_val == "kegg":
                kegg_path=self.outdir+'/kegg'
                os.makedirs(kegg_path,mode=0o755, exist_ok=True)
                os.makedirs(self.outdir+'/All-Unigene.fa_map',mode=0o755, exist_ok=True)
                for i in range(1, num):
                    blast_cmd+="{diamond} blastx {parameter} -d {nr_db} -q {split_fa} " \
                             "-o {outdir}/All-Unigene.fa.{num}.blast.kegg --threads 5 --outfmt  6 " \
                             "--seg no   --max-target-seqs 5 --more-sensitive -b 0.2 --salltitles\n".format(

                        diamond=self.diamond,
                        parameter=self.parameter["Annotation_Options"],
                        nr_db=dbclass.KEGG,
                        split_fa=split_fasta_path+"/All-Unigene.fa."+str(i),
                        outdir=kegg_path,
                        num=i

                    )

                blast_process_cmd+="cat {kegg_path}/*blast.kegg > {outdir}/All-Unigene.fa.blast.kegg;" \
                                   "perl {scriptbin}/get_annot_info.pl -tophit 5 -topmatch 1 -id {kegg_id} -input " \
                                   "{outdir}/All-Unigene.fa.blast.kegg -out {outdir}/All-Unigene.fa.blast.kegg.xls;" \
                                   "perl {scriptbin}/blast2ko.pl -input {Unigene} -output {outdir}/All-Unigene.fa.ko -blastout {outdir}/All-Unigene.fa.blast.kegg -kegg {kegg_db};" \
                                   "perl {scriptbin}/pathfind.pl -kegg {kegg_db} -maptitle {maptitle} -komap {komap} -fg {outdir}/All-Unigene.fa.ko -output {outdir}/All-Unigene.fa.path;" \
                                   "perl {scriptbin}/keggMap_nodiff.pl -komap {komap} -mapdir {mapdir} -ko {outdir}/All-Unigene.fa.ko -outdir {outdir}/All-Unigene.fa_map;" \
                                   "perl {scriptbin}/genPathHTML.pl -indir {outdir};" \
                                   "perl {scriptbin}/unigene.drawKEGG.pl -path {outdir}/All-Unigene.fa.path -outprefix {outdir}/All-Unigene.fa -idCol 3 -level1Col 4 -level2Col 5 -geneCol 6\n".format(
                    scriptbin=self.scriptbin,
                    outdir=self.outdir,
                    kegg_id=dbclass.KEGG_GENEID,
                    Unigene=Unigene,
                    kegg_db=dbclass.KEGG,
                    maptitle=dbclass.KEGG_MAP_TITLE,
                    komap=dbclass.KEGG_KOMAP,
                    mapdir=dbclass.KEGG_MAP_DIR,
                    kegg_path=kegg_path
                )
                output.append(self.outdir + "/All-Unigene.fa.blast.kegg.xls")
                stat_cmd += " -kegg " + self.outdir + "/All-Unigene.fa.blast.kegg.xls "
            if db_val == "kog":
                kog_path=self.outdir+'/kog'
                os.makedirs(kog_path, mode=0o755, exist_ok=True)
                for i in range(1, num):
                    blast_cmd+="{diamond} blastx {parameter} -d {nr_db} -q {split_fa} " \
                             "-o {outdir}/All-Unigene.fa.{num}.blast.kog --threads 5 --outfmt  6 " \
                             "--seg no --max-target-seqs 5 --more-sensitive -b 0.2 --salltitles\n".format(

                        diamond=self.diamond,
                        parameter=self.parameter["Annotation_Options"],
                        nr_db=dbclass.KOG_ALL,
                        split_fa=split_fasta_path+"/All-Unigene.fa."+str(i),
                        outdir=kog_path,
                        num=i

                    )
                blast_process_cmd +="cat {kog_path}/*blast.kog > {outdir}/All-Unigene.fa.blast.kog;" \
                                    "perl {scriptbin}/get_annot_info.pl -tophit 5 -topmatch 1 -id {kog_id} -input {outdir}/All-Unigene.fa.blast.kog -out {outdir}/All-Unigene.fa.blast.kog.xls;" \
                                    "perl {scriptbin}/cog_parser.pl {merge} {fun} {outdir}/All-Unigene.fa.blast.kog.xls ;" \
                                    "perl {scriptbin}/cog_R.pl -catalog {outdir}/All-Unigene.fa.KOG2Gene.xls -sample All-Unigene.fa -outdir {outdir}\n".format(
                    kog_path=kog_path,
                    scriptbin=self.scriptbin,
                    outdir=self.outdir,
                    merge=dbclass.KOG_WHOG,
                    fun=dbclass.KOG_FUN,
                    kog_id=dbclass.KOG_GENEID
                )
                output.append(self.outdir + "/All-Unigene.fa.blast.kog.xls")
                stat_cmd += " -cog " + self.outdir + "/All-Unigene.fa.blast.kog.xls "
            if db_val == "cog":
                cog_path=self.outdir+'/cog'
                os.makedirs(cog_path, mode=0o755, exist_ok=True)
                for i in range(1, num):
                    blast_cmd+="{diamond} blastx {parameter} -d {cog_db} -q {split_fa} " \
                             "-o {outdir}/All-Unigene.fa.{num}.blast.cog --threads 5 --outfmt  6 " \
                             "--seg no   --max-target-seqs 5 --more-sensitive -b 0.2 --salltitles\n".format(

                        diamond=self.diamond,
                        parameter=self.parameter["Annotation_Options"],
                        cog_db=dbclass.COG_ALL,
                        split_fa=split_fasta_path+"/All-Unigene.fa."+str(i),
                        outdir=cog_path,
                        num=i

                    )
                blast_process_cmd +="cat {kog_path}/*blast.cog > {outdir}/All-Unigene.fa.blast.cog;" \
                                    "perl {scriptbin}/get_annot_info.pl -tophit 5 -topmatch 1 -id {cog_id} -input {outdir}/All-Unigene.fa.blast.cog -out {outdir}/All-Unigene.fa.blast.cog.xls;" \
                                    "perl {scriptbin}/cog_parser.pl {merge} {fun} {outdir}/All-Unigene.fa.blast.cog.xls ;" \
                                    "perl {scriptbin}/cog_R.pl -catalog {outdir}/All-Unigene.fa.COG2Gene.xls -sample All-Unigene.fa -outdir {outdir}\n".format(
                    kog_path=cog_path,
                    scriptbin=self.scriptbin,
                    outdir=self.outdir,
                    merge=dbclass.COG_WHOG,
                    fun=dbclass.COG_FUN,
                    cog_id=dbclass.COG_GENEID
                )
                output.append(self.outdir + "/All-Unigene.fa.blast.cog.xls")
                stat_cmd += " -cog " + self.outdir + "/All-Unigene.fa.blast.cog.xls "
            if db_val == "go":
                go_path=self.outdir+'/go'
                os.makedirs(go_path, mode=0o755, exist_ok=True)
                os.makedirs(self.outdir+'/rna_data', mode=0o755, exist_ok=True)
                blast_process_cmd+="perl {scriptbin}/blast_m7_parser.pl {outdir}/All-Unigene.fa.blast.nr {go_outdir}/All-Unigene.fa.blast.nr.tab; " \
                       "python {scriptbin}/BLAST3GO.py --NR_tab {go_outdir}/All-Unigene.fa.blast.nr.tab --Accession_go {accssion_go} " \
                       " --output {go_outdir}/All-Unigene.fa.blast.nr.annot; " \
                       "perl {scriptbin}/annot2goa.pl {obo} {go_outdir}/All-Unigene.fa.blast.nr.annot {outdir}/rna_data/species;" \
                       "perl {scriptbin}/blast_m7_m8.pl -input {nr_outdir}/*blast.nr -output {outdir}/rna_data/species.nr.m8; " \
                       "perl {scriptbin}/getNrDesc.pl -input {outdir}/rna_data/species.nr.m8  -rank 1 -nr {nr_db} {outdir}/rna_data/species.nr.desc;" \
                       "cat {outdir}/rna_data/species.[CFP] | cut -f2 |sort -u >{outdir}/gene.list; " \
                       "perl {scriptbin}/drawGO.pl -list {outdir}/gene.list  -goprefix {outdir}/rna_data/species -goclass {goclass} -outprefix " \
                       "{outdir}/All-Unigene.fa\n".format(
                    scriptbin=self.scriptbin,
                    outdir=self.outdir,
                    go_outdir=go_path,
                    nr_db=dbclass.NR,
                    obo = dbclass.OBO,
                    accssion_go=dbclass.ACCESSION2GO,
                    goclass=dbclass.GOCLASS,
                    nr_outdir=self.outdir+"/nr"
                )
                output.append(self.outdir + "/All-Unigene.fa.Gene2GO.xls")
                stat_cmd += " -go " + self.outdir + "/All-Unigene.fa.Gene2GO.xls -obo " + dbclass.OBO

        tmp_infile=",".join(output)
        tmp_name=",".join(db_list)

        os.makedirs(self.outdir+"/Venny", mode=0o755, exist_ok=True)
        stat_cmd += " -outxls {outdir}/annotation.xls -outstat {outdir}/annotation_stat.xls;" \
                    "perl {scriptbin}/venny.pl -infile {infile} -name {name} -header -color -outdir " \
                    "{outdir}/Venny -imgname Annotation.venn;" \
                    "perl {scriptbin}/addDesc.pl -input \" {outdir}/Venny/*.xls\" -annot {outdir}/annotation.xls;" \
                    "cp {outdir}/All-Unigene.fa.ko {outdir}/rna_data/species.ko".format(
            outdir=self.outdir,
            scriptbin=self.scriptbin,
            infile = tmp_infile,
            name = tmp_name
        )

        cmd.append(blast_cmd)
        cmd.append(blast_process_cmd)
        cmd.append(stat_cmd)
        output.append(self.outdir+"/annotation.xls")
        output.append(self.outdir + "/annotation_stat.xls")
        output.append(self.outdir+ "/rna_data/species.ko")
        output.append(self.outdir + "/species.C")
        return cmd,output

    def makedefault(self,inputfq):

        trinity_o = trinity_assemble()
        trinity_o.species = self.species
        trinity_o.outdir=self.outdir.replace("Annotation_Blast", "Denovo_Trinity")

        spe = trinity_o.species["RNAdenovo"][0]
        self.parameter["Annotation_dbClass"] = spe

        Unigene = trinity_o.makedefault(inputfq)["output"][0]
        output=[]
        output.append(self.outdir+"/annotation.xls")
        output.append(self.outdir + "/annotation_stat.xls")
        output.append(self.outdir+ "/rna_data/species.ko")
        output.append(self.outdir + "/rna_data/species.C")

        default={
            'input': Unigene,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "5G,5CPU",
            'output': output
        }

        return default

class ssr(common):
    def __init__(self):
        super(ssr,self).__init__()
        self.outdir="RNAdenovo/SSR_MISA"
        self.parameter="1-12,2-6,3-5,4-5,5-4,6-4 100 150"
        self.ssr=self.getsoftware().SSR
        self.blastall=self.getsoftware().BLASTALL
        self.formatdb=self.getsoftware().FORMATDB
        self.program = [self.ssr, self.blastall, self.formatdb]


    def makeCommand(self,inputfq):
        output=[]
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)

        Unigene=inputfq
        cmd="cd {0} && ln -s {1} ./All-Unigene.fa;".format(self.outdir,Unigene)
        cmd+="{0}/misa.pl All-Unigene.fa {1} All-Unigene {2};".format(self.ssr,self.outdir,self.parameter)
        cmd+="{0}/2primer_designer.pl All-Unigene.ssr.out.xls All-Unigene.raw_primer.out.xls All-Unigene.primer_results.out.xls; ".format(self.ssr)
        cmd+="{0}/3filter_primer.pl All-Unigene.primer_results.out.xls All-Unigene.rescreen.out.xls All-Unigene.blastin; ".format(self.ssr)
        cmd+="{0} -i All-Unigene.fa -p F; ".format(self.formatdb)
        cmd += "{0} -i All-Unigene.blastin -d All-Unigene.fa -p blastn -F F -a 15 -o All-Unigene.blast.out.xls; ".format(self.blastall)
        cmd += "{0}/5blast_parse.pl All-Unigene.blast.out.xls All-Unigene.query_sbjct.out.xls All-Unigene.stactic.out.xls 4; ".format(self.ssr)
        cmd += "{0}/6generate_primer_pair_new.pl All-Unigene.query_sbjct.out.xls All-Unigene.primer.tab.xls; ".format(self.ssr)
        cmd += "{0}/7final_pair.pl All-Unigene.primer.tab.xls All-Unigene.rescreen.out.xls All-Unigene.inter_primer.tab.xls; ".format(self.ssr)
        cmd += "{0}/8product_ssr_check.pl All-Unigene.fa All-Unigene.inter_primer.tab.xls All-Unigene.product_file.out.xls; ".format(self.ssr)
        cmd += "{0}/ssr_finder.pl --ssr 12 All-Unigene.product_file.out.xls All-Unigene.product_ssr.out.xls; ".format(self.ssr)
        cmd += "{0}/9filter_ssr.pl All-Unigene.product_ssr.out.xls All-Unigene.inter_primer.tab.xls All-Unigene.only_primer.tab.xls; ".format(self.ssr)
        cmd +="{0}/10get_fitprimer.pl All-Unigene.only_primer.tab.xls All-Unigene.product_ssr.out.xls All-Unigene.rescreen.out.xls All-Unigene.final_primer.out; ".format(self.ssr)
        cmd += "sort -k 1,1 All-Unigene.final_primer.out >All-Unigene.final_primer.out.xls; "
        cmd += "{0}/stat_draw.pl All-Unigene.statistics.xls {1} ".format(self.ssr,self.outdir)

        output.append(self.outdir+"/All-Unigene.final_primer.out.xls")
        return [cmd],output

    def makedefault(self,inputfq):
        trinity_o = trinity_assemble()
        trinity_o.species=self.species
        trinity_o.fqLink=self.fqLink
        trinity_o.outdir = self.outdir.replace("SSR_MISA", "Denovo_Trinity")

        Unigene = trinity_o.makedefault(inputfq)["output"][0]
        output=[]
        output.append(self.outdir + "/All-Unigene.final_primer.out.xls")
        default={
            'input': Unigene,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,2CPU",
            'output': output
        }

        return default

class cds_predict(common):
    def __init__(self):
        super(cds_predict,self).__init__()
        self.outdir="RNAdenovo/CDSpredict_TransDecoder"
        self.diamond=self.getsoftware().DIAMOND
        self.hmmscan=self.getsoftware().HMMSCAN
        self.transdecoder=self.getsoftware().TRANSDECODER
        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/CDSpredict"
        self.parameter={
            "Annotation_Options":"--evalue  1e-5  --threads 5  --seg no --max-target-seqs 1 --more-sensitive -b 0.5 --salltitles ",
            "Annotation_dbClass":"pl"
        }
        self.program=[self.diamond,self.hmmscan]

    def makeCommand(self,inputfq):
        output=[]
        cmd=[]
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        Unigene=inputfq
        cds_shell="cd {0};/usr/bin/perl {1}/TransDecoder.LongOrfs -t {2}; ".format(self.outdir,self.transdecoder,Unigene)
        self.parameter["Annotation_dbClass"]=self.species["RNAdenovo"][0]

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
            cds_shell+="{diamond} blastp -d {db} -q {outdir}/All-Unigene.fa.transdecoder_dir/longest_orfs.pep {para} -o " \
                 "{outdir}/blastp.outfmt6;".format(
                    diamond=self.diamond,
                    db=dbclass.SWISSPROT,
                    para=self.parameter["Annotation_Options"],
                    outdir=self.outdir
                )
            cds_shell+="{hmmscan} --cpu 5  --domtblout pfam.domtblout {pfam} {outdir}/All-Unigene.fa.transdecoder_dir/longest_orfs.pep;".format(
                hmmscan=self.hmmscan,
                pfam=dbclass.PFAM,
                outdir=self.outdir
            )
            cds_shell+="/usr/bin/perl {program}/TransDecoder.Predict -t {Unigene} --retain_pfam_hits {outdir}/pfam.domtblout --retain_blastp_hits {outdir}/blastp.outfmt6;".format(
                program=self.transdecoder,
                Unigene=Unigene,
                outdir=self.outdir
            )

            cds_shell +="perl {0}/fa_quality.pl -len -Head -N -gc {1}/All-Unigene.fa.transdecoder.cds;".format(self.scriptbin,self.outdir)
            cds_shell +="perl {0}/barplot.pl {1}/All-Unigene.fa.transdecoder.cds.quality.xls  All-Unigene.fa.cds;".format(self.scriptbin,self.outdir)
            cds_shell +="perl {scriptbin}/Fasta_stat.pl {outdir}/All-Unigene.fa.transdecoder.cds 2>{outdir}/PredictSummary.xls".format(
                scriptbin=self.scriptbin,
                outdir=self.outdir
            )
            output.append(self.outdir+"/All-Unigene.fa.transdecoder.cds")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
            cds_shell+="{diamond} blastp -d {db} -q {outdir}/All-Unigene.fa.transdecoder_dir/longest_orfs.pep {para} -o " \
                 "{outdir}/blastp.outfmt6;".format(
                    diamond=self.diamond,
                    db=dbclass.SWISSPROT,
                    para=self.parameter["Annotation_Options"],
                    outdir=self.outdir
                )
            cds_shell+="{hmmscan} --cpu 5  --domtblout pfam.domtblout {pfam} {outdir}/All-Unigene.fa.transdecoder_dir/longest_orfs.pep;".format(
                hmmscan=self.hmmscan,
                pfam=dbclass.PFAM,
                outdir=self.outdir
            )
            cds_shell+="/usr/bin/perl {program}/TransDecoder.Predict -t {Unigene} --retain_pfam_hits {outdir}/pfam.domtblout --retain_blastp_hits {outdir}/blastp.outfmt6;".format(
                program=self.transdecoder,
                Unigene=Unigene,
                outdir=self.outdir
            )

            cds_shell +="perl {0}/fa_quality.pl -len -Head -N -gc {1}/All-Unigene.fa.transdecoder.cds;".format(self.scriptbin,self.outdir)
            cds_shell +="perl {0}/barplot.pl {1}/All-Unigene.fa.transdecoder.cds.quality.xls  All-Unigene.fa.cds;".format(self.scriptbin,self.outdir)
            cds_shell +="perl {scriptbin}/Fasta_stat.pl {outdir}/All-Unigene.fa.transdecoder.cds 2>{outdir}/PredictSummary.xls".format(
                scriptbin=self.scriptbin,
                outdir=self.outdir
            )
            output.append(self.outdir+"/All-Unigene.fa.transdecoder.cds")
        elif self.parameter["Annotation_dbClass"] == "fg":
            pass
        else:
            sys.exit(1)
        cmd.append(cds_shell)
        return cmd,output

    def makedefault(self,inputfq):
        trinity_o=trinity_assemble()
        trinity_o.species=self.species
        trinity_o.outdir=self.outdir.replace("CDSpredict_TransDecoder","Denovo_Trinity")
        Unigene = trinity_o.makedefault(inputfq)["output"][0]
        self.parameter["Annotation_dbClass"] = trinity_o.species["RNAdenovo"][0]
        output=[]
        output.append(self.outdir + "/All-Unigene.fa.transdecoder.cds")
        default={
            'input': Unigene,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,5CPU",
            'output': output
        }

        return default

class geneexp(common):

    def __init__(self):
        super(geneexp,self).__init__()
        self.outdir= "RNAdenovo/GeneExp"
        self.parameter = "-q --phred64 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant  -p 1 -k 200"
        self.rsempreparereference=self.getsoftware().RSEM+"rsem-prepare-reference"
        self.rsem_calculate_expression=self.getsoftware().RSEM+"rsem-calculate-expression"
        self.samtools = self.getsoftware().SAMTOOLS
        self.bowtie2 = self.getsoftware().BOWTIE2
        self.scriptbin="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/GeneExp"
        self.program = [self.samtools,self.bowtie2,self.rsempreparereference]

    def makeCommand(self,inputfq):
        os.makedirs(self.outdir, mode=0o755, exist_ok=True)

        cmd=[]
        output=[]
        BamDict = {}
        ExpDict = {}
        rsemshell=""

        refmrna = inputfq[0]
        gene2mark= inputfq[1]
        CleanDataDict=inputfq[2]

        stat_shell=""

        buildindexshell = "cd {outdir};{rsem_prepare_reference} {gene_fa} {outdir}/refMrna.fa --bowtie2 --bowtie2-path " \
                          "{bowtie2} --transcript-to-gene-map {gene2tr}".format(
            outdir=self.outdir,
            rsem_prepare_reference=self.rsempreparereference,
            gene_fa=refmrna,
            bowtie2=self.bowtie2,
            gene2tr=gene2mark
        )
        cmd.append(buildindexshell)
        os.makedirs(self.outdir + "/GeneExpression", mode=0o755, exist_ok=True)
        if self.PEorSE == "PE":
            for SampleID, Cleandata in CleanDataDict.items():
                os.makedirs(self.outdir + '/' + SampleID, mode=0o755, exist_ok=True)
                cleanFqA = Cleandata["clean_fq1"]
                cleanFqB = Cleandata["clean_fq2"]
                rsemshell += "{bowtie2} {bowtie2_para} -x {gene_index} -1 {fq1} -2 {fq2} | " \
                            "{samtools} view -S -b -o {outdir}/{sampleid}.bam - ;" \
                             "perl {scriptbin}/BowtieMapStat.pl -bam {outdir}/{sampleid}.bam -key {outdir}/{sampleid}.Bowtie2Gene -seqType PE -samtools {samtools} -gene2tr {gene2tr};" \
                            "{rsem_calculate_expression} --paired-end -p 8 --bam {bampath} {refdir}/refMrna.fa {outdir}/{sampleid} ; " \
                            "awk '{{if($7!=0.00) print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}' " \
                            "{outdir}/{sampleid}.genes.results |grep -v '^ERCC' > {outdir}/{sampleid}.gene.fpkm.xls;" \
                             "cp {outdir}/{sampleid}.gene.fpkm.xls {refdir}/GeneExpression\n".format(
                    bowtie2=self.bowtie2+"bowtie2",
                    scriptbin=self.scriptbin,
                    gene2tr=gene2mark,
                    bowtie2_para=self.parameter,
                    gene_index=self.outdir + "/refMrna.fa",
                    fq1=cleanFqA, fq2=cleanFqB,
                    samtools=self.samtools,
                    bampath=self.outdir + '/' + SampleID +"/"+SampleID + '.bam',
                    rsem_calculate_expression=self.rsem_calculate_expression,
                    refdir=self.outdir,
                    outdir=self.outdir + '/' + SampleID, sampleid=SampleID
                )
                BamPath = "{outDir}/{sampleid}.bam" \
                          "".format(outDir=self.outdir + '/' + SampleID, sampleid=SampleID)
                BamDict[SampleID] = BamPath
                ExpDict[SampleID] = self.outdir + '/' + SampleID + '/' + SampleID + '.gene.fpkm.xls'
        elif self.PEorSE == "SE":
            for SampleID, Cleandata in CleanDataDict.items():
                cleanFq = Cleandata["clean_fq1"]
                rsemshell = "{bowtie2} {bowtie2_para} -x {gene_index} -1 {fq1} 2>{outdir}/{sampleid}.Map2GeneStat.xls | " \
                            "{samtools} view -S -b -o {outdir}/{sampleid}.bam - ;" \
                            "perl {scriptbin}/BowtieMapStat.pl -bam {outdir}/{sampleid}.bam -key {outdir}/{sampleid}.Bowtie2Gene -seqType SE -samtools {samtools} -gene2tr {gene2tr};" \
                           "{rsem_calculate_expression} --paired-end -p 8 --bam {bampath} {refdir}/refMrna.fa {outdir}/{sampleid}; " \
                            "awk '{{if($7!=0.00) print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}\' " \
                            "{outdir}/{sampleid}.genes.results |grep -v '^ERCC' > {outdir}/{sampleid}.gene.fpkm.xls\n".format(
                    bowtie2=self.bowtie2,
                    scriptbin=self.scriptbin,
                    gene2tr=gene2mark,
                    bowtie2_para=self.parameter,
                    gene_index=self.outdir + "/refMrna.fa",
                    fq1=cleanFq,
                    samtools=self.samtools,
                    bampath=self.outdir + '/' + SampleID +"/"+ SampleID + '.bam',
                    rsem_calculate_expression="",
                    refdir=self.outdir,
                    outdir=self.outdir + '/' + SampleID, sampleid=SampleID
                )
                BamPath = "{outDir}/{sampleid}.bam" \
                          "".format(outDir=self.outdir+'/' + SampleID, sampleid=SampleID)
                BamDict[SampleID] = BamPath
                cmd.append(rsemshell)
                ExpDict[SampleID] = self.outdir + '/' + SampleID + '.gene.fpkm.xls'
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)

        stat_shell+="perl {scriptbin}/AllGeneStat.pl {outdir} {outdir}/All.GeneExpression.FPKM.xls;" \
                   "perl {scriptbin}/drawstacked_bar.pl -indir {outdir}/GeneExpression -outdir {outdir}/GeneExpression/;" \
                   "perl {scriptbin}/drawbox_density-plot.pl -i {outdir}/GeneExpression -o {outdir}/GeneExpression;" \
                   "perl {scriptbin}/AllMapStat.pl {outdir} {outdir}/MappingSummary.xls".format(
            scriptbin=self.scriptbin,
            outdir=self.outdir

        )
        all_gene_exp=self.outdir+"/All.GeneExpression.FPKM.xls"
        cmd.append(rsemshell)
        cmd.append(stat_shell)
        output.append(BamDict)
        output.append(ExpDict)
        output.append(all_gene_exp)
        return cmd,output

    def makedefault(self,inputfq):
        trinity_o = trinity_assemble()
        trinity_o.outdir=self.outdir.replace("GeneExp","Denovo_Trinity")
        refmrna = trinity_o.makedefault(inputfq)["output"][1]
        gene2mark= trinity_o.makedefault(inputfq)["output"][2]

        filter_o = filter()
        filter_o.outdir=self.outdir.replace("GeneExp","Filter_SOAPnuke")
        filter_o.fqLink=self.fqLink
        filter_o.species=self.species
        CleanDataDict=filter_o.makedefault(inputfq)["output"][0]

        input=[]
        input.append(refmrna)
        input.append(gene2mark)
        input.append(CleanDataDict)

        output=[]
        BamDict = {}
        ExpDict = {}

        if self.PEorSE == "PE":
            for SampleID, Cleandata in CleanDataDict.items():
                BamPath = "{outDir}/{sampleid}.bam" \
                          "".format(outDir=self.outdir + '/' + SampleID, sampleid=SampleID)
                BamDict[SampleID] = BamPath
                ExpDict[SampleID] = self.outdir + '/' + SampleID + '/' + SampleID + '.gene.fpkm.xls'
        elif self.PEorSE == "SE":
            for SampleID, Cleandata in CleanDataDict.items():
                BamPath = "{outDir}/{sampleid}.bam" \
                          "".format(outDir=self.outdir, sampleid=SampleID)
                BamDict[SampleID] = BamPath
                ExpDict[SampleID] = self.outdir + '/' + SampleID + '.gene.fpkm.xls'
        else:
            logging.warning("Your Library Fastq File maybe some ERROR!")
            sys.exit(1)

        output.append(BamDict)
        output.append(ExpDict)
        output.append(self.outdir+"/All.GeneExpression.FPKM.xls")

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
        self.parameter = {}
        self.program = "DEGseq,DEseq2,EBseq,NOIseq,PossionDis"
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/GeneDiffExp"
        self.outdir = "RNAdenovo/GeneDiffExp_Allin"

    def makeCommand(self, inputfq):
        ExpDict = inputfq

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

                degshell += "perl {scriptbin}/DEGseq.pl  -list {CompareList} -diff {diffcom} -GeneIDColumn 1 -GeneExpColumn 4 " \
                           "-GeneLenColumn 3 -method MARS -threshold 5 -pValue 1e-3 " \
                           "-zScore 4 {deg_para} -outdir {outdir}\n".format(
                    scriptbin=self.scriptbin,
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
                degshell += "perl {GeneDiffExpBin}/EBseq.pl -list {Samplelist} -diff {CompareList} -group {Grouplist} " \
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
                degshell += "{GeneDiffExpBin}/NOIseq.pl -list {Samplelist} -diff {CompareList} -group {Grouplist} " \
                           "{deg_para} -outdir {outdir}\n".format(
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
                degshell += "perl {GeneDiffExpBin}/PossionDis.pl  -list {CompareList}  {deg_para}  -outdir {outdir}\n".format(
                    GeneDiffExpBin=self.scriptbin,
                    CompareList=CompareList,
                    deg_para=self.parameter["PossionDis_Filter"],
                    outdir=self.outdir + 'PossionDis'
                )
        cmd.append(degshell)
        return cmd, output

    def makedefault(self,inputfq):

        gxp = geneexp()
        gxp.species=self.species
        gxp.fqLink=self.fqLink
        gxp.outdir=self.outdir.replace("GeneDiffExp_Allin","GeneExp")
        ExpDict = gxp.makedefault(inputfq)["output"][1]
        output=[]

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
        self.outdir = "RNAdenovo/GO_Hypergeometric/GO"
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/Enrichment"

        self.parameter = {
            "Annotation_dbClass": "pl"
        }

        self.program = "/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2016a/Enrichment/go2.pl"

    def makeCommand(self,inputfq):
        GeneDiffExpFilter = inputfq[0]
        Annotation=inputfq[1]
        prefix_name=os.path.basename(Annotation).split('.')[0]
        prefix_path = os.path.split(Annotation)[0]
        GO_prefix = prefix_path+"/"+prefix_name
        Gene2Tr = inputfq[2]

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)


        cmd=[]
        output=[]
        GODict={}
        go_shell=""

        go_tmpdir=self.outdir+"/tmp_file"
        os.makedirs(go_tmpdir, mode=0o755, exist_ok=True)

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            go_shell+="export PATH=/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/perl-V5/bin/:$PATH; " \
                     "awk '{{if($5>0) print $1\"\\t\"$5\"\\tup\";else print $1\"\\t\"$5\"\\tdown\"}}'" \
                     " {diff_list} >{tmpdir}/{keyname}.glist; " \
                     "perl {scriptbin}/drawGO.pl -list {tmpdir}/{keyname}.glist -goclass {goclass} " \
                     "-goprefix {prefix} -outprefix {outdir}/{keyname}; ".format(
                    diff_list=diff_list,
                    keyname=diff_id,
                    tmpdir=go_tmpdir,
                    prefix=GO_prefix,
                    scriptbin=self.scriptbin,
                    outdir=self.outdir,
                    goclass=dbclass.GOCLASS
            )
            GODict[diff_id]=[self.outdir+"/"+diff_id+"_C.txt",self.outdir+"/"+diff_id+"_F.txt",self.outdir+"/"+diff_id+"_P.txt"]
        go_shell+="perl {scriptbin}/go.pl -gldir {tmpdir} -sdir `dirname {prefix}` -species `basename {prefix}` -outdir {outdir};" \
                  "perl {scriptbin}/topGO.pl -gldir {tmpdir} -godir {outdir} -outdir {outdir} -prefix  {prefix} -list {gene2tr} -outdir {outdir}".format(
            tmpdir=go_tmpdir,
            prefix=GO_prefix,
            scriptbin=self.scriptbin,
            outdir=self.outdir,
            gene2tr=Gene2Tr
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

        annotation_o = annotation()
        annotation_o.outdir=self.outdir.replace("GO_Hypergeometric/GO","Annotation_Blast")
        annotation_o.species=self.species
        annotation_o.fqLink=self.fqLink
        GO_Annotation = annotation_o.makedefault(inputfq)["output"][-1]
        trinity_o = trinity_assemble()
        trinity_o.outdir = self.outdir.replace("GO_Hypergeometric/GO","Denovo_Trinity")
        trinity_o.species=self.species
        trinity_o.fqLink=self.fqLink
        Gene2Tr = trinity_o.makedefault(inputfq)["output"][-1]

        input=[]
        input.append(GeneDiffExpFilter)
        input.append(GO_Annotation)
        input.append(Gene2Tr)

        self.parameter["Annotation_dbClass"] = deg.species["RNAdenovo"][0]

        GODict={}
        output = []
        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            GODict[diff_id] = [self.outdir + "/" + diff_id + "_C.txt", self.outdir + "/" + diff_id + "_F.txt",self.outdir + "/" + diff_id + "_P.txt"]
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
        self.outdir = "RNAdenovo/Pathway_Hypergeometric/Pathway"
        self.parameter = {
            "Annotation_dbClass":"pl"
        }
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/Enrichment"
        self.program = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/Enrichment/pathfind.pl"

    def makeCommand(self,inputfq):
        GeneDiffExpFilter=inputfq[0]
        bg_ko = inputfq[1]

        kegg_tmpdir = self.outdir + "/tmp_file"
        os.makedirs(kegg_tmpdir,mode=0o755, exist_ok=True)

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
            os.system("grep -v \'Human Diseases\' "+ dbclass.KEGG_MAP_TITLE +" >"+kegg_tmpdir+"/map_title.tab")
            dbclass.KEGG_MAP_TITLE =kegg_tmpdir+"/map_title.tab"
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)

        cmd=[]
        output=[]
        KEGGDict={}
        kegg_shell=""
        keggdraw_shell=""

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            os.makedirs(self.outdir + '/'+diff_id+'_map',mode=0o755, exist_ok=True)
            kegg_shell+="awk '{{print $1\"\\t\"$5}}' {diff_list} >{tmpdir}/{keyname}.glist; " \
                       "/usr/bin/perl {scriptbin}/getKO.pl -glist {tmpdir}/{keyname}.glist -bg {bg_ko} -outdir {tmpdir}; " \
                       "/usr/bin/perl {scriptbin}/pathfind.pl -kegg {kegg_fa} -komap {ko_map} -maptitle {map_title} -fg {tmpdir}/{keyname}.ko " \
                       "-bg {bg_ko} -output {outdir}/{keyname}.path; " \
                       "awk '{{if($5>0) print $1\"\\t\"$5\"\\tup\";else print $1\"\\t\"$5\"\\tdown\"}} ' {diff_list}" \
                       " > {tmpdir}/{keyname}.glist.temp; " \
                       "perl {scriptbin}/keggMap.pl -ko {tmpdir}/{keyname}.ko -diff {tmpdir}/{keyname}.glist -komap " \
                       "{ko_map} -mapdir {mapdir} -outdir {outdir}/{keyname}_map; " \
                       "perl {scriptbin}/drawKEGG.pl -path {outdir}/{keyname}.path -outprefix {outdir}/{keyname} " \
                       "-idCol 6 -level1Col 7 -level2Col 8 -geneCol 9 -list {tmpdir}/{keyname}.glist.temp\n" \
                       "".format(
                diff_list=diff_list,
                tmpdir=kegg_tmpdir,
                keyname=diff_id,
                bg_ko=bg_ko,
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
        cmd.append(kegg_shell)
        cmd.append(keggdraw_shell)
        output.append(KEGGDict)
        return cmd,output

    def makedefault(self,inputfq):
        deg=genediffexp()
        deg.fqLink=self.fqLink
        deg.species=self.species
        deg.outdir=self.outdir.replace("RNAdenovo/Pathway_Hypergeometric/Pathway","RNAdenovo/GeneDiffExp_Allin")
        GeneDiffExpFilter = deg.makedefault(inputfq)["output"]

        annotation_o = annotation()
        annotation_o.species=self.species
        annotation_o.fqLink=self.fqLink
        annotation_o.outdir=self.outdir.replace("RNAdenovo/Pathway_Hypergeometric/Pathway","RNAdenovo/Annotation_Blast")
        Annotation = annotation_o.makedefault(inputfq)["output"]
        bg_ko =Annotation[-2]

        input=[]
        input.append(GeneDiffExpFilter)
        input.append(bg_ko)

        KEGGDict={}
        output=[]
        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            KEGGDict[diff_id] = [self.outdir + "/" + diff_id + ".path"]
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
        super(ppi, self).__init__()
        self.outdir = "RNAdenovo/PPI_Interaction"
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/PPI"

        self.parameter = {
            "Annotation_dbClass": "pl"
        }
        self.convert=self.getsoftware().CONVERT
        self.diamond=self.getsoftware().DIAMOND
        self.ppi="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/PPI/RNA_DENOVO_PPI.py"
        self.R = self.getsoftware().RSCRIPT325
        self.program=[self.convert,self.diamond,self.R,self.ppi]


    def makeCommand(self, inputfq):

        GeneDiffExpFilter = inputfq[0]
        annotation_xls = inputfq[2]
        Unigene = inputfq[1]

        if self.parameter["Annotation_dbClass"] == "pl":
            dbclass = DataBasePath(dbclass="pl")
        elif self.parameter["Annotation_dbClass"] == "an":
            dbclass = DataBasePath(dbclass="an")
        elif self.parameter["Annotation_dbClass"] == "fg":
            dbclass = DataBasePath(dbclass="fg")
        else:
            print("Your dbClass has some error!")
            sys.exit(1)

        cmd=[]
        output=[]
        ppi_shell="grep  -vw 'NA' {annotation} >{outdir}/annotation.xls; {diamond}   blastx --evalue  1e-5  --outfmt  6  -d " \
                  "{db} -q {Unigene} -o {outdir}/ppi.outfmt6  --seg no  --threads 5 " \
                  " --max-target-seqs 1 --more-sensitive -b 0.5 --salltitles;".format(
            annotation=annotation_xls,
            diamond=self.diamond,
            db=dbclass.STRING,
            Unigene=Unigene,
            outdir=self.outdir
        )

        os.makedirs(self.outdir, mode=0o755, exist_ok=True)

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            ppi_shell+="/usr/bin/python {scriptbin}/RNA_DENOVO_PPI.py --m6 {outdir}/ppi.outfmt6 --db {db_relation} --protein_list" \
                       " {annotation_xls} --diff_list {diff_list} --result {outdir}/{diff_id}.network.relation.txt;" \
                       "awk 'NR>1{{print $0}}' {outdir}/{diff_id}.network.relation.txt|sort -k5nr |head -100 >" \
                       "{outdir}/{diff_id}.network.relation.temp;/usr/bin/python {scriptbin}/fd_relaion.py " \
                       "{outdir}/{diff_id}.network.relation.temp {diff_list} >{outdir}/{diff_id}.network.relation.fd.temp; " \
                       "{R} {scriptbin}/PPI.R {outdir}/{diff_id}.network.relation.temp" \
                       " {outdir}/{diff_id}.network.relation.fd.temp {outdir}/{diff_id}.network.pdf;" \
                       "{convert} -density 300 {outdir}/{diff_id}.network.pdf {outdir}/{diff_id}.network.png;".format(
                scriptbin=self.scriptbin,
                annotation_xls=annotation_xls,
                outdir=self.outdir,
                diff_id=diff_id,
                R=self.R,
                convert=self.convert,
                db_relation=dbclass.STRING_RELATION,
                diff_list=diff_list
            )
            output.append(self.outdir + "/"+diff_id+".network.pdf")
            output.append(self.outdir + "/" + diff_id + ".network.png")
        cmd.append(ppi_shell)

        return cmd,output

    def makedefault(self,inputfq):
        deg = genediffexp()
        deg.species=self.species
        deg.fqLink=self.fqLink
        deg.outdir=self.outdir.replace("PPI_Interaction","GeneDiffExp_Allin")
        GeneDiffExpFilter = deg.makedefault(inputfq)["output"]

        annotation_o = annotation()
        annotation_o.species=self.species
        annotation_o.fqLink=self.fqLink
        annotation_o.outdir=self.outdir.replace("PPI_Interaction","Annotation_Blast")
        annotation_xls = annotation_o.makedefault(inputfq)["output"][-4]

        trinity_o = trinity_assemble()
        trinity_o.outdir=self.outdir.replace("PPI_Interaction","Denovo_Trinity/")
        trinity_o.fqLink=self.fqLink
        trinity_o.species=self.species
        Unigene = trinity_o.makedefault(inputfq)["output"][0]

        input=[]
        output=[]
        input.append(GeneDiffExpFilter)
        input.append(Unigene)
        input.append(annotation_xls)

        for diff_list in GeneDiffExpFilter:
            diff_id = ".".join(os.path.basename(diff_list).split('.')[0:2])
            output.append(self.outdir + "/"+diff_id+".network.pdf")
            output.append(self.outdir + "/" + diff_id + ".network.png")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "2G,5CPU",
            'output': output
        }
        return default

class snp(common):
    def __init__(self):
        super(snp, self).__init__()
        self.outdir = "RNAdenovo/SNP_GATK"
        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/SNP"

        self.parameter = {
            "Annotation_dbClass": "pl"
        }

        self.gatk=self.getsoftware().GATK
        self.java18=self.getsoftware().JAVA18
        self.java=self.getsoftware().JAVA
        self.hisat2=self.getsoftware().HISAT2+"hisat2"
        self.hisat_build=self.getsoftware().HISAT2+"hisat2-build"
        self.samtools=self.getsoftware().SAMTOOLS
        self.picard=self.getsoftware().PICARD
        self.diamond=self.getsoftware().DIAMOND
        self.R=self.getsoftware().RSCRIPT325
        self.program = [self.gatk,self.java18,self.java,self.hisat2,self.hisat_build,self.samtools,self.picard,self.diamond,self.R]

    def makeCommand(self,inputfq):
        Unigene=inputfq[0]
        CleanDataDict=inputfq[1]


        os.makedirs(self.outdir, mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/index-build", mode=0o755, exist_ok=True)

        index_shell=""
        cmd=[]
        output=[]

        index_shell+="cd {outdir}/index-build;perl {scriptbin}/sort_fa.pl {Unigene} refMrna.fa;".format(
            outdir=self.outdir,
            scriptbin=self.scriptbin,
            Unigene=Unigene
        )

        index_shell+="{hisat_build} refMrna.fa refMrna.fa;{samtools} faidx refMrna.fa;{picard}/CreateSequenceDictionary.jar R=refMrna.fa O=refMrna.dict" \
                     "".format(
            hisat_build=self.hisat_build,
            samtools=self.samtools,
            picard=self.picard
        )

        cmd.append(index_shell)

        gatk_shell=""
        depth_shell=""
        snp_file=""
        snp_basic_stat_shell = "{java18}  -Xmx2G -Djava.io.tmpdir={outdir}/java_tmp_snp  -jar {gatk} -T CombineVariants -R {refMrna} ".format(
            java18=self.java18,
            outdir=self.outdir,
            gatk=self.gatk,
            refMrna=self.outdir+"/index-build/refMrna.fa"
        )

        for SampleID, Cleandata in CleanDataDict.items():
            sample_dir=self.outdir+"/"+SampleID
            os.makedirs(sample_dir, mode=0o755, exist_ok=True)
            os.makedirs(sample_dir+"/java_tmp", mode=0o755, exist_ok=True)
            cleanFqA = Cleandata["clean_fq1"]
            cleanFqB = Cleandata["clean_fq2"]
            gatk_shell += "cd {sample_dir};{hisat2} --phred64 --sensitive --no-discordant --no-mixed -I 1 -X 1000 -x " \
                          "{refMrna} -1 {fq1} -2 {fq2} -S {sampleid}.sam;" \
                          "{java} -Djava.io.tmpdir={sample_dir}/java_tmp -jar {picard}/AddOrReplaceReadGroups.jar " \
                          "I={sampleid}.sam O={sampleid}.addRG.bam RGID={sampleid} RGLB={sampleid}_library RGPL=illumina RGPU=machine RGSM={sampleid} VALIDATION_STRINGENCY=SILENT;" \
                          "{java} -Djava.io.tmpdir={sample_dir}/java_tmp -jar {picard}/ReorderSam.jar I={sampleid}.addRG.bam O={sampleid}.addRG.Reorder.bam R={refMrna} VALIDATION_STRINGENCY=SILENT;" \
                          "{java} -Djava.io.tmpdir={sample_dir}/java_tmp -jar {picard}/SortSam.jar I={sampleid}.addRG.Reorder.bam O={sampleid}.addRG.Reorder.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT;" \
                          "{java} -Djava.io.tmpdir={sample_dir}/java_tmp -jar {picard}/MarkDuplicates.jar REMOVE_DUPLICATES=false I={sampleid}.addRG.Reorder.Sort.bam O={sampleid}.addRG.Reorder.Sort.markDup.bam METRICS_FILE={sampleid}.addRG.Reorder.Sort.markDup.metrics VALIDATION_STRINGENCY=SILENT;" \
                          "{samtools} index {sampleid}.addRG.Reorder.Sort.markDup.bam;" \
                          "{java18} -Djava.io.tmpdir={sample_dir}/java_tmp  -jar {gatk} -T SplitNCigarReads -R {refMrna}  -I {sampleid}.addRG.Reorder.Sort.markDup.bam -o {sampleid}.addRG.Reorder.Sort.markDup.splitN.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS;" \
                          "{samtools} index Control_1.addRG.Reorder.Sort.markDup.splitN.bam;" \
                          "{java18} -Xmx20G  -Djava.io.tmpdir={sample_dir}/java_tmp -jar {gatk} -T HaplotypeCaller -R {refMrna} -I {sampleid}.addRG.Reorder.Sort.markDup.splitN.bam -allowPotentiallyMisencodedQuals -stand_call_conf 20.0 -stand_emit_conf 20.0 -o {sampleid}.gatk.vcf;" \
                          "{java18} -Xmx20G  -Djava.io.tmpdir={sample_dir}/java_tmp -jar {gatk} -T SelectVariants -R {refMrna} -selectType SNP -V {sampleid}.gatk.vcf -o {sampleid}.gatk.select_snp.vcf;" \
                          "{java18} -Djava.io.tmpdir={sample_dir}/java_tmp -jar {gatk} -T VariantFiltration -R {refMrna} -V {sampleid}.gatk.select_snp.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o {sampleid}.raw.snp.vcf;" \
                          "awk '(/^#/ || $7 == \"PASS\")' {sampleid}.raw.snp.vcf >{sampleid}.snp.vcf\n".format(
                sample_dir=sample_dir,
                samtools=self.samtools,
                hisat2=self.hisat2,
                refMrna=self.outdir+"/index-build"+"/refMrna.fa",
                fq1=cleanFqA,fq2=cleanFqB,
                sampleid=SampleID,
                java=self.java,
                picard=self.picard,
                gatk=self.gatk,
                java18=self.java18
            )

            depth_shell +="{samtools} depth {sample_dir}/{sampleid}.addRG.Reorder.Sort.bam >{sample_dir}/{sampleid}.depth\n".format(
                samtools=self.samtools,
                sample_dir=sample_dir,
                sampleid=SampleID
            )

            snp_basic_stat_shell+="--variant: {sampleid} {sample_dir}/{sampleid}.snp.vcf ".format(
                sampleid=SampleID,
                sample_dir=sample_dir
            )

            snp_file+="{sample_dir}/{sampleid}.snp.vcf,".format(
                sample_dir=sample_dir,
                sampleid=SampleID
            )

            output.append(sample_dir+"/"+SampleID+".snp.vcf")
            output.append(sample_dir+"/"+SampleID+".depth")


        cmd.append(gatk_shell)
        cmd.append(depth_shell)

        snp_basic_stat_shell+=""
        os.makedirs(self.outdir+"/java_tmp_snp", mode=0o755, exist_ok=True)

        depth_list=[]
        for i in range(len(output)):
            if i%2 == 0:
                pass
            else:
                depth_list.append(output[i])

        snp_basic_stat_shell+=" -o {outdir}/All.snp.combine.vcf -genotypeMergeOptions UNIQUIFY;perl {scriptbin}/snp_statistics.pl -snp {snp_file} -outdir {outdir}".format(
            outdir=self.outdir,
            scriptbin=self.scriptbin,
            snp_file=snp_file
        )
        depth_file=",".join(depth_list)
        snp_population_stat_shell ="{perl} {scriptbin}/snp_population.pl -snp {snp_file} -depth {depth_file} -outdir {outdir}".format(
            scriptbin=self.scriptbin,
            perl="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/perl-V5/bin/perl",
            snp_file=snp_file,
            depth_file=depth_file,
            outdir=self.outdir
        )

        cmd.append(snp_basic_stat_shell)
        cmd.append(snp_population_stat_shell)

        return cmd,output

    def makedefault(self,inputfq):

        trinity_o = trinity_assemble()
        trinity_o.species = self.species
        trinity_o.fqLink = self.fqLink
        trinity_o.outdir = self.outdir.replace("SNP_GATK", "Denovo_Trinity")
        Unigene = trinity_o.makedefault(inputfq)["output"][0]

        filter_o = filter()
        filter_o.fqLink = self.fqLink
        filter_o.species = self.species
        filter_o.outdir = self.outdir.replace("SNP_GATK", "Filter_SOAPnuke")
        CleanDataDict = filter_o.makedefault(inputfq)["output"][0]

        input=[]
        output=[]
        input.append(Unigene)
        input.append(CleanDataDict)

        for SampleID, Cleandata in CleanDataDict.items():
            sample_dir=self.outdir+"/"+SampleID
            output.append(sample_dir+"/"+SampleID+".snp.vcf")
            output.append(sample_dir+"/"+SampleID+".depth")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "18G,1CPU",
            'output': output
        }
        return default

class tf(common):
    def __init__(self):
        super(tf, self).__init__()
        self.outdir = "RNAdenovo/TFpredict"
        self.diamond=self.getsoftware().DIAMOND
        self.getorf=self.getsoftware().GETORF
        self.rscript = self.getsoftware().RSCRIPT
        self.convert = self.getsoftware().CONVERT
        self.program=[self.diamond,self.getorf,self.rscript,self.convert]

        self.parameter = {
            "Annotation_dbClass": "pl"
        }

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/TF"

    def makeCommand(self, inputfq):

        all_unigene_fa = inputfq[0]
        all_gene_exp=inputfq[-1]
        os.makedirs(self.outdir, mode=0o775, exist_ok=True)
        database = DataBasePath(dbclass=self.parameter["Annotation_dbClass"])

        tf_sh = ""
        cmd = []
        output = []

        if self.parameter["Annotation_dbClass"] == "an":
            db = database.TF
            tf_sh += "{diamond} blastp --evalue 1e-5 --threads 5 --outfmt 6 -d {db} -q {unigene} -o {outdir}/denovo_tf_result.outfmt6 --seg no --max-target-seqs 1 --more-sensitive -b 0.5 --salltitles;" \
                     "awk \'!uniq[$1]++\' {outdir}/denovo_tf_result.outfmt6 >{outdir}/denovo_tf_best_result.outfmt6;" \
                     "echo -e \"gene_id\\tEnsembl ID\\tFamily\\twebsite\" >{outdir}/Animal_TF_result.xls;" \
                     "for i in `ls {scriptbin}/Animal-TF-tab/*xls`;do awk -F\'\\t\' \'{{print $2\"\\t\"$(NF-1)\"\\t\"$NF}}\' $i; " \
                     "done|sort|uniq|awk -F\'\\t\' \'NR==FNR{{a[$1]=$0;next}}{{split($2,array,/:/);if (array[1] in a) print $1\"\\t\"a[array[1]];else if(array[2] in a) print $1\"\\t\"a[array[2]]}}\' - {outdir}/denovo_tf_best_result.outfmt6 >{outdir}/Animal_TF_result.xls;" \
                     "/usr/bin/python {scriptbin}/animal_tf.py --tf {outdir}/Animal_TF_result.xls |sed \"1i TF_family\\tNumber of Genes\\tGenes\" >{outdir}/All-Unigene.TF2Gene.xls; " \
                     "{Rscript} {scriptbin}/draw_bar.R {outdir}/All-Unigene.TF2Gene.xls {outdir}/All-Unigene.TF_family.pdf;" \
                     "{convert} -density 300 -resize 30% {outdir}/All-Unigene.TF_family.pdf {outdir}/All-Unigene.TF_family.png;" \
                     "perl {scriptbin}/TF_heatmap.R.pl {All_GeneExpression_FPKM} {outdir}/Animal_TF_result.xls {outdir}".format(
                Rscript=self.rscript,
                convert=self.convert,
                diamond=self.diamond,
                unigene=all_unigene_fa,
                outdir=self.outdir,
                All_GeneExpression_FPKM=all_gene_exp,
                scriptbin=self.scriptbin,
                db=db
            )
            output.append(self.outdir + "/TF_heatmap.png")
            output.append(self.outdir + "/All-Unigene.TF2Gene.xls")
            output.append(self.outdir + "/All-Unigene.TF_family.png")
        elif self.parameter["Annotation_dbClass"] == "pl":
            tf_sh += "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{scriptbin}/../software;" \
                     "{getorf} -minsize 150 -sequence {unigene} -outseq {outdir}/All-Unigene.orf;" \
                     "perl {scriptbin}/select_orf.pl -seq {unigene} -input {outdir}/All-Unigene.orf -output {outdir}/All-Unigene.pep;" \
                     "perl {scriptbin}/TFCodingGene_Predict.pl {outdir}/All-Unigene.pep All-Unigene {outdir};" \
                     "perl {scriptbin}/TF_heatmap.R.pl {All_GeneExpression_FPKM} {outdir}/All-Unigene.TFCodingGene.xls {outdir};" \
                     "".format(
                Rscript=self.rscript,
                convert=self.convert,
                diamond=self.diamond,
                unigene=all_unigene_fa,
                outdir=self.outdir,
                All_GeneExpression_FPKM=all_gene_exp,
                scriptbin=self.scriptbin,
                getorf=self.getorf
            )
            output.append(self.outdir + "/TF_heatmap.png")
            output.append(self.outdir + "/All-Unigene.TF2Gene.xls")
            output.append(self.outdir + "/All-Unigene.TFCodingGene.xls")
            output.append(self.outdir + "/All-Unigene.TF_family.png")
        elif self.parameter["Annotation_dbClass"] == "fg":
            pass
        else:
            sys.exit(1)
        cmd.append(tf_sh)


        return cmd,output

    def makedefault(self, inputfq):
        trinity_o = trinity_assemble()
        trinity_o.species=self.species
        trinity_o.fqLink=self.fqLink
        trinity_o.outdir=self.outdir.replace("TFpredict","Denovo_Trinity")
        all_unigene_fa = trinity_o.makedefault(inputfq)["output"][0]

        gxp = geneexp()
        gxp.outdir=self.outdir.replace("TFpredict","GeneExp")
        gxp.fqLink=self.fqLink
        gxp.species=self.species
        all_trans_fpkm = gxp.makedefault(inputfq)["output"][-1]

        input=[]
        output=[]

        self.parameter["Annotation_dbClass"]=self.species["RNAdenovo"][0]

        input.append(all_unigene_fa)
        input.append(all_trans_fpkm)
        if self.parameter["Annotation_dbClass"] == "an":
            output.append(self.outdir + "/TF_heatmap.png")
            output.append(self.outdir + "/All-Unigene.TF2Gene.xls")
            output.append(self.outdir + "/All-Unigene.TF_family.png")
        elif self.parameter["Annotation_dbClass"] == "pl":
            output.append(self.outdir+"/TF_heatmap.png")
            output.append(self.outdir+"/All-Unigene.TF2Gene.xls")
            output.append(self.outdir+"/All-Unigene.TFCodingGene.xls")
            output.append(self.outdir+"/All-Unigene.TF_family.png")
        else:
            print ("Your db set is error!")
            os._exit(1)
        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "1G,1CPU",
            'output': output
        }

        return  default

class prg(common):
    def __init__(self):
        super(prg, self).__init__()
        self.outdir = "RNAdenovo/PRG_Blastx"
        self.diamond=self.getsoftware().DIAMOND
        self.program=[self.diamond]

        self.parameter = {
            "Annotation_dbClass": "pl"
        }

        self.scriptbin = "/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAdenovo/PRG/bin"

    def makeCommand(self, inputfq):

        all_unigene_fa = inputfq[0]

        os.makedirs(self.outdir, mode=0o775, exist_ok=True)

        database = DataBasePath(dbclass=self.parameter["Annotation_dbClass"])

        prg_sh = ""
        cmd = []
        output = []

        if self.parameter["Annotation_dbClass"] == "an":
            prg_sh +="echo \"#Your sample is animal, you don't need to run this step.\" >{outdir}/Unigene2PRG.result.xls".format(
                    outdir=self.outdir
            )
        elif self.parameter["Annotation_dbClass"] == "pl":
            prg_sh += "ln -sf {unigene} {outdir}/reference.fa;" \
                     "grep \'>\' {unigene} |sed \'s/>//g\'|awk \'{{print $1\"\\t\"$1}}\' > {outdir}/reference.gene2tr;" \
                     "{diamond} blastx --query {outdir}/reference.fa --db {prg_db} --evalue 1e-5 --threads 4 --outfmt 6 --max-target-seqs 1 --more-sensitive -b 0.5 --out {outdir}/degtr2PRG.out;" \
                     "perl {scriptbin}/get_PRG.pl -trseq {outdir}/reference.fa  -blresult {outdir}/degtr2PRG.out -gene2tr {outdir}/reference.gene2tr  -anno {prg_db_an} -cov 50 -iden 40 -output {outdir}/Unigene2PRG.result.xls;"\
                .format(
                diamond=self.diamond,
                unigene=all_unigene_fa,
                prg_db=database.PRG,
                prg_db_an = database.PRG_AN,
                outdir=self.outdir,
                scriptbin=self.scriptbin,
            )
        elif self.parameter["Annotation_dbClass"] == "fg":
            prg_sh += "echo \"#Your sample is fg, you don't need to run this step.\" >{outdir}/Unigene2PRG.result.xls".format(
                        outdir=self.outdir
            )
        else:
            sys.exit(1)
        cmd.append(prg_sh)
        output.append(self.outdir+"/Unigene2PRG.result.xls")
        return cmd,output

    def makedefault(self, inputfq):
        trinity_o = trinity_assemble()
        trinity_o.species=self.species
        trinity_o.fqLink=self.fqLink
        trinity_o.outdir=self.outdir.replace("PRG_Blastx","Denovo_Trinity")
        all_unigene_fa = trinity_o.makedefault(inputfq)["output"][0]

        input=[]
        output=[]

        self.parameter["Annotation_dbClass"]=self.species["RNAdenovo"][0]

        input.append(all_unigene_fa)
        output.append(self.outdir+"/Unigene2PRG.result.xls")

        default={
            'input': input,
            'parameter': self.parameter,
            'program': self.program,
            'resource': "1G,4CPU",
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
        os.makedirs(self.outdir+"/2.Assembly", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/3.Annotation", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir+"/4.Structure", mode=0o755, exist_ok=True)
        os.makedirs(self.outdir + "/5.Quantify", mode=0o755, exist_ok=True)
        cpshell=""
        cmd=[]
        output=[]
        cpshell +="mkdir -p {outdir}/4.Structure/CDSpredict/CDSpredict;" \
                  "mkdir -p {outdir}/4.Structure/PHI;" \
                  "mkdir -p {outdir}/4.Structure/PRG;" \
                  "mkdir -p {outdir}/4.Structure/SNP;" \
                  "mkdir -p {outdir}/4.Structure/SSR;" \
                  "mkdir -p {outdir}/4.Structure/TFpredict;" \
                  "mkdir -p {outdir}/5.Quantify/DifferentiallyExpressedGene/Functional_Enrichment/Pathway;" \
                  "mkdir -p {outdir}/5.Quantify/DifferentiallyExpressedGene/Functional_Enrichment/GO;" \
                  "mkdir -p {outdir}/5.Quantify/DifferentiallyExpressedGene/HierarchicalCluster;" \
                  "mkdir -p {outdir}/5.Quantify/DifferentiallyExpressedGene/PPI;" \
                  "mkdir -p {outdir}/5.Quantify/GeneExpression/Clustering_Mfuzz;" \
                  "mkdir -p {outdir}/5.Quantify/PCA;" \
                  "mkdir -p {outdir}/5.Quantify/Venny;" \
                  "cp {filteroutdir}/*/*.filter.stat.xls {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.RawReadsClass.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.base.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/*/*.qual.png {outdir}/1.CleanData/;" \
                  "cp {filteroutdir}/FilterSummary.xls {outdir}/1.CleanData/;" \
                  "cp {trinityoutdir}/All-Unigene.fa {outdir}/2.Assembly;" \
                  "cp {trinityoutdir}/Unigene.fa {outdir}/2.Assembly;" \
                  "cp {trinityoutdir}/All-Unigene.gene2mark {outdir}/2.Assembly;" \
                  "cp {annotationoutdir}/{{*xls,*pdf,*png,*path,*.htm,*.ko,*map}} {outdir}/3.Annotation/;" \
                  "cp {ssroutdir}/{{All-Unigene.misa.xls,All-Unigene.statistics.xls,SSR_statistics.png,SSR_statistics.xls}} {outdir}/4.Structure/SSR/;" \
                  "cp {cdsoutdir}/{{PredictSummary.xls,All-Unigene.fa.transdecoder.cds*,All-Unigene.fa.transdecoder.pep,All-Unigene.fa.transdecoder.gff3,All-Unigene.fa.transdecoder.bed,All-Unigene.cds.length.png,All-Unigene.cds.length.pdf}} {outdir}/4.Structure/CDSpredict/CDSpredict/;" \
                  "cp {geneexpoutdir}/*/*.gene.fpkm.xls {outdir}/5.Quantify/GeneExpression/;" \
                  "cp {geneexpoutdir}/*/*.Bowtie2Gene.MapReadsStat.xls {outdir}/5.Quantify/GeneExpression/;" \
                  "cp {geneexpoutdir}/{{All.GeneExpression.FPKM.xls,MappingSummary.xls}} {outdir}/5.Quantify/GeneExpression/;" \
                  "cp {geneexpoutdir}/GeneExpression/* {outdir}/5.Quantify/GeneExpression/;" \
                  "cp {genediffexpoutdir}/*/*GeneDiffExp*xls {genediffexpoutdir}/*/*.MA-plot.* {genediffexpoutdir}/*/*.Scatter-plot.* {genediffexpoutdir}/*/*.Volcano-plot.* {outdir}/5.Quantify/DifferentiallyExpressedGene/;" \
                  "cp {goenrichmentoutdir}/* {outdir}/5.Quantify/DifferentiallyExpressedGene/Functional_Enrichment/GO/;" \
                  "cp {pathwayenrichmentoutdir}/{{*.xls,*.htm,*.pdf,*.png,*.path,*map}} {outdir}/5.Quantify/DifferentiallyExpressedGene/Functional_Enrichment/Pathway/;" \
                  "cp {ppioutdir}/{{*.network.pdf,*.network.png,*.network.relation.txt}} {outdir}/5.Quantify/DifferentiallyExpressedGene/PPI/;" \
                  "cp {snpoutdir}/snp_population.xlsx {outdir}/4.Structure/SNP/;" \
                  "cp {snpoutdir}/{{SnpSummary*,All.snp.combine.vcf}} {outdir}/4.Structure/SNP/;" \
                  "cp {snpoutdir}/*/*.snp.vcf {outdir}/4.Structure/SNP/;" \
                  "cp {tfoutdir}/{{*pdf,*png,*xls}} {outdir}/4.Structure/TFpredict/;" \
                  "cp {phioutdir}/Unigene2PHI.result.xls {outdir}/4.Structure/PHI/;" \
                  "cp {prgoutdir}/Unigene2PRG.result.xls {outdir}/4.Structure/PRG/;".format(
            filteroutdir=outd+"/Filter_SOAPnuke",
            trinityoutdir=outd+"/Denovo_Trinity",
            annotationoutdir=outd+"/Annotation_Blast",
            ssroutdir=outd+"/SSR_MISA",
            cdsoutdir=outd+"/CDSpredict_TransDecoder",
            alignmentoutdir=outd+"/GenomeMapping_HISAT",
            geneexpoutdir=outd+"/GeneExp",
            genediffexpoutdir=outd+"/GeneDiffExp_Allin",
            goenrichmentoutdir=outd+"/GO_Hypergeometric/GO",
            pathwayenrichmentoutdir=outd+"/Pathway_Hypergeometric/Pathway",
            ppioutdir=outd+"/PPI_Interaction",
            snpoutdir=outd+"/SNP_GATK",
            tfoutdir=outd+"/TFpredict",
            phioutdir=outd+"/PHI_Blastx",
            prgoutdir=outd+"/PRG_Blastx",
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
        common.__init__(self)
        self.step = [["filter"], ["trinity_assemble"], ["geneexp","annotation","cds_predict","ssr","snp"], ["genediffexp"], ["goenrichment", "pathwayenrichment","ppi","tf"]]

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
                    outshell = open(step+".sh", mode='w')
                    astep = eval(step)
                    astepo = astep()
                    astepo.fqLink = self.fqLink
                    stepdict = json.dumps(astepo.makedefault(self.fqList))
                    out.write("\"%s\":%s,\n" % (step, stepdict))
                    outshell.close()
                    os.system("sed -i s/;/\n/" +step+".sh")
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

