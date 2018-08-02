#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

""" 
@author:yueyao 
@file: DataBasePath.py 
@time: 2018/05/11 
"""

class DataBasePath(object):

    def __init__(self,dbclass="pl"):

        self.OBO = "/ifs4/BC_PUB/biosoft/db/Pub/go/RNA/20171220/gene_ontology.1_2.obo"
        self.GOCLASS = "/ifs4/BC_PUB/biosoft/db/Pub/go/RNA/20171220/go.class"
        self.GENE2GO = "/ifs4/BC_PUB/biosoft/db/Pub/go/RNA/20171220/gene2go"
        self.ACCESSION2GO = "/ifs4/BC_PUB/biosoft/db/Pub/go/RNA/20171220/accession2go"
        self.KEGG_MAP_DIR = "/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/map"
        self.PFAM = "/ifs4/BC_PUB/biosoft/db/Pub/Pfam/Pfam-A.hmm"
        self.STRING="/ifs4/BC_PUB/biosoft/db/Pub/PPI/STRING/protein.sequences.v10.fa"
        self.STRING_RELATION="/ifs4/BC_PUB/biosoft/db/Pub/PPI/STRING/uniq.protein.links.v10.txt"
        self.INTERPROT = "/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/software/interproscan5"
        self.COG_ALL = "/ifs4/BC_PUB/biosoft/db/Pub/cog/RNA/20090331/cog_clean.fa"
        self.COG_GENEID = "/ifs4/BC_PUB/biosoft/db/Pub/cog/RNA/20090331/cog_clean.fa.id"
        self.COG_WHOG = "/ifs4/BC_PUB/biosoft/db/Pub/cog/RNA/20090331/whog"
        self.COG_FUN = "/ifs4/BC_PUB/biosoft/db/Pub/cog/RNA/20090331/fun.txt"
        self.KOG_ALL = "/ifs4/BC_PUB/biosoft/db/Pub/kog/RNA/20090331/data/kog_clean.fa"
        self.KOG_GENEID = "/ifs4/BC_PUB/biosoft/db/Pub/kog/RNA/20090331/data/kog_clean.fa.id"
        self.KOG_WHOG = "/ifs4/BC_PUB/biosoft/db/Pub/kog/RNA/20090331/data/twog.merge"
        self.KOG_FUN = "/ifs4/BC_PUB/biosoft/db/Pub/kog/RNA/20090331/data/fun.txt"
        self.SWISSPROT_GENEID = "/ifs4/BC_PUB/biosoft/db/Pub/swissprot/RNA/release-2017_09/uniprot_sprot.id.annot.xls"
        self.PRG = "/ifs4/BC_PUB/biosoft/db/Pub/PRG/PRG_20160913/all_protein_all.fa"
        self.PRG_AN = "/ifs4/BC_PUB/biosoft/db/Pub/PRG/PRG_20160913/all_annotation_all.tsv"
        self.PHI="/ifs4/BC_PUB/biosoft/db/Pub/PHI/PHI_v4.1/PHI.pep"

        if dbclass == "pl":
            self.NT="/ifs4/BC_PUB/biosoft/db/Pub/nt/20170924/Plants.fa"
            self.NR="/ifs4/BC_PUB/biosoft/db/Pub/nr/RNA/20170924/Plants.fa"
            self.KEGG="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/plant.fa"
            self.KEGG_KOMAP="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/komap/plant_ko_map.tab"
            self.KEGG_GENEID="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/plant.id.annot.xls"
            self.KEGG_MAP_TITLE="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/map_title.tab"
            self.SWISSPROT="/ifs4/BC_PUB/biosoft/db/Pub/swissprot/RNA/release-2017_09/Plants.fa"

        elif dbclass == "an":
            self.TF="/ifs4/BC_PUB/biosoft/db/Pub/AnimalTFDB/Animal-TF.fa"
            self.NT="/ifs4/BC_PUB/biosoft/db/Pub/nt/20170924/animal.fa"
            self.NR="/ifs4/BC_PUB/biosoft/db/Pub/nr/RNA/20170924/animal.fa"
            self.KEGG="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/animal.fa"
            self.KEGG_KOMAP="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/komap/animal_ko_map.tab"
            self.KEGG_GENEID="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/animal.id.annot.xls"
            self.KEGG_MAP_TITLE="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/map_title.tab"
            self.SWISSPROT="/ifs4/BC_PUB/biosoft/db/Pub/swissprot/RNA/release-2017_09/Animal.fa"

        elif dbclass == "fg":
            self.NT="/ifs4/BC_PUB/biosoft/db/Pub/nt/20170924/fungi.fa"
            self.NR="/ifs4/BC_PUB/biosoft/db/Pub/nr/RNA/20170924/fungi.fa"
            self.KEGG="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/fungi.fa"
            self.KEGG_KOMAP="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/komap/fungi_ko_map.tab"
            self.KEGG_GENEID="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/fungi.id.annot.xls"
            self.KEGG_MAP_TITLE="/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/map_title.tab"
            self.SWISSPROT="/ifs4/BC_PUB/biosoft/db/Pub/swissprot/RNA/release-2017_09/Fungi.fa"

        else:
            print ("The dbclass only can be pl,an or fg")
            sys.exit(1)

    def get_config(self,species="GT"):
        if species == "GT":
            self.GENOME="/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/GT.genome.fa"
            self.CDS="/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/GT.gene-v2.cds"
            self.GO_PREFIX="/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/GT"
            self.Gene2Symbol=""
            self.Gene2Tr="/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/GT.gene2mark"
            self.HISAT_INDEX="/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/GT.genome.fa"
            self.mRNA_GTF="/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/GT.gtf"
            self.lncRNA_GTF=""
            self.Protein_Fasta="/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/GT.gene-v2.pep"
            self.NR_Desc=""
            self.KO="/ldfssz1/ST_BIGDATA/USER/yueyao/17.Train/database/canu.Gene.pep.ko"
            self.DNA_VCF=""
        elif species == "human":
            self.GENOME="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/GenomeGatkIndex/chrALL.sort.fa"
            self.CDS="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/GeneBowtie2Index/refMrna.fa"
            self.GO_PREFIX="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg79/hg19"
            self.Gene2Symbol="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/human.gene2symbol.txt"
            self.Gene2Tr="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg76/refMrna.fa.gene2mark"
            self.HISAT_INDEX="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/GenomeHisat2Index/chrALL"
            self.mRNA_GTF="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation/hg19_filter.psl.gtf"
            self.lncRNA_GTF=""
            self.Protein_Fasta="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/CPC_format/HUMAN.protein.fa"
            self.NR_Desc="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg79/hg19.nr.desc"
            self.KO="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg79/hg19.ko"
            self.DNA_VCF="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2015a/example/dna_vcf.list"
        else:
            print("Your species is not in our database!")
            sys.exit(1)