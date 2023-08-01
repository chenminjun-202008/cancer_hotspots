#! /usr/bin/env python3
# -*- coding=utf-8 -*-
#================================================================

# Project       : snv annotation
# Description   : get revel  and splice socre etc
# Usage         : python anno.py -h
# Author        : chen.minjun
# Time          : 2023-07-17

#================================================================
import re
import sys
import os
import argparse
import json
import time
import pandas as pd


bindir = os.path.abspath(os.path.dirname(__file__))
author='chen.minjun'
mail= 'chen.minjun@genecast.com.cn'
doc= 'the decription of program'

#============function=============================================

#读取maf，获得snv/indel的蛋白突变信息和物理坐标关系
def Read_maf(infile):
    dict_maf={}
    dict_trans={}
    with open(infile,"r") as inpt:
        for line in inpt:
            if line.startswith("#"):
                pass
            elif line.startswith("Hugo_Symbol"):
                header=line.strip().split("\t")
            else:
                items = dict(zip(header,line.strip().split("\t")))
                if "Amino_Acid_Change" not in items: continue
                Amino_Acid_Change = items["Amino_Acid_Change"]
                gene = items["Hugo_Symbol"]
                AA = gene + ":" + Amino_Acid_Change
                hgvs = "\t".join([items[i] for i in ["Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","HGVSc"]])
                if AA not in dict_maf: dict_maf[AA] = []
                if hgvs not in dict_maf[AA]: dict_maf[AA].append(hgvs)

                trans = items["Transcript_ID"]
                if gene not in dict_trans: dict_trans[gene]=[]
                if trans not in  dict_trans[gene]: dict_trans[gene].append(trans)

    return(dict_maf,dict_trans)

#=========================================================
def main():
        parser=argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog='author:\t{0}\nmail:\t{1}'.format(author,mail))
        parser.add_argument("-m",dest="maf",help="cancerhotspots.v2.maf",required=True,type=str)
        parser.add_argument("-s",dest="hotspots",help="hotspots_v2.xls",required=True,type=str)
        parser.add_argument("-o",dest="outfile",help="out.hotspots_v2.addPOS.txt",required=True,type=str)
        args = parser.parse_args()

        maf_sites,transcript = Read_maf(args.maf)

        with open(args.outfile+".transcript","w") as trs:
            for gen in transcript:
                trs.write(gen+"\t"+"|".join(transcript[gen])+"\n")

        with open (args.hotspots,"r")  as hot, open(args.outfile,"w") as out, open(args.outfile+".err","w") as err,open(args.outfile+".err_for_transvar","w") as err2:

            addtile = ["Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","HGVSc"]
            header = hot.readline().strip("\n").split("\t")
            out.write("\t".join(addtile+header) + "\n")
            err.write("\t".join(addtile+header) + "\n")

            for line in hot:
                content = line.strip("\n").split("\t")
                items = dict( zip(header,content) )
                Amino_Acid_Change = items["Amino_Acid_Change"] 
                if Amino_Acid_Change in maf_sites:
                    for addinfos in maf_sites[Amino_Acid_Change]: 
                        new_line = addinfos + "\t" + line 
                        out.write(new_line)
                else:
                    new_line =  "\t".join(["-"]*6) + "\t" + line
                    for lst in transcript[items["Hugo_Symbol"]]:
                      transvar_use = lst + "\t" + Amino_Acid_Change +"\n"
                      err2.write(transvar_use)
                    err.write(new_line)


        








if __name__ == "__main__":
        print('Start : ',time.ctime())
        main()
        print('End : ',time.ctime())
            

