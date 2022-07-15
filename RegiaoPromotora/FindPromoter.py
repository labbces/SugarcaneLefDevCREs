#!usr/bin/env python3

import argparse
import pyranges as pr
import pandas as pd

parser= argparse.ArgumentParser(description='Insert the lenth of the Expected Promoter')
parser.add_argument('integers', metavar='N', type=int,
                            help='insert the lenth of the Pormoter')
Promoter = parser.parse_args()
print(Promoter.integers)

lenPromoter = Promoter.integers

gff = pr.read_gff3("/mnt/d/joaon/TCC/Programação_TCC/GFFs/TESTE_CTBE.gff3", as_df=True) #lê o arquivo gff3

gffgeneplus=gff[(gff.Feature == "gene") & (gff.Strand == "+")] #separa fitas positivas
gffgeneplus2=gffgeneplus.assign(PromoterStart=lambda x: x.Start-1-lenPromoter)
gffgeneplus2=gffgeneplus2.assign(PromoterEnd=lambda x: x.Start-1)

gffgenepluscomparison = gffgeneplus.copy()
gffgenepluscomparison.ID = "Promoter_" + gffgenepluscomparison.ID.astype(str)
PromoterStartplus = gffgeneplus2.loc[:,"PromoterStart"]
gffgenepluscomparison.loc[:,"Start"] = PromoterStartplus
PromoterEndplus = gffgeneplus2.loc[:,"PromoterEnd"]
gffgenepluscomparison.loc[:,"End"] = PromoterEndplus

PRgplus = pr.PyRanges(gffgeneplus)
PRgpluscomparison = pr.PyRanges(gffgenepluscomparison)

print(PRgplus, PRgpluscomparison)
regiaopromotora = PRgpluscomparison.subtract(PRgplus, strandedness="same")  
SeqPromotora = pr.get_fasta(regiaopromotora, "/mnt/d/joaon/TCC/Programação_TCC/Genomas/TESTE_GCA_CTBE.fna")
regiaopromotora.seq=SeqPromotora
print(regiaopromotora)
print(SeqPromotora)

with open ("TesteFastaSeqPromotoraPlus.txt","w") as TesteFastaSeqPromotora:
    for line in SeqPromotora:
        TesteFastaSeqPromotora.write(line+"\n"+"\n")



gffgeneminus=gff[(gff.Feature == "gene") & (gff.Strand == "-")] #separa fitas negativas
gffgeneminus2=gffgeneplus.assign(PromoterEnd=lambda x: x.End+1)
gffgeneminus2=gffgeneplus2.assign(PromoterStart=lambda x: x.End+1+lenPromoter)

