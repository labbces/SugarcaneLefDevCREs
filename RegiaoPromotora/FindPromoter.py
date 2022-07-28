#!usr/bin/env python3

import argparse
import pyranges as pr
import pandas as pd


parser= argparse.ArgumentParser(description='Insert the lenth of the Expected Promoter')
parser.add_argument('integers', metavar='N', type=int,
                            help='insert the lenth of the Pormoter')
Promoter = parser.parse_args()
#print(Promoter.integers)

lenPromoter = Promoter.integers



gff = pr.read_gff3(r"/mnt/d/joaon/TCC/Programação_TCC/GFFs/TESTE_CTBE.gff3", as_df=True) #lê o arquivo gff3

gffgeneplus=gff[(gff.Feature == "gene") & (gff.Strand == "+")] #separa fitas positivas
gffgeneplus2=gffgeneplus.assign(PromoterStart=lambda x: x.Start-1-lenPromoter)
gffgeneplus2=gffgeneplus2.assign(PromoterEnd=lambda x: x.Start)

gffgenepluscomparison = gffgeneplus.copy()
gffgenepluscomparison.ID = "Promoter_" + gffgenepluscomparison.ID.astype(str)
PromoterStartplus = gffgeneplus2.loc[:,"PromoterStart"]
gffgenepluscomparison.loc[:,"Start"] = PromoterStartplus
PromoterEndplus = gffgeneplus2.loc[:,"PromoterEnd"]
gffgenepluscomparison.loc[:,"End"] = PromoterEndplus

PRgplus = pr.PyRanges(gffgeneplus)
PRgpluscomparison = pr.PyRanges(gffgenepluscomparison)
PRgpluscomparison.Feature = "Promoter"

#1print(PRgplus, PRgpluscomparison)
regiaopromotora = PRgpluscomparison.subtract(PRgplus, strandedness="same")  
print(regiaopromotora)
#ler o arquivo da regiao promotora e vê se o ID correspondente no arquivo gffgeneplus "menos Promoter_" está a 1 index de distância, se não deleta ele do arquivo regiao promotora
# agora que achei o gene que repete, da para fazer um subset bseado no ID e armazenar em uma variável e apagar os que repetem do regiaopromotora. COmpara o subset com o gffgeneplus e ver qual é o nearest upstream (ou downstream para comparar e concatenar novamente no regiaopromotora.

IDanterior="1"
for idpromoter in regiaopromotora.ID: 
    if idpromoter == IDanterior:
        reppromoter = regiaopromotora[regiaopromotora.ID == idpromoter]
        regiaopromotora = regiaopromotora[regiaopromotora.ID != idpromoter]
        Pgeneid = idpromoter.replace("Promoter_", "")
        promotorcorreto = reppromoter.nearest(PRgplus[PRgplus.ID == Pgeneid], strandedness="same", how="downstream")
        promotorcorreto = promotorcorreto[promotorcorreto.Distance == 1]
        regiaopromotora = pr.concat([regiaopromotora, promotorcorreto])
    IDanterior=idpromoter

#SeqPromotora = pr.get_fasta(regiaopromotora, "/mnt/d/joaon/TCC/Programação_TCC/Genomas/TESTE_GCA_CTBE.fna")
#regiaopromotora.seq=SeqPromotora

with open ("TesteFastaSeqPromotoraPlus.txt","w") as TesteFastaSeqPromotora:
    for ID in regiaopromotora.ID:
        seqID = pr.get_fasta(regiaopromotora[regiaopromotora.ID == ID], "/mnt/d/joaon/TCC/Programação_TCC/Genomas/TESTE_GCA_CTBE.fna")
        for line in seqID:
            TesteFastaSeqPromotora.write(f'>{ID}\n{line}\n')
        

########## para a fita negativa
gffgeneminus=gff[(gff.Feature == "gene") & (gff.Strand == "-")] #separa fitas negativas
gffgeneminus2=gffgeneplus.assign(PromoterEnd=lambda x: x.End+1)
gffgeneminus2=gffgeneplus2.assign(PromoterStart=lambda x: x.End+1+lenPromoter)

