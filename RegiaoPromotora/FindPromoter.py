#!usr/bin/env python3

import argparse
import pyranges as pr
import pandas as pd


parser= argparse.ArgumentParser(description='Insert the lenth of the Expected Promoter')
parser.add_argument('integers', metavar='N', type=int,
                            help='insert the lenth of the Pormoter')
Promoter = parser.parse_args()

lenPromoter = Promoter.integers

gff = pr.read_gff3(r"/mnt/d/joaon/TCC/Programação_TCC/GFFs/TESTE_CTBE.gff3", as_df=False) #lê o arquivo gff3

gffgeneplus=gff[(gff.Feature == "gene") & (gff.Strand == "+")] #separa fitas positivas
gffgeneplus2=gffgeneplus.assign("PromoterStart", lambda x: x.Start-lenPromoter)
gffgeneplus2=gffgeneplus2.assign("PromoterEnd", lambda x: x.Start)

gffgeneminus=gff[(gff.Feature == "gene") & (gff.Strand == "-")] #separa fitas negativas
gffgeneminus2=gffgeneminus.assign("PromoterStart", lambda x: x.End)
gffgeneminus2=gffgeneminus2.assign("PromoterEnd", lambda x: x.End+lenPromoter)

gffgenecomparison=pr.concat([gffgeneplus2, gffgeneminus2])

gffgenecomparison.ID = "Promoter_" + gffgenecomparison.ID.astype(str)
gffgenecomparison.Start=gffgenecomparison.PromoterStart
gffgenecomparison.End=gffgenecomparison.PromoterEnd

gffgenecomparison.Feature = "Promoter"

regiaopromotora = gffgenecomparison.subtract(gff, strandedness="same")  

#ler o arquivo da regiao promotora e vê se o ID correspondente no arquivo gffgeneplus "menos Promoter_" está a 1 index de distância, se não deleta ele do arquivo regiao promotora
# agora que achei o gene que repete, da para fazer um subset bseado no ID e armazenar em uma variável e apagar os que repetem do regiaopromotora. COmpara o subset com o gffgeneplus e ver qual é o nearest upstream (ou downstream para comparar e concatenar novamente no regiaopromotora.

IDanterior="1"
for idpromoter in regiaopromotora.ID: 
    if idpromoter == IDanterior:
        reppromoter = regiaopromotora[regiaopromotora.ID == idpromoter]
        regiaopromotora = regiaopromotora[regiaopromotora.ID != idpromoter]
        Pgeneid = idpromoter.replace("Promoter_", "")
        promotorcorreto = reppromoter.nearest(gff[gff.ID == Pgeneid], strandedness="same", how="downstream")
        promotorcorreto = promotorcorreto[promotorcorreto.Distance == 1]
        regiaopromotora = pr.concat([regiaopromotora, promotorcorreto])
    IDanterior=idpromoter

#Escrevendo arquivo fasta com as regiões promotoras
with open ("TesteFastaSeqPromotora.txt","w") as TesteFastaSeqPromotora:
    for ID in regiaopromotora.ID:
        seqID = pr.get_fasta(regiaopromotora[regiaopromotora.ID == ID], "/mnt/d/joaon/TCC/Programação_TCC/Genomas/TESTE_GCA_CTBE.fna")
        for line in seqID:
            TesteFastaSeqPromotora.write(f'>{ID}\n{line}\n')
        

