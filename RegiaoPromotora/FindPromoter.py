#!usr/bin/env python3

import argparse
import pyranges as pr
import pandas as pd
import os

parser= argparse.ArgumentParser(description='Insert the lenth of the Expected Promoter')
parser.add_argument('length', type=int, 
                            help='insert the length of the Pormoter')

parser.add_argument('genoma', type=str,
                            help='genoma sequence in fasta format')

parser.add_argument('gff', type=str, 
                            help='genoma anotation in gff3 format')

parser.add_argument('outputfile', type=str, 
                            help='output file')

parser.add_argument('--fixid', action='store_true', 
                            help='store true to adjust the gene id adding version to the promoter id')
parser.add_argument('--v', type=int, action="store",
                            help='chromosome version in fasta file, use with fixid. Ex: 1,2,3...')
parser.add_argument('--mcores', type=int, default=1,
                            help='number of computing cores the job requests')

args = parser.parse_args()

genoma = args.genoma
gff = args.gff
lenPromoter = args.length
outputfile = args.outputfile
versao = args.v
mcores = args.mcores


#conferir se o arquivo existe no disco

if not os.path.isfile(genoma) or not os.path.isfile(gff):
    if not os.path.isfile(genoma) and os.path.isfile(gff):
        print('''Genoma file does not exist
Please check the file name/directory again''')
    if not os.path.isfile(gff) and os.path.isfile(genoma):
        print('''GFF file does not exist
Please check the file name/directory again''')
    if not os.path.isfile(genoma) and not os.path.isfile(gff):
        print('''Genoma file and GFF file do not exist
Please check the file name/directory again''')
    
    quit()
    

print("Arquivos conferidos") #Feedback 1 - arquivos conferidos


gff = pr.read_gff3(gff, as_df=False) #lê o arquivo gff3

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

regiaopromotora = gffgenecomparison.subtract(gff, strandedness="same", nb_cpu=mcores)  

#ler o arquivo da regiao promotora e vê se o ID correspondente no arquivo gffgeneplus "menos Promoter_" está a 1 index de distância, se não deleta ele do arquivo regiao promotora
# agora que achei o gene que repete, da para fazer um subset bseado no ID e armazenar em uma variável e apagar os que repetem do regiaopromotora. COmpara o subset com o gffgeneplus e ver qual é o nearest upstream (ou downstream para comparar e concatenar novamente no regiaopromotora.

IDanterior="1"
for idpromoter in regiaopromotora.ID: 
    if idpromoter == IDanterior:
        reppromoter = regiaopromotora[regiaopromotora.ID == idpromoter]
        regiaopromotora = regiaopromotora[regiaopromotora.ID != idpromoter]
        Pgeneid = idpromoter.replace("Promoter_", "")
        promotorcorreto = reppromoter.nearest(gff[gff.ID == Pgeneid], strandedness="same", how="downstream", nb_cpu=mcores)
        promotorcorreto = promotorcorreto[promotorcorreto.Distance == 1]
        regiaopromotora = pr.concat([regiaopromotora, promotorcorreto])
    IDanterior=idpromoter


print("Overlaps e repetições conferidos e corrigidos") #Feedback 2 - Correção dos overlaps e repetições

if args.fixid:
        regiaopromotora.Chromosome = regiaopromotora.Chromosome.astype(str) + (f'.{versao}')

#Escrevendo arquivo fasta com as regiões promotoras
with open (outputfile,"w") as TesteFastaSeqPromotora:
    for ID in regiaopromotora.ID:
        seqID = pr.get_fasta(regiaopromotora[regiaopromotora.ID == ID], genoma)
        for line in seqID:
            TesteFastaSeqPromotora.write(f'>{ID}\n{line}\n')
        

