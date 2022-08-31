#!usr/bin/env python3

import argparse
import pyranges as pr
import pandas as pd
import os
import datetime
import numpy as np
import sys

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
    

print("Arquivos conferidos em", datetime.datetime.now()) #Feedback 1 - arquivos conferidos
sys.stdout.flush()

gff = pr.read_gff3(gff, as_df=False) #lê o arquivo gff3

print("Arquivo gff importado dentro do pyranges em", datetime.datetime.now()) #Feedback 2 - Importa o arquivo gff no pyranges 
sys.stdout.flush()

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

print("Arquivo preparado para o subtract", datetime.datetime.now(), "- Próximo passo: Subtract") #Feedback 3 - Arquivo preparado para o subtract
sys.stdout.flush()

regiaopromotora = gffgenecomparison.subtract(gff, strandedness="same", nb_cpu=mcores)  

print("Subtract efetuado em", datetime.datetime.now(), "- Próximo passo: excluir overlaps") #Feedback 4 - Subtract efetuado, próximo passo - excluir overlaps
sys.stdout.flush()

IDanterior="1"
for idpromoter in regiaopromotora.ID: 
    if idpromoter == IDanterior:
        reppromoter = regiaopromotora[regiaopromotora.ID == idpromoter]
        regiaopromotora = regiaopromotora[regiaopromotora.ID != idpromoter]
        Pgeneid = idpromoter.replace("Promoter_", "")
        promotorcorreto = reppromoter.nearest(gff[gff.ID == Pgeneid], strandedness="same", how="downstream")
        promotorcorreto = promotorcorreto[promotorcorreto.Distance == 1]
        promotorcorretoconcat=pr.PyRanges(chromosomes=promotorcorreto.Chromosome, starts=promotorcorreto.Start, ends=promotorcorreto.End, strands=promotorcorreto.Strand) #Necessita-se excluir algumas colunas que são denominadas "*_b", pois há um erro que insere novas colunas em IDs com o mesmo cromossomo
        promotorcorretoconcat.ID=promotorcorreto.ID
        promotorcorretoconcat.Source=promotorcorreto.Source
        promotorcorretoconcat.Feature=promotorcorreto.Feature
        promotorcorretoconcat.Name=promotorcorreto.Name
        promotorcorretoconcat.Score=promotorcorreto.Score
        promotorcorretoconcat.Frame=promotorcorreto.Frame
        regiaopromotora = pr.concat([regiaopromotora, promotorcorretoconcat])
    IDanterior=idpromoter


print("Overlaps e repetições conferidos e corrigidos em", datetime.datetime.now(), "- Próximo passo: escrever arquivo fasta com as regiões promotoras") #Feedback 5 - Correção dos overlaps e repetições e indicação da próxima etapa
sys.stdout.flush()

if args.fixid:
        regiaopromotora.Chromosome = regiaopromotora.Chromosome.astype(str) + (f'.{versao}')

#Necessita-se zerar o Start caso haja valores negativos, pois o get_fasta não ajusta esses valores e não os reconhece
#End não precisa ser ajustado

print("  -Zerando os valores de start negativos",  datetime.datetime.now())
sys.stdout.flush()

zero=np.int32(0) #Zero com inteiro no mesmo formato do pyranges para evitar mensagem de erro

zerarpromotor = regiaopromotora[regiaopromotora.Start < 0]
regiaopromotora = regiaopromotora[regiaopromotora.Start >= 0]
zerarpromotor.Start = zero
regiaopromotora = pr.concat([regiaopromotora, zerarpromotor])

print("  -Início da conversão para o arquivo fasta", datetime.datetime.now())
sys.stdout.flush()

#Escrevendo arquivo fasta com as regiões promotoras
with open (outputfile,"w") as TesteFastaSeqPromotora:
    for ID in regiaopromotora.ID:
        seqID = pr.get_fasta(regiaopromotora[regiaopromotora.ID == ID], genoma)
        for line in seqID:
            TesteFastaSeqPromotora.write(f'>{ID}\n{line}\n')

print("Termino do script - arquivo fasta com a regiões promotoras finalizado em", datetime.datetime.now()) #Feedback 6 - Correção dos overlaps e repetições e indicação da próxima etapa
sys.stdout.flush()

