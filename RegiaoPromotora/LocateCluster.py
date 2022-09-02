#!usr/bin/env python3

import argparse
import pyranges as pr
import pandas as pd
import os
import datetime
import numpy as np
import re

parser = argparse.ArgumentParser(description='co-expression dataframe')
parser.add_argument('csv_data', type=str, 
                                    help='gene and cluster in cvs format')

args = parser.parse_args()

csv_data=args.csv_data

genecluster=pd.read_csv(csv_data, sep=";")

print(genecluster)

dictgenecluster={}
for i,r in genecluster.iterrows():
    if r['Cluster'] not in dictgenecluster.keys():
        dictgenecluster[r['Cluster']]=[]
    dictgenecluster[r['Cluster']].append(r['GeneID'])

for cluster in dictgenecluster:
    fastaoutput = "Cluster_" + str(cluster) + "_coexpressao_cana_RianoSouza.fa"
    with open (fastaoutput, "w") as clusterfile:
        for ID in dictgenecluster[cluster]:
            clusterfile.write(ID+"\n")


