#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:30:21 2020

@author: diegogotex
"""


# -*- coding: utf-8 -*-

#carregando as bibliotecas
from Bio import SeqIO



#carregando o arquivo do genbank
#o ACC MK518395 tem um erro e a informação Assembly-Data tem que ser removida 
gbFile = "/Users/diegogotex/Desktop/MANEL/seqs/sequence_pept.gp"
#02/03/2020
#loop de teste para saber se deu certo carregar o arquivo
count = 0
for seq_record in SeqIO.parse(gbFile, "genbank"):
    #print(seq_record.id)
    count +=1
    
print(count)

#criando um arquivo de escrita
file = open("/Users/diegogotex/Desktop/MANEL/seqs/sequence.tab", "w+")

#loop sério 
for seq_record in SeqIO.parse(gbFile, "genbank"): 
    source = seq_record.features[0] #puxando as features de cada ACCESSION
    organism = str(source.qualifiers.get('organism')) #pega o valor com informação do organismo
    organism = organism.replace(" ", "_") #removendo os espaços
    organism = str(organism)[2:-2]#removendo a aspa e o colchete
    db_xref = str(source.qualifiers.get('db_xref')) #puxando o taxonID
    db_xref = db_xref.replace("taxon:", "") #removendo os characteres "taxon:"
    db_xref = str(db_xref)[2:-2]#removendo a aspa e o colchete
    file.write("%s\t%s\t%s\t%s" % (seq_record.id, organism, db_xref, seq_record.seq))
    file.write("\n")
    
file.close

    
import pandas as pd

df = pd.read_table('/Users/diegogotex/Desktop/MANEL/seqs/sequence.tab')
df = pd.DataFrame(file)
print(df)
    
    
    
    
    