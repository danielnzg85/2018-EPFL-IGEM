#!/usr/bin/env python
# coding: utf-8

# In[4]:


#neoRANK - ranks the output of netMHC according to their rank and outputs only the strong binders and weakbinders for each allele type. 
#The input file of the code is a XLS file created by netMHC tool.
#The software outputs the results in a more user fiendly TXT file, with the header containing the HLA type.
#More information about usage and options can be 

__author__ = "2018 EPFL IGEM Team (Develop by Daniel Nakhaee-Zadeh Gutierrez)"
__license__ = "Public Domain"
__version__ = "1.0"

from pylab import *
import re
import numpy as np
import sys
import pandas as pd
import xlrd
import csv
import argparse

input_file = 'NetMHC_out_66.xlsx'
output_file = 'neoSearch_66.csv'

#Open xls file
df = xlrd.open_workbook(input_file)
df4 = df.sheets()[0]

rows_n = df4.nrows
cols_n = df4.ncols

target_col =[];
for i in range(df4.nrows):
    row = df4.row(i)
    for p, cell in enumerate(row):
        if cell.value == "Rank" :
            target_col.append(p)
            

heads = df4.row_values(0, start_colx=0, end_colx=None)
header_index = [];
heads_new =[];
for i in range (0,len(heads)):
    if heads[i]:
        header_index.append(i)
        heads_new.append(heads[i])

data = [df4.row_values(i) for i in range(df4.nrows)]
headers = data[1:2]
data = data[2:]

lines = []
for i in range (0,len(data)):
    lines.append(data[i][0:3])

#Adjust lines and replicate lines
master=[];
index_head=[];
for i in range(0,len(header_index)):
    heads_adj_len = len(lines[0])+ header_index[1]-header_index[0];
    heads_adj = ['']*heads_adj_len;
    heads_adj.insert(0,heads_new[i])
    index_head.append(len(master))
    master.append(heads_adj)
    master.append(headers[0][0:heads_adj_len])
    for p in range (0,len(data)):
        master.append(lines[p]+data[p][header_index[i]:header_index[i]+(header_index[1]-header_index[0])])
        
 
#Sort rows for each HLA type
master_sorted =[];
for g in range (0,len(index_head)):
    if g == len(index_head)-1:
        master_tmp=master[index_head[g]+2:]
    else:
        master_tmp=master[index_head[g]+2:index_head[g+1]]
    index_rank = master[index_head[g]+1].index("Rank")
    master_tmp.sort(key=lambda x: x[index_rank])
    master_sorted.append(master[index_head[g]:index_head[g]+2]+master_tmp[0:])
    
    
#Delete not binding peptides and add SB and WK label
peptide_f = []
for d in range (0,len(master_sorted)):
    for i in range (0,len(master_sorted[0])):
        if i in index_head or i in [x+1 for x in index_head]:
            peptide_f.append(master_sorted[d][i]+[''])
        elif 0 < master_sorted[d][i][index_rank] < 0.5:
            peptide_f.append(master_sorted[d][i]+['SB'])
        elif 0.5 < master_sorted[d][i][index_rank] < 2:
            peptide_f.append(master_sorted[d][i]+['WB'])

#Create a file with output sequences        
with open(output_file, 'w') as f:
   writer = csv.writer(f, delimiter=',')
   writer.writerows(peptide_f)


# In[ ]:





# In[ ]:




