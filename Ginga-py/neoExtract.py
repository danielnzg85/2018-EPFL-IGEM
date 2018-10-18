#!/usr/bin/env python
# coding: utf-8

# In[1]:


#neoExtract -  sorts the amino-acid sequences etracted from annovar gene annotation and only contains the mutated amino-acid of specify length
#Typically Neoantigens are 8-11 ammino acid long although can often reach up to 14-mer
#This script inputs a fasta file containing the coding sequences and outputs another fasta file with the sequences 
#WILDTYPE proteins are removed
#Lines starting with '>' are mantain as the header of the protein sequence

from pylab import *
import re
import numpy as np
import sys
import argparse


input_file = 'filename'
output_file = 'outputfile'
length_peptide = 'length'

filename=open(input_file,'r');
m=filename.read();
filename.close();
#input_invalid = 1
#while input_invalid:
#    length_peptide = int(float(input("Please enter maximum amino acid length of the neoantigens (from 8 to 14 mer): ")))
#    if length_peptide <= 14:
#        input_invalid = 0
#        length_peptide = length_peptide-1
#    else:
#        print("Please enter a valid integer")



#Divide into headers
d=m.replace('\n','').replace('>','@>');
D={};
D=d.split('@');
del D[0]

#Remove WildType Sequences
D1 = [];
for k in range(0,len(D)):
    if not re.search('WILDTYPE', D[k]):
        D1.append(D[k])
    
#Separete sequence from header    
C3=[];
for i in range (0, len(D1)):
    p = re.split('\)', D1[i])
    p[0] = p[0]+')'
    C3.append(p[0])
    C3.append(p[1])


index_mut = [];
#Separate depending on the type of mutation and extract indexes of mutation
for i in range (0, len(C3)):
    if re.search('insertion', C3[i]):
        result = re.search('position (.*) has', C3[i])
        index_mut.append(result.group(1))
    elif re.search('changed', C3[i]):
        result = re.search('position (.*) changed', C3[i])
        index_mut.append(result.group(1))
    elif re.search('amino', C3[i]):
        result = re.search('position (.*) changed', C3[i])
        index_mut.append('0')
        

#Separate initial and end indexes and create master matrix
master_matrix=[];
p=0;
rt=0;
for i in range (0,len(index_mut)):
    if re.search('-',index_mut[i]):    
        h1 = re.split('-', index_mut[i])
    else:
        rt = 1+rt;
        h1 = [index_mut[i],'0'];
    master_matrix.append([h1[0], h1[1], C3[p],C3[p+1]]);
    p=p+2;
    
    
#Sort peptide sequences accordingly
just_seq =[];
final_seq=[];

for i in range (0, len(master_matrix)):
    if len(master_matrix[i][2]) >= 2000:
        line_index = master_matrix[i][2][0:2000];
    else:
        line_index = master_matrix[i][2];
        
    neat_seq = master_matrix[i][3];
    if re.search('to',master_matrix[i][2]):
        prem_seq=re.search(' to (.*)\)', master_matrix[i][2])
        if prem_seq:
            exact_seq = prem_seq.group(1)
        else:
            exact_seq = '';
    else:
        prem_seq=re.search('insertion (.*)\)', master_matrix[i][2])
        if prem_seq:
            exact_seq = prem_seq.group(1)
        else:
            exact_seq='';

    if int(float(master_matrix[i][0])) == 0:
        cut_seq = master_matrix[i][3];
    elif int(float(master_matrix[i][1])) == 0:     
        cut_seq = master_matrix[i][3];
        if re.search('immediate-stoploss',master_matrix[i][2]):
            cut_seq = neat_seq[len(neat_seq)-length_peptide-len(exact_seq):len(neat_seq)]
        elif float(master_matrix[i][0]) > len(neat_seq):
            cut_seq = neat_seq;
        elif (float(master_matrix[i][0])-length_peptide) < 1 and (float(master_matrix[i][0])+length_peptide+len(exact_seq)) > len(neat_seq):
            cut_seq = neat_seq[0:len(neat_seq)]
        elif float(master_matrix[i][0])-length_peptide < 1:
            cut_seq = neat_seq[0:int(float(master_matrix[i][0]))+length_peptide+len(exact_seq)]
        elif float(master_matrix[i][0])+length_peptide+len(exact_seq) > len(neat_seq):
            cut_seq = neat_seq[int(float(master_matrix[i][0]))-length_peptide:len(neat_seq)]
        else:
            cut_seq = neat_seq[int(float(master_matrix[i][0]))-length_peptide+1:int(float(master_matrix[i][0]))+length_peptide+len(exact_seq)]
    else:
        if len(neat_seq)==len(exact_seq):
            cut_seq = neat_seq;
        elif float(master_matrix[i][0])-length_peptide-1 < 1 and  float(master_matrix[i][1])+length_peptide > len(neat_seq):
            cut_seq = neat_seq[0:len(neat_seq)]
        elif float(master_matrix[i][0])-length_peptide-1 < 1:
            cut_seq = neat_seq[0:int(float(master_matrix[i][1]))+length_peptide]
        elif float(master_matrix[i][1])+length_peptide > len(neat_seq):
            cut_seq = neat_seq[int(float(master_matrix[i][0]))-length_peptide:len(neat_seq)]
        elif float(master_matrix[i][1])+length_peptide+len(exact_seq) > len(neat_seq):
            cut_seq = neat_seq[int(float(master_matrix[i][0]))-length_peptide:int(float(master_matrix[i][0]))-1+len(exact_seq)]
        else:
            cut_seq = neat_seq[int(float(master_matrix[i][0]))-length_peptide:int(float(master_matrix[i][0]))+len(exact_seq)+length_peptide]
    
    final_seq.append(line_index)
    final_seq.append(cut_seq)
    just_seq.append(cut_seq)
    
    
#Calculate average length of the peptides just for validation
length_index = [];
for i in range (1,len(just_seq)):
    length_index.append(len(just_seq[i]))
ave_length= sum(length_index)/len(length_index)

#Write output to file
with open(output_file, 'w') as f:
    for item in final_seq:
        f.write("%s\n" % item)


# In[ ]:




