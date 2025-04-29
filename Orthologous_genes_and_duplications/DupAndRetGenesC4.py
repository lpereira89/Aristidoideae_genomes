#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 11:46:26 2024

@author: larapereiragarcia
"""
import os
print(os.getcwd())
os.chdir('/Users/larapereiragarcia/Desktop/Sheff/Aristidoideae-genomes/Aristidoideae-resultsNov')
print(os.getcwd())


input_file='N7.tsv'
species_file='species-C4.txt'
output_file='dupANDretained_C4.txt'
maize_file='maize_genes.txt'

with open(species_file, 'r') as f:
    species = [line.rstrip('\n') for line in f]

output = open(output_file,"w+")
genes_tokeep = open(maize_file,"w+")

with open(input_file, 'r') as handle:
    for duplication in handle:
        dup=duplication.split('\t')
        OG=dup[0]
        genes1=dup[5]
        genes2=dup[6]
        count1=0
        count2=0
        for s in species:
            if s in genes1:
                count1 += 1
        for s in species:
            if s in genes2:
                count2 += 1
        if count1 == 15 and count2 == 15:
            output.writelines(OG + ';' + genes1 + ';' + genes2)
            gene_list1=genes1.split(',')
            gene_list2=genes2.split(',')
            maize_genes=[]
            print(OG)
            for gene1 in gene_list1:
                if 'Zea' in gene1:
                    maize_genes.append(gene1)
            for gene2 in gene_list2:
                if 'Zea' in gene2:
                    maize_genes.append(gene2)
            genes_tokeep.writelines(OG + '\t' + '\t'.join(str(gene) for gene in maize_genes) + '\n')
            
                

    