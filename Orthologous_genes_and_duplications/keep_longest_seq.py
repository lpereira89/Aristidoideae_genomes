#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

os.chdir('/Users/larapereiragarcia/Desktop/Sheff/Aristidoideae-genomes/Aristidoideae-resultsNov/trees/trees_onlymatch')

list_files='OG.list'

# Put the fasta files to be processed in a list
with open(list_files, 'r') as f:
    OG = [line.rstrip('\n') for line in f]

for ID in OG:
    file_name = 'OG_seqs/' + ID + '.fa'
    # For each file, generate a list with all the records using SeqRecord Biopython class
    records = []
    for seq_record in SeqIO.parse(file_name, "fasta"):
        records.append(seq_record)
    # Initialize a variable to count the number of initial records, and two dictionaries to store the sequence length and the record index
    count = 0
    seq_length = {}
    seq_index = {}
    for index in range(len(records)):
        # If this sequence is already contained in the dictionary, compare sequence length and either save the index (if longer) or ignore (if shorter)
        if records[index].id in seq_length.keys():
            if len(records[index].seq) > seq_length[records[index].id]:
                count += 1
                seq_length[records[index].id] = len(records[index].seq)
                seq_index[records[index].id] = index
                #print(f'Seq {records[index].id} included because its length {len(records[index].seq)} is greater than previous records ({seq_length[records[index].id]})')
            else:
                count += 1
                #print(f'Seq {records[index].id} discarded because its length {len(records[index].seq)} is smaller than previous records ({seq_length[records[index].id]})')
        # If it's the first time that sequence appear, store it as 'selected' in the seq_index dictionary and store its lenght to compare with other matches if any
        else:
            count += 1
            seq_length[records[index].id] = len(records[index].seq)
            seq_index[records[index].id] = index
            #print(f'Seq {records[index].id} with length {len(records[index].seq)} included - first ocurrence')
    output_file = ID + '_filt.fa'
    print(f'{len(seq_index)} out of {count} to be exported to filtered fasta file')
    with open(output_file, 'w') as handle:
        # Use the stored indexes in the dictionary seq_index to write the largest ocurrence of each sequence in a new fasta file
        for record_index in seq_index.values():
            SeqIO.write(records[record_index], handle, 'fasta')
