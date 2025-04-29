# Use the maize genes to blast sequences and manually construct gene trees

#unwrap FASTA
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' dbs/Zea_mays.fa > Zea_mays-wrapped.fa

#save CDS sequence for each gene in its own file
cat Zmgenes_fortrees.txt | while read line; do grep -A 1 ${line} Zea_mays-wrapped.fa >> ${line}.fa; done

#merge all the CDS files for all the species and construct a merged database
source activate blast
cd dbs
cat *.fa > merged_db.fasta
makeblastdb -in merged_db.fasta -dbtype nucl

#run blastn for the maize genes and filter by match length = 300
cd ..
mkdir blast_results
cat Zmgenes_fortrees.txt | while read gene; do blastn -query ${gene}.fa -db dbs/merged_db.fasta \
  -outfmt '6 sseqid length sseq pident evalue' | awk '$2 >= 300' > blast_results/${gene}_blast_filtered.txt ; done

#select the match sequence for each blast hit and put them in a fasta file
cat ${wd}/Zmgenes_fortrees.txt | while read gene; do cat ${wd}/blast_results/${gene}_blast_filtered.txt | while read match; do \
  echo ${match} | awk '{print ">" $1 "\n" $3}' >> ${wd}/files_to_align/${gene}.fa; done; done

#merge the files that correspond to the same OG
cat /mnt/parscratch/users/bo1lpg/Aristidoideae/OG.list | while read OG; do cat ../${OG}.txt | sed 's/Zea_mays_//g' | while read gene; do \
  cat files_to_align/${gene}.fa >> files_to_align/${OG}.fa; done; done
cat OG.list | while read OG; do sed -i 's/-//g' files_to_align/${OG}.fa; done

#modify names because some don't work later for the tree
sed -i -E 's/Phragmites_australis_PaChr[0-9]+_[0-9]+_[0-9]+_evm.model.//g' OG0*_filt.fa
sed -i -E 's/Dendrocalamus_sinicus_Dsi[0-9]+[ABC]_[0-9]+_[0-9]+_//g' OG0*_filt.fa
sed -i -E 's/Dendrocalamus_sinicus_DsiCtg[0-9]+[a-z]?_[0-9]+_[0-9]+_//g' OG0*_filt.fa
sed -i -E 's/_-?1//g' OG0*_filt.fa
sed -i -E 's/:/_/g' OG0*_filt.fa

# Use python to select only the longest match per gene --> script keep_longest_seq.py
#align the sequences - run in batch mode
source activate mafft
mkdir alignments
cat OG.list | while read OG; do cat ../${OG}.txt | sed 's/Zea_mays_//g' | while read gene; do \
  mafft --thread -1 --add --auto files_to_align/${OG}.fa queries/${gene} > alignments/${OG}_aln.fa; done; done

#Change the alignment format from fasta to phylip
#Location of the script in the following line, but not working in the cluster - copy to the wd
#fasta_to_phylip=/shared/dunning_lab/Shared/scripts/perl/Fasta2Phylip.pl
cat OG.list | while read OG; do perl ../trees/Fasta2Phylip.pl alignments/${OG}_aln.fa alignments/${OG}_aln.phy; done

#Construct trees
#Location of the script in the following line, but not working in the cluster - copy to the wd
#sms=/shared/dunning_lab/Shared/programs/sms-1.8.1/./sms.sh
mkdir trees
cd trees
#Very slow, do it with an SGE job array
cat OG.list | while read OG; do ../sms-1.8.1/./sms.sh -i ../alignments/${OG}_aln.phy -d nt -t -b 100; done
