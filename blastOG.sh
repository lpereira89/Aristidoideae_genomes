#### For de novo genes and gene losses in BOP, we blast one gene per orthogroup (we don't take into account the paralogs)

# Directory with one file per OG, with all the sequences from all the species
seq=/shared/dunning_lab/Shared/Ahmed/OrthoFinder/primary_transcripts/OrthoFinder/Results_Oct31/Orthogroup_Sequences
wd=/mnt/parscratch/users/bo1lpg/Aristidoideae/2025Feb/

# Wrap the OG sequences fasta files from OrthoFinder
cat de_novo_OGs.txt | while read OG; do awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' \
    ${seq}/${OG}.fa > ${seq}/${OG}_wrapped.fa; done
cat gene_loss_OGs.txt | while read OG; do awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' \
    ${seq}/${OG}.fa > ${seq}/${OG}_wrapped.fa; done

# Go over the OGs and get one sequence per OG in a fasta file
mkdir de-novo-queries
cat de_novo_OGs.txt | while read OG; do sel_gene=$(cat ${seq}/${OG}.fa | grep ">" | shuf -n 1); \
    grep -A 1 ${sel_gene} ${seq}/${OG}_wrapped.fa > de-novo-queries/${OG}_random.fa; done
mkdir gene-loss-queries
cat gene_loss_OGs.txt | while read OG; do sel_gene=$(cat ${seq}/${OG}.fa | grep ">" | shuf -n 1); \
    grep -A 1 ${sel_gene} ${seq}/${OG}_wrapped.fa > gene-loss-queries/${OG}_random.fa; done

# Add the OG id into the fasta header
cat de_novo_OGs.txt | while read OG; do sed -i "s/>/>${OG}_/g" de-novo-queries/${OG}_random.fa; done
cat gene_loss_OGs.txt | while read OG; do sed -i "s/>/>${OG}_/g" gene-loss-queries/${OG}_random.fa; done

# Put all the genes in one fasta file
cd de-novo-queries
cat OG00* >> de-novo-queries.fa
cd ../gene-loss-queries
cat OG00* >> gene-loss-queries.fa

# Translate the cds using emboss (load environment)
transeq -sequence de-novo-queries.fa -outseq de-novo-queries-prot.fa
transeq -sequence gene-loss-queries.fa -outseq gene-loss-queries-prot.fa

# Make database with swissprot (download 17-02-2025)
wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot

# Blast them against swissprot (load environment)
blastp -query de-novo-queries-prot.fa -db ../database/uniprot_sprot.fasta -max_target_seqs 1 -outfmt "6 qseqid sallseqid sallgi qstart qend sstart send evalue bitscore salltitles" > de-novo-genes-blast1.txt
blastp -query gene-loss-queries-prot.fa -db ../database/uniprot_sprot.fasta -max_target_seqs 1 -outfmt "6 qseqid sallseqid sallgi qstart qend sstart send evalue bitscore salltitles" > gene-loss-genes-blast1.txt

# Extract maize genes for GO enrichment analysis #checked number of maize genes in the GeneCounts table and it matches
cat de_novo_OGs.txt | while read OG; do cat ${seq}/${OG}.fa | grep ">Zm" | shuf -n 1 | sed 's/>//g' >> maize_de-novo-genes.txt; done
cat gene_loss_OGs.txt | while read OG; do cat ${seq}/${OG}.fa | grep ">Zm" | shuf -n 1 | sed 's/>//g' >> maize_loss-genes.txt; done

#### For the duplications, we do it at the gene level (we take into account the paralogs)

# For each duplication event in PACMAD clade with >0.5 support, select one random gene from the gene1 column for blast # checked number and matches (482)
mkdir dup-queries
while read line; do
  OG=$(echo ${line} | cut -f 1 -d ' ');
  gene1=$(echo ${line} | sed 's/, /,/g' | cut -f 6 -d ' ' | sed 's/,/\n/g' | shuf -n 1);
  if [ -n "${gene1}" ]; then
    echo "${OG} ${gene1}" >> dup-queries/dup_random_genes.txt
  else
    gene2=$(echo ${line} | sed 's/, /,/g' | cut -f 7 -d ' ' | sed 's/,/\n/g' | shuf -n 1);
    echo "${OG} ${gene2}" >> dup-queries/dup_random_genes.txt;
  fi
done < filtered_N7.txt

# Take the species name added at the beginning of the line for each gene name
while read species; do
  sed -i "s/${species}//g" dup-queries/dup_random_genes.txt;
done < species.txt

# Wrap the sequences in the OG sequences file from OrthoFinder
cat dup-queries/dup_random_genes.txt | cut -f 1 -d ' ' | while read OG; do awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' \
    ${seq}/${OG}.fa > ${seq}/${OG}_wrapped.fa; done

# Get the cds sequence for the selected genes into a fasta file
cat dup-queries/dup_random_genes.txt | while read line; do OG=$(echo ${line} | cut -f 1 -d ' '); \
    sel_gene=$(echo ${line} | cut -f 2 -d ' '); \
    grep -A 1 ${sel_gene} ${seq}/${OG}_wrapped.fa >> dup-queries/${OG}_random.fa; done

# Put all the genes in one fasta file
####! Three extra sequences (485 instead of 482) -> two transcripts for D. radicosa genes, I took them out
cd dup-queries
cat OG00* >> dup-queries.fa
bioawk -c fastx '!($name == "Dra05G0194600.2" || $name == "Dra07G0172200.2" || $name == "Dra09G0389800.2") {print ">"$name"\n"$seq}' dup-queries.fa > dup-queries-uniq.fa

# Translate the cds using emboss (load environment)
transeq -sequence dup-queries-uniq.fa -outseq dup-queries-prot.fa

# Blast them against swissprot (load environment)
blastp -query dup-queries-prot.fa -db ../database/uniprot_sprot.fasta -max_target_seqs 1 -outfmt "6 qseqid sallseqid sallgi qstart qend sstart send evalue bitscore salltitles" > dup-genes-blast1.txt

# Add the OG id before the gene ID --> checked 1577 lines before and after

while read line; do
  OG=$(echo ${line} | cut -f 1 -d ' ');
  gene=$(echo ${line} | cut -f 2 -d ' ');
  # Escape special characters in gene for use with sed - if not it doesn't match the pattern
  gene_escaped=$(echo "${gene}" | sed 's/[&/\]/\\&/g')
  # Use sed to modify dup-genes-blast.txt and save the result to a temporary file
  #sed "s,${gene_escaped},${OG}_${gene},g" dup-genes-blast1.txt > tmpfile && mv tmpfile dup-genes-blast2.txt
  # Using temporary file was actually the problem! It substitutes one gene at a time but with every iteration it deletes the previous one
  sed -i "s,${gene_escaped},${OG}_${gene},g" dup-genes-blast2.txt
done < dup_random_genes.txt

# Take out the 'Phragmites_australis' species name from the gene ID
sed -i 's/Phragmites_australis_//g' dup-genes-blast1.txt
