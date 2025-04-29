# Extract maize genes for GO enrichment analysis - one gene per duplication event (in some cases, several per OG)
while read line; do
  OG=$(echo ${line} | cut -f 1 -d ' ');
  gene1=$(echo ${line} | sed 's/, /,/g' | cut -f 6 -d ' ' | sed 's/,/\n/g' | grep "Zea_mays_" | shuf -n 1);
  if [ -n "${gene1}" ]; then
    echo "${OG} Gene1=${gene1}" >> maize_dup_genes.txt;
  else
    gene2=$(echo ${line} | sed 's/, /,/g' | cut -f 7 -d ' ' | sed 's/,/\n/g' | grep "Zea_mays_" | shuf -n 1);
    echo "${OG} Gene2=${gene2}" >> maize_dup_genes.txt;
  fi
done < filtered_N7.txt
grep "Zea_mays_" maize_dup_genes.txt | cut -f 2 -d ' ' | cut -f 3 -d '_' > maize_dup_genes_IDs.txt

#### To analyse expression patterns of candidate genes, we need to extract gene IDs for sorghum, maize and rice

mkdir gene-expression
cd gene-expression

# For duplicated genes, we get gene lists from maize and sorghum from the duplications file. All genes for each duplication event
mkdir dups
cd dups

while read line; do
  OG=$(echo ${line} | cut -f 1 -d ' ');
  gene1=$(echo ${line} | sed 's/, /,/g' | cut -f 6 -d ' ' | sed 's/,/\n/g' | grep "Zea_mays_");
  gene2=$(echo ${line} | sed 's/, /,/g' | cut -f 7 -d ' ' | sed 's/,/\n/g' | grep "Zea_mays_");
  echo "${gene1}" >> maize_${OG}.txt;
  echo "${gene2}" >> maize_${OG}.txt;
done < ../../filtered_N7.txt
ls maize_OG* | while read file; do sed -i 's/Zea_mays_//g' ${file}; done
# Merge the individual OG files and keep only unique records
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do
  cat maize_${OG}.txt | while read gene; do
    echo -e "${OG}\t${gene}" >> maize_dups_temp.txt;
  done
done
cat maize_dups_temp.txt | sort | uniq > maize_dups.txt

while read line; do
  OG=$(echo ${line} | cut -f 1 -d ' ');
  gene1=$(echo ${line} | sed 's/, /,/g' | cut -f 6 -d ' ' | sed 's/,/\n/g' | grep "Sorghum_bicolor_");
  gene2=$(echo ${line} | sed 's/, /,/g' | cut -f 7 -d ' ' | sed 's/,/\n/g' | grep "Sorghum_bicolor_");
  echo "${gene1}" >> sorghum_${OG}.txt;
  echo "${gene2}" >> sorghum_${OG}.txt;
done < ../../filtered_N7.txt
ls sorghum_OG* | while read file; do sed -i 's/Sorghum_bicolor_//g' ${file}; done
# Merge the individual OG files and keep only unique records
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do
  cat sorghum_${OG}.txt | while read gene; do
    echo -e "${OG}\t${gene}" >> sorghum_dups_temp.txt;
  done
done
cat sorghum_dups_temp.txt | sort | uniq > sorghum_dups.txt

# Search for the Orthogroups that do not have a maize/sorghum representative
awk '{
    if ($2 != "") {          # If the second column (gene_id) is not empty
        seen[$1] = 1;        # Mark the orthogroup (first column) as having a gene
    } else {
        missing[$1] = 1;     # Mark the orthogroup (first column) as missing a gene
    }
}
END {
    # Now print out the orthogroups that have no gene ID
    for (og in missing) {
        if (!(og in seen)) {
            print og;        # Print orthogroups with no associated gene ID
        }
    }
}' maize_dups.txt

# For rice genes, we can search for the paralogs present at the orthogroup level
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do cat ${seq}/${OG}.fa | grep ">LOC_Os" | sed 's/>//g' >> rice_${OG}.txt; done
cat rice_OG* > rice_dups.txt
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do
  cat rice_${OG}.txt | while read gene; do
    echo -e "${OG}\t${gene}" >> rice_dups.txt;
  done
done
rm sorghum_dups_temp.txt maize_dups_temp.txt
cd ..

# For de novo genes and gene losses, we get them at the orthogroup level from the orthogroup sequences files only for maize and sorghum (absent in rice)
mkdir losses
cd losses
# First create one file per OG
cat ../../gene_loss_OGs.txt | while read OG; do cat ${seq}/${OG}.fa | grep ">Zm" | sed 's/>//g' >> maize_${OG}.txt; done
cat ../../gene_loss_OGs.txt | while read OG; do cat ${seq}/${OG}.fa | grep ">SORBI" | sed 's/>//g' >> sorghum_${OG}.txt; done
# Then merge them
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do
  cat maize_${OG}.txt | while read gene; do
    echo -e "${OG}\t${gene}" >> maize_gene_loss.txt;
  done
done
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do
  cat sorghum_${OG}.txt | while read gene; do
    echo -e "${OG}\t${gene}" >> sorghum_gene_loss.txt;
  done
done
cd ..

mkdir de-novo
cd de-novo
# First create one file per OG
cat ../../de_novo_OGs.txt | while read OG; do cat ${seq}/${OG}.fa | grep ">Zm" | sed 's/>//g' >> maize_${OG}.txt; done
cat ../../de_novo_OGs.txt | while read OG; do cat ${seq}/${OG}.fa | grep ">SORBI" | sed 's/>//g' >> sorghum_${OG}.txt; done
# Then merge them
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do
  cat maize_${OG}.txt | while read gene; do
    echo -e "${OG}\t${gene}" >> maize_de_novo.txt;
  done
done
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do
  cat sorghum_${OG}.txt | while read gene; do
    echo -e "${OG}\t${gene}" >> sorghum_de_novo.txt;
  done
done

#### Once the GO term overrepresentation test from PANTHER is done, extract the GO term IDs
sed -n 's/.*(GO:\([0-9]*\)).*/GO:\1/p' Overrepr_dup.txt > Overrepr_dup_clean.txt
sed -n 's/.*(GO:\([0-9]*\)).*/GO:\1/p' Underrepr_dup.txt > Underrepr_dup_clean.txt
sed -n 's/.*(GO:\([0-9]*\)).*/GO:\1/p' Underrepr_denovo.txt > Underrepr_denovo_clean.txt

#### Get the maize gene ID for the Orthogroups (for losses and de novo genes) that were retained in all C4 species
mkdir retained
cd retained
cat OG_retained.txt | while read OG; do cat ${seq}/${OG}.fa | grep ">Zm" | sed 's/>//g' > maize_${OG}.txt; done
ls maize_OG* | sed 's/maize_//g' | sed 's/.txt//g' | while read OG; do
  gene=$(tr '\n' ',' < maize_${OG}.txt);
  echo -e "${OG}\t${gene}" >> maize_retained.txt
done
sed -i 's/,$//g' maize_retained.txt
