## Gene expression of candidate C4 precursors

We used previously published datasets to check whether the candidate C4 precursors were differentially expressed between mesophyll and bundle sheath cells in two C4 species (maize and sorghum) and one C3 species (rice). We first extracted the maize, sorghum and rice gene IDs from our candidate orthogroups [script **ExtractGeneIDs.sh**]. 

Then we used R to find these genes in the original datasets and tested whether genes duplicated at the base of PACMAD were more likely to be differentially expressed between mesophyll and bundle sheath cells than genome-wide [script **Gene_expression.R**].
