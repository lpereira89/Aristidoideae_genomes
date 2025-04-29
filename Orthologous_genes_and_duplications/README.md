## Identifying orthologous genes and gene duplications

We used Orthofinder to identify groups of orthologs and duplication events in grasses. We included 34 grass genomes in our analyses, including the two newly sequenced Aristidoideae genomes.

```./orthofinder -f primary_transcripts/ -S diamond -d -M msa -A muscle```

Then, using Orthofinder results, we selected the duplications originated at the base of the PACMAD clade with >0.5 support. We also defined a core set of duplicates that were present in all C4 species included in the analyses [script **DupAndRetGenesC4.py**].

For the core set, we built individual gene trees to verify the duplication [script **trees.sh**].
