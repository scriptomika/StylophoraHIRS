symTyper.py --verbose -t 8 subtype -s all.samples.ids -H ~/MPL/SymTyper.full2/hmmer_hits \
-b BLASTOUTDIR -r blastResults \
-f fasta --blastdb ~/software/symTyper/database/blast_DB/ITS2_Database_04_23_13.fas

symTyper.py --verbose -t 8 resolveMultipleHits -s all.samples.ids \
-m ~/MPL/SymTyper.full2/blastResults/MULTIPLE -c resolved_clusters 

symTyper.msp.py buildPlacementTree -c resolved_clusters/correctedMultiplesHits/corrected \
-n /home/unhML/katm/software/symTyper/database/clades_phylogenies/ -o placements

#output subclade counts per sample as tsv file in placements/C/, placements/D/ etc
