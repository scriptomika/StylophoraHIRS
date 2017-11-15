
# use modified symTyper.py script to take folder of sample fastas instead 
# of single concatenated file (reducing run time by several hours)
symTyper.msp.py --verbose -t 16 clade -f fasta -s all.samples.ids \
--hmmdb  ~/software/symTyper/database/HMMER_ITS2_DB/All_Clades.hmm

#counts for each clade assignment per sample in hmmer_parsedOutput/DETAILED_counts.tsv
