## PART III:
#       use R to generate list of tips in chlorophyllide clade:
#       tre<-read.tree("RAxML_result.allseqs"); plot(tre); drop<-extract.clade(tre, 317); tre$tip.label[!(tre$tip.label %in% drop$tip.label)]
        #save good tips as list: goodnifH

#       BLAST nr (culturable bacteria) with 'good' novo nifH oligotypes (no chlorophyllides)
        #save to: nifH_nrhits.txt
        perl -pi -e 's/\ >.*$//g' nifH_nrhits.fa 
        perl -pi -e 's/\ /_/g' nifH_nrhits.fa 

# remake TREE, with nr hits, excluding chlorophyllide clade, with bootstraps
##  use bona fide nifH sequences, inferred from initial ML tree ('goodnifH')
##  by using sequences not in the clade with chlorophyllide reductases etc


cat nifH.aa.fa| perl -p -e 's/\|\\./_/g' | fasta_formatter -w 0 > temp
cat nifH_nrhits.fa <(grep -A 1 -f goodnifH temp) > nifH.aa.good.fa
perl -pi -e 's/--\n//g' nifH.aa.good.fa 

# re-align
mafft nifH.aa.good.fa > nifH.aa.good.aln

#re-tree
~/software/seqConverter.pl -dnifH.aa.good.aln -if -ope
#ML tree
raxmlHPC-PTHREADS-AVX -T 8 -s nifH.aa.good.phylip -m PROTGAMMAAUTO -n nifH_nr_good -p 1234
#boostraps, test for sufficient boot  sampling
raxmlHPC-PTHREADS-AVX -T 8 -s nifH.aa.good.phylip -m PROTGAMMALG -n boots -p 12345 -# 50 -x 1234
raxmlHPC-PTHREADS-AVX -T 8 -s nifH.aa.good.phylip -m PROTGAMMALG -n boots2 -p 12345 -# 50 -x 1234
raxmlHPC-PTHREADS-AVX -T 8 -s nifH.aa.good.phylip -m PROTGAMMALG -n boots3 -p 12345 -# 100 -x 1234
raxmlHPC-PTHREADS-AVX -T 8 -s nifH.aa.good.phylip -m PROTGAMMALG -n boots4 -p 12345 -# 300 -x 1234
cat RAxML_bootstrap* > allboots
raxmlHPC-PTHREADS-AVX -T 2 -m PROTGAMMALG -z allboots -I autoMRE -n testboots -p 1234
#map boot bipartions onto ML tree
raxmlHPC-PTHREADS-AVX -T 2 -m PROTGAMMALG -p 12345 -f b -t RAxML_bestTree.nifH_nr_good -z allboots -n bootslabel
