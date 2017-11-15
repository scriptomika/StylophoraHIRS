# Part II: Predict ORFs, get AA for new Oligotypes and well-labeled DB seqs
#translate oligotypes
cp ../taxadiva_HIRS_4/MED/NODE-REPRESENTATIVES.fasta .
perl -pi -e 's/\|size.*$//g' NODE-REPRESENTATIVES.fasta
cat NODE-REPRESENTATIVES.fasta | perl -p -e 's/-//g' | perl -p -e 's/>00000/>nifH_novo_/g' > NODES.fasta; 
cdhit-est -i NODES.fasta -o NODES.cdhit.fasta -c 0.95
TransDecoder.LongOrfs -t NODES.cdhit.fasta 

#blast to predict correct ORF
blastp -query NODES.cdhit.fasta.transdecoder_dir/longest_orfs.pep -subject nifHdb.cdhit.fa.transdecoder.pep  -max_target_seqs 1 -outfmt 6 -evalue 1e-20 > blastp2
TransDecoder.Predict -t NODES.cdhit.fasta --retain_blastp_hits blastp2 --single_best_orf
rm *gff3 *bed *cds; rm -rf *.transdecoder_dir 

#select only AA seqs from oligotyping that had a blast match
cut -f3 -d':' blastp2 > seqs2keep
selectSeqs.pl -p -f seqs2keep NODES.cdhit.fasta.transdecoder.pep > new_oligotypes.aa.fas
rm seqs2keep

#align new oligotypes with database nifH
cat nifHdb.cdhit.fa.transdecoder.pep new_oligotypes.aa.fas | perl -p -e 's/::g.*$|XXXXXXXXXXXXXXXXXX|\*//g' | perl -p -e 's/>Gene\.[0-9]+::/>/g' > nifH.aa.fa

#combine with chlorophyllide/ferrodoxin reductases to id non-nifH oligotypes
cat chlorophyllides.fa nifH.aa.fa > allseqs.aa.fa

mafft allseqs.aa.fa > allseqs.aa.aln
~/software/seqConverter.pl -dallseqs.aa.aln -if -ope

#tree
raxmlHPC-PTHREADS-AVX -T 8 -s allseqs.aa.phylip -m PROTGAMMAAUTO -n allseqs -p 1234 

