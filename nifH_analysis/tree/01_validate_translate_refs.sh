# Part I: retreive AA sequences for NifH representatives

#consolidate Refs
cdhit-est -i ~/software/taxadiva/fDb.fasta -o nifHdb.cdhit.fa -c 0.8
perl -pi -e 's/\;|\.|-/_/g' nifHdb.cdhit.fa
perl -pi -e 's/\;|\.|-/_/g' nifHdb.cdhit.fa.clstr

#1966  finished        136  clusters

  # download full GB entry to find prot id for all nuc acc numbers
  # cat nifHdb.cdhit.fa | grep '>' |cut -f3 -d':' |cut -f1 -d';'|cut -f2 -d'>' 
  # grep 'protein_id' sequence.gb.txt |cut -f2 -d'"' > nifHprotacc
  # download protein fasta for all retrieved protein acc numbers
  # align in SeaView, remove outliers
  # = nifH.ncbi.vetted.fa 
  # perl -pi -e 's/-|\.|\[|\]//g' nifH.ncbi.vetted.fa
#translate Refs  (better annotations in fasta than NCBI peptides)
TransDecoder.LongOrfs -t nifHdb.cdhit.fa -S
blastp -query nifHdb.cdhit.fa.transdecoder_dir/longest_orfs.pep -subject nifH.ncbi.vetted.fa  -max_target_seqs 1 -outfmt 6 -evalue 1e-20 > blastp1
TransDecoder.Predict -t nifHdb.cdhit.fa --retain_blastp_hits blastp1 --single_best_orf
rm *gff3 *bed *cds; rm -rf *.transdecoder_dir 
