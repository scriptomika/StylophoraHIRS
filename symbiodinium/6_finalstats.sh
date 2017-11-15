#cat fasta/*fasta > all.samples.fasta
symTyper.py stats -i all.samples.fasta --outputs_dir . --out_file ./stats_output
symTyper.py makeTSV  --outputs_dir .
