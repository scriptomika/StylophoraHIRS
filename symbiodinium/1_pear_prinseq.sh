
#!/bin/bash
# Use PEAR to merge the paired-end reads and use prinseq to quality filter them.
# usage: ./1_pear_prinseq.sh <directoryname>
# reqs: specify directory on command line; directory contains subdirectories with fastq.gz read files

#set up log file to record discards
echo Library$'\t'PEAR.ReadsDiscarded$'\t'PEAR.R1unassembled$'\t'PEAR.R2unassembled$'\t'PrinSeq.badreads > ${1}.pear_prinseq_discard.log

for dir in ${1}/*
do

	mkdir -p ./PEAR_merged
	mkdir -p ./prinseq_filtered
	if [ "$(ls -A $dir)" ]; then
		ffilez=$(ls ${dir}/*R1_001.fastq.gz)
		rfilez=$(ls ${dir}/*R2_001.fastq.gz)
		#echo $ffilez; echo $rfilez
		out1=${ffilez%_L00*}
		out=$(basename $out1)

		if [ ! -e "./prinseq_filtered/${out}_good.fasta" ]; then
			# Decompress raw data
			gunzip $ffilez
			gunzip $rfilez

			# Loop through the sample forward and reverse reads files and use PEAR to merge them
			pear-0.9.10-bin-64 -f ${ffilez%.gz} -r ${rfilez%.gz} \
			-o ./PEAR_merged/${out}.merged.fastq -p 1.0 -m 450 -n 250 -y 4000M -j 4 | tee -a PEAR_stringent_log.txt

			run=$(grep -c '^+$' ./PEAR_merged/${out}.merged.fastq.unassembled.reverse.fastq)
			fun=$(grep -c '^+$' ./PEAR_merged/${out}.merged.fastq.unassembled.forward.fastq)
			dis=$(grep -c '^+$' ./PEAR_merged/${out}.merged.fastq.discarded.fastq)

			rm -f ./PEAR_merged/${out}.merged.fastq.unassembled* 
			rm -f ./PEAR_merged/${out}.merged.fastq.discarded* 

			# Loop through the sample files of merged sequences and use prinseq to quality filter them
			prinseq-lite -fastq ./PEAR_merged/${out}.merged.fastq.assembled.fastq \
			-out_format 1 -out_good ./prinseq_filtered/${out}_good -out_bad ./prinseq_filtered/${out}_bad \
			-min_len 250 -max_len 450 -min_qual_score 20 -min_qual_mean 30 -noniupac -log

			bad=$(grep -c '>' ./prinseq_filtered/${out}_bad.fasta)

			#record discard stats from PEAR and prinseq
			echo $out $dis $fun $run $bad | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5}' >> ${1}.pear_prinseq_discard.log

			# Recompress raw data for storage
			gzip ${ffilez%.gz}; gzip ${rfilez%.gz}

			#chimera removal step

			#symtyper
		fi
	fi
done
