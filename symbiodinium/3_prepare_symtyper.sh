#!/bin/bash
# usage: ./prepare_symptyper.sh 

# input requirements for symTyper 
# sample.ids: file listing unique id per sample(sampleID)
# fasta: single fasta with all samples' sequences. header: sampleID::number


# generate sample ID file
for f in $(ls /home/unhML/katm/MPL/non_chimeras/*.fasta)  #directory with fastas
do
name=${f%_good.nonchimera.fasta}   # common suffix to remove
sample=$(basename $name)
echo $sample >> samples_symtyper.ids

if [ ! -e "${f}.renum" ]; then
#renumber fasta headers with sampleid
count=0
while read line; do
if [[ "$line" == ">"* ]];then
((count++))
newline=`echo $line |sed -e "s/>.*/>$sample::$count/"`
else
newline=$line
fi
echo $newline >> ${f}.renum
done < ${f}
fi

done
#for orignal symTyper script, must generate  single fasta input
cat /home/unhML/katm/MPL/non_chimeras/*renum > samples_symtyper.fasta

## for modified msp script, must provide all fastas (ending with *fasta) in a folder called 'fasta'
## cp *renum files to new folder, then go there and run:
#rename 's/_good\.nonchimera\.fasta\.renum$/\.fasta/' *.renum

