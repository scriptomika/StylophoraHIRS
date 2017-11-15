

taxadiva.pl -s SampleListHIRSnifH.txt -o taxadiva_HIRS_4 \
-d /home/unhML/spank/software/taxadiva/fDb.fasta \
-t /home/unhML/spank/software/taxadiva/fTax.db.tsv \
--pear "-v 20 -m 450 -n 300 -p 1.0 -j 12" \
-y -k -r 26 -l 29 -j 12 --keepc4 -g 500 --med-metadata MEDmetafile.txt

#if fails at oligotype step because MEDmetafile incorrect sample names
# fix and re-run manually at failed command:
#decompose -o taxadiva_HIRS_4/MED -E MEDmetafile.txt taxadiva_HIRS_4/MED/taxadiva_HIRS_4.reads.med.withPadgaps
