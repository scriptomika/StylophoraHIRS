# StylophoraHIRS

This repository contains the scripts used to analyze diazotroph (nifH) and *Symbiodinium* (ITS) diversity in *Stylophora pistillata*.

Raw data is available at NCBI SRA for [BioProject PRJNA390752](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA390752)

##NifH analysis contains:
*shell script to run TaxaDiva pipeline. nifH database and executables available from [TaxaDiva github](https://github.com/lavanyarishishwar/taxadiva/blob/master/README.md) 

*Metadata info for each sample

*phylogenetic analyses. Scripts require taxadiva nifH database. Selected outputs are provided:
  -allseqs.aa.phylip : alignment of all novel predicted NifH and taxadiva reference sequences. Also includes chlorophyllide reductases (Cluster V nifH-like paralogs)
  -RAxML_bestTree.allseqs: ML phylogeny of all novel + reference NifH sequences
  -RAxML_bipartitions.bootslabel: ML phylogeny of novel + reference NifH sequences from clusters I-III, with bootstrap supports


*R script for statistical analysis and plots

## *Symbiodinium* ITS analysis files
* shell scripts to run SymTyper pipeline.
  -SymTyper pipeline executables available [here](https://github.com/UH-Bioinformatics/symTyper/tree/master/commands)
  -includes modified version of core symtyper.py script

* R script to analyze *Symbiodinium* phylotype counts from *Stylophora* samples
  -selected outputs from SymTyper included
