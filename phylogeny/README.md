## Description

In order to determime whether AMPs are expressed in similar levels, they were characterized by Expression Levels and a phylogenetic tree was built based on 356 sequences.

## Workflow 
  1. The AMP Fasta file is first processed to simplify the fasta headers
  2. Multiple Sequence Alignment is performed via MAFFT
  3. Newick File is created from the aligned fasta file using MEGA X
  4. [R script (ggtree,ape,tidyverse)](phylogeny.R) is used to visualize the tree 
