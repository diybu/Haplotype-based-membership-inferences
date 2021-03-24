# Haplotype-based-membership-inferences
The source code to conduct haplotype-based membership inferences (KHF-inference, KHB-inference and NHF-inference).

## Prerequisites
The vcf files and panel file from [1000 genomes project](https://www.internationalgenome.org/data#download) is used as the original data.

The haplotypes from the original data are obtained by the softerware [Haploview](https://www.broadinstitute.org/haploview/haploview) 

## Getting started
1. collect_haplotype: collect known haplotypes from 1000 genomes and select the target haplotypes.

2. build_private_database: randomly extract a group of genomes to build a private database, save the frequency table and beacon table.

## Performing membership inferences
KHF-inference: infer known haplotypes from the frequency table.

KHB-inference: infer known haplotypes from the beacon table.

NHF-inference: reconstruct novel (unknown) haplotypes from the frequency table.
