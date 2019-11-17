# README

## dmr_pw_counts_pe.py
Does pairwise comparisons of methylation from allC files for a given set of regions and outputs to a file; results include per-sample number of methylated reads, unmethylated reads, and weighted methylation for each region

## dmr_pw_ztesting_pe.py
Computes the one-sided z-test pairwise between branches to identify significant differences in methylation. Uses the result of `dmr_counts_pe.py` as the input

## bed_gene_overlap.py
Assigns genomic feature to regions

## region_haplo_filter.py
Filter regions based on read coverage and number of covered positions in the region; remaining regions were used to determine methylation pseudo-alleles

## region_haplo_allele_v4.py
Determine the methylation haplotype of regions for a given sample