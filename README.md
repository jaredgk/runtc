# msh-python

Compute tc/alpha values for SNPs in a given VCF file. There are three components to this process: reversing the original VCF file,
generating the lengths of the most shared haplotypes using the Burrows-Wheeler transform, and using those lengths to estimate tc. 
