# RUNTC

Generate estimates of first coalescent time, tc.  

This documentation covers the usage of runtc.py and phase_singletons.py.

## VCF Requirements

Typically a vcf file should contain phased genotypes at all variable positions of a sample of individuals' genomes for a single chromosome. 

If the file has not been phased at singleton sites,  then it should be first run through phase_singletons.py

For this analysis, all sites with missing data, multiple derived alleles, or no variation will be ignored. 


## phase_singletons.py

This program takes as its main argument a vcf or vcf.gz file and returns to stdout a new file with new placements for singleton alleles as estimated from msh lengths. 
In the course of its operation it also generates two temporary msh (or msh.gz) files which can be saved if desired, but are deleted by default. 

#### Command line usage
##### positional arguments:
* vcfname              vcf or vcf.gz input file

##### optional arguments:
*  -h, --help           show this help message and exit
*  --gzip               gzip temporary msh files
*  --keep               save temporary msh and reversed vcf files

Example.  To generate a new vcf.gz file from population1chr22.vcf.gz   
```
    python phase_singletons.py population1chr22.vcf.gz > population1chr2_singletonsphased.vcf.gz 
```    
	
## runtc.py	

runtc.py requires an input vcf file that has been phased for singletons (as described above).  
It returns to stdout a table of estimated tc values.   In the course of a run, msh files and a reversed vcf file are generated.
The most time consuming part of most analyses is the generation of msh files.  To speed this up, the program will use pypy3 if it has been installed and available
(https://pypy.org/download.html). 

The program can be used to simply generate msh files that are to be saved for futher analysis using the --msh-only option.  
A subsequent run would then use the --resuse option.

#### Command line usage
usage: runtc.py [-h] [--alledges] [--gzip] [--k1] [--k-all]  
                [--k-range K_RANGE K_RANGE] [--map MAPNAME] [--msh-only]  
                [--mut MUT_RATE] [--n0 N0_MODEL] [--outfn OUTFN]  
                [--outmsh OUTMSH] [--output-all-est] [--positions POSNAME]  
                [--rec REC_RATE] [--rev-vcf REVNAME] [--reuse]  
                [--randn RANDOM_N] [--seed RANDOM_N_SEED] [--sub SUBNAME]  
                vcfname  
  
##### positional arguments:  
* vcfname               vcf or vcf.gz input file  
  
##### optional arguments:  
*  -h, --help            show this help message and exit    
*  --alledges            If set, will run invariant site estimator on all external edges.  
*  --gzip                Compress msh files  
*  --k1                  Estimate tc values for singletons  
*  --k-all               Estimate tc values for all variants  
*  --k-range K_RANGE K_RANGE  Estimate tc values for inclusive range of k  
*  --map MAPNAME         Genetic map to be used for calculating genetic distances (either --rec or --map must be used)  
*  --msh-only            If set, will stop after generating files with msh values  
*  --mut MUT_RATE        Mutation rate per base per generation (default = 1e-8)  
*  --n0 N0_MODEL         Effective population size of sampled population (default: 10000)  
*  --outfn OUTFN         Output name for file of estimates  
*  --outmsh OUTMSH       Output label for msh files  
*  --output-all-est      Output all estimates instead of mean/max, use with --alledges or for non-singletons  
*  --positions POSNAME   List of base positions for which to generate estimates. Use with --alledges  
*  --rec REC_RATE        Recombination rate per base per generation (either --rec or --map must be used)  
*  --rev-vcf REVNAME     Path for reversed VCF file (by default this file will not be overwritten, even if --reuse is not invoked)  
*  --reuse               Will reuse available msh files or available reversed VCF file named *_reversed.vcf  
*  --randn RANDOM_N      use random_n chromosomes selected at random  
*  --seed RANDOM_N_SEED  seed value to use with --randn  
*  --sub SUBNAME         name of file with list of chromosome numbers to which analyses are limited (values are 0-based)  
				
#### Required Options

All runs require the following:

* vcf filename 
* --mut
* one of [--rec, --map]  
* one of [--alledges, --k1,  --k-range, --k-all]

If --msh-only is being used and the program is being used to generate msh files, then in addition to the above arguments, the following are required:
* --msh-only
* --outmsh  

if --reuse is being used to run on previously generated msh files then following are required
* --resuse 
* --outmsh  (followed by the same name used in the previous run using --msh-only) 

#### Examples
Examples using vcf file runtcexample.vcf which had been previously run through phase_singletons.py.  

To generate msh files for a run that will generate estimates at all SNPs with a genetic map in runtcexample.map and a mutation rate of 1e-8
```
python runtc.py  runtcexample.vcf --msh-only --mut 1e-8 --rec 1e-8 --k-all --outmsh example1 

```
To estimate tc values for all snps using previously saved msh files 
```
python runtc.py runtcexample.vcf --reuse --mut 1e-8 --rec 1e-8 --k-all --outmsh example1  > example2_kall.out 
```

Genetic map files have two columns of values, with the first column being ordered base positions (small to large) and the second column being 
the distance to that same point in centimorgans.  

To estimate tc values for singletons using a genetic map  
```
python runtc.py runtcexample.vcf --mut 1e-8 --map runtcexample.map --k1   > example3_k1.out 
```












