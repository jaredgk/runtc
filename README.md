
# RUNTC

Generate estimates of first coalescent time, tc  (Platt et al., 2019). 

This documentation covers the usage of *runtc* and *phase_singletons*.

## VCF Requirements

A vcf file should contain phased genotypes at all variable SNPs of a sample of aligned genomes for a single chromosome. 

If the file has not been phased at singleton sites,  then it should be first run through *phase_singletons*.

The *runtc* program ignores all sites with missing data or multiple derived alleles.


## phase_singletons

*phase_singletons* is a python program that takes as its main argument a vcf or vcf.gz file and returns to stdout a new file with new placements for singleton alleles as estimated from msh lengths. In the course of its operation it also generates two temporary msh (or msh.gz) files which can be saved if desired, but are deleted by default. 

#### Command line usage
##### positional arguments:
* the name of the vcf file (filename should end in '.vcf' or '.vcf.gz')

##### optional arguments:
*  -h, --help           show this help message and exit
*  --gzip               gzip temporary msh files
*  --keep               save temporary msh and reversed vcf files

Example.  To generate a new vcf.gz file from population1chr22.vcf.gz   
```
    python phase_singletons.py population1chr22.vcf.gz > population1chr2_singletonsphased.vcf.gz 
```    
	
## runtc	

*runtc* is a python program that takes as input a vcf file that has been phased for singletons (as described above).  It returns to stdout a table of estimated tc values.   In the course of a run, msh files and a reversed vcf file are generated.
The most time consuming part of most analyses is the generation of msh files.  To speed this up, the user can compile the C++ program *msh_vcf* and have *runtc* use this for generating msh files.  If *msh_vcf* is not used, then the program will use pypy3 for making msh files, if it has been installed and available (https://pypy.org). 

If desired, the program can be run so as to only generate msh files using the --msh-only option.  A subsequent run would then use the --resuse option.

### Command line usage
usage: runtc.py [-h] [--alledges | --k1 | --k-all | --k-range K_RANGE K_RANGE]
                --map MAPNAME | --rec REC_RATE] [--gzip] [--msh-only]
                [--mut MUT_RATE] [--n0 N0_MODEL] [--outfn OUTFN]
                [--outmsh OUTMSH] [--output-all-est] [--positions POSNAME]
                [--rev-vcf REVNAME] [--reuse] [--randn RANDOM_N]
                [--seed RANDOM_N_SEED] [--sub SUBNAME] [--c-msh]
                vcfname

  
##### positional arguments:  
* the name of the vcf file (filename should end in '.vcf' or '.vcf.gz')
  
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
*    --c-msh               Use C++ version of MSH code, must be compiled in
                        script directory

				
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

* --reuse 
* --outmsh  (followed by the same name used in the previous run using --msh-only) 

### Overview: generating tc estimates with *runtc*

A full run of *runtc*, starting with only a vcf file (or vcf.gz file), will carry out the following steps:

 * make a copy of the vcf file with the order of lines reversed
 * make the left msh file (\*\_left_msh.txt) from the vcf file   
 * make the reversed-right msh file (\*\_reversed_right_msh.txt) from the reversed vcf file  
 * reverse the reversed-right msh file to make the right msh file
 * generate tc estimates from the msh files


#### Working with large vcf files

It can be useful,  particularly when working with large vcf files  (e.g. 100's of Gigabytes or larger), to do the vcf reversal and the making of the msh files as separate jobs.  For example,  to reverse a file without doing any other operations, use reversefile.py.  To generate msh files as a standalone job invoke --msh-only and --outmsh when running *runtc*, or use the compiled executable *msh_vcf* directly.  Then the tc calculations can be done by running *runtc* while invoking --reuse and --outmsh  to use the previously generated reversed vcf file and previously generated msh files.  All of these jobs should be run in the same folder as the vcf file. 

Generating msh files can often be the most time consuming step.  Use of the C++ program *msh_vcf*  can speed this up considerably.  This can be run in standalone mode,  or if a compiled executable is in the same folder as the python scripts,  then invoking --c-msh at the command line for *runtc*  will use this program. 

Users should be aware of available disk space, as the reversed vcf file will be as large as the original and the msh files can be larger than the vcf files.  The *runtc* program can start with a vcf.gz file,  however all file-reversal and msh-related operations required uncompressed files, so uncompressed files will be generated in the course of the run.  

### Examples
Examples using vcf file runtcexample.vcf which had been previously run through *phase_singletons*.  

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

To estimate tc values for alleles that occur 2 thru 10 times using a genetic map and using the msh_vcf executable to make msh files  

```
python runtc.py runtcexample.vcf --mut 1e-8 --map runtcexample.map --k-range 2 10 --c-msh > example4_k1.out
```

## References
Platt A, Pivirotto A, Knoblauch J, Hey J. 2019. An estimator of first coalescent time reveals selection on young variants and large heterogeneity in rare allele ages among human populations. PLoS Genetics 15:e1008340.


