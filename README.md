# msh-python

Compute tc/alpha values for SNPs in a given VCF file. There are three components to this process: reversing the original VCF 
file,generating the lengths of the most shared haplotypes using the Burrows-Wheeler transform, and using those lengths to 
estimate tc. As of now, there is no installer for this package, so usage is simply: 
```
python path/to/scripts/runtc_real.py [options]
```

The only required argument is a VCF file, either uncompressed or gzipped. Using the default options, the runtc_real script will filter out any non-biallelic sites, indels, and sites with missing data, then
generate tc estimates for every singleton in a VCF file, while not counting singletons toward the termination of MSH tracts. 
Both haplotypes from an individual with the singleton at a site will have tc estimated, with the larger estimate returned as
the chosen one. The output file has two columns: physical position and estimate. 

An example estimates file, at 'test/test_alpha_estimates.txt', contains estimates for alpha values at each SNP in the
provided VCF file, 'test/test_head.vcf'. The commandline used to generate a file for comparison:
```
python msh_est/runtc_real.py test/test_head.vcf --alpha --mut 1e-9 --outest test_estimates.txt
```
This will reverse the VCF file (creating a 'test/test_head_reversed.vcf' file), generate MSH lengths for every haplotype
at every SNP in the VCF file (--alpha), and generate an alpha estimate for each SNP using a chosen mutation rate. A diff 
between this file and the test/test_alpha_estimates.txt files will likely show differences, but the estimates themselves should be similar.

### Model Options
--n0 [pop_size]: sample population size (default 1000)

--mut [mutation_rate] : Mutation rate per base pair (default 1e-8)

--rec [recomb_rate]: Recombination rate per base pair (default 1e-8)

--exp-model [growth rate]: Use estimator with exponential population growth model, providing value for growth rate
per generation

--twophase-model [growth_rate] [growth_start]: Use estimator with a constant population size followed by an exponential 
phase growing at (growth_rate) starting (growth_start) generations ago.

### Analysis Options

--alpha: Instead of only generating tc estimates for singletons, will identify MSH for every haplotype at all sites and 
output the geometric mean of the tc estimates generated per site.

--include-singletons: Allow singleton sites to terminate MSH lengths.

--k-range [low_bound] [high_bound]: For variants with allele counts within the specified range, find every derived allele's longest MSH shared with a haplotype with the reference allele, then use composite estimator to generate tc estimate for that site.

--k-all: Same as --k-range, but will do for all sites.

### Input options

--pypy: Use Pypy for reversing files and calculating lengths. Defaults to pypy3. If any errors are encountered, will run in standard python mode. 

--rev [reversed_vcfname]: Filename for reversed VCF file, if reversal is done outside of pipeline. This will prevent automatic generation of 
reversed VCF file.

--map [map_filename]: Filename for genetic map file, with column 1 having physical positions and column 2 having corresponding genetic 
positions.

--reuse: If length/reversed VCF files are present for given output prefix, will not attempt to regenerate them.

### Output options:

--outmsh [out_prefix]: Prefix for MSH length files (defaults to name of input VCF stripped of .vcf extension).

--outest [est_filename]: Name for output estimates file (default is stdout).

--gzip: Will compress length files.

--lengths-only: Will stop pipeline once lengths are generated.


