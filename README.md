# msh-python

Compute tc/alpha values for SNPs in a given VCF file. There are three components to this process: reversing the original VCF 
file,generating the lengths of the most shared haplotypes using the Burrows-Wheeler transform, and using those lengths to 
estimate tc. As of now, there is no installer for this package, so usage is simply: 
```
python path/to/files/runtc_real.py [options]
```
If --pypy is to be used for faster length calculation, the scripts must all be in the directory from which the runtc_real 
script is called.

The only required argument is a VCF file, either uncompressed or gzipped. Using the default options, the runtc_real script will 
generate tc estimates for every singleton in a VCF file, while not counting singletons toward the termination of MSH tracts. 
Both haplotypes from an individual with the singleton at a site will have tc estimated, with the larger estimate returned as
the chosen one. The output file has two columns: physical position and estimate. 

### Model Options
--n0: sample population size (default 1000)

--mut: Mutation rate per base pair (default 1e-8)

--rec: Recombination rate per base pair (default 1e-8)

### Analysis Options

--alpha: Instead of only generating tc estimates for singletons, will identify MSH for every haplotype at all sites and 
output the geometric mean of the tc estimates generated per site.

--include-singletons: Allow singleton sites to terminate MSH lengths.

--k-range: For variants with allele counts within the specified range, find every derived allele's longest MSH shared with a haplotype with the reference allele, then use composite estimator to generate tc estimate for that site.

--k-all: Same as --k-range, but will do for all sites.

### Input options

--rev: Filename for reversed VCF file, if reversal is done outside of pipeline. This will prevent automatic generation of 
reversed VCF file.

--map: Filename for genetic map file, with column 1 having physical positions and column 2 having corresponding genetic 
positions.

--reuse: If length/reversed VCF files are present for given output prefix, will not attempt to regenerate them.

### Output options:

--outmsh: Prefix for MSH length files (defaults to name of input VCF stripped of .vcf extension)

--outest: Name for output estimates file (default is stdout)

--gzip: Will compress length files.

--lengths-only: Will stop pipeline once lengths are generated.


