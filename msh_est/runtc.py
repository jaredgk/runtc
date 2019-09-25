#-------------------------------------------------------------------------------
# Name:        runtc.py
# Authors:     Jared Knoblauch, Alex Platt, Jody Hey
# Created:     11/07/2019
# Copyright:   (c) Jody Hey 2019
#-------------------------------------------------------------------------------

import os
import sys
import argparse
from os.path import isfile
import subprocess
import gzip
import datetime
from msh_from_vcf import getmsh
from reversefile import reverse_file
from aae_work import run_estimator

PYPY_VERSION="pypy3"

def createParser():
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    # exclusive groups
    rec_group = parser.add_mutually_exclusive_group()    # --rec and --map
    k_edges_group = parser.add_mutually_exclusive_group() # --alledges,  --k1,  --k-all and --k-range
    subgroup = parser.add_mutually_exclusive_group()  # --randn and --sub
    modelgroup = parser.add_mutually_exclusive_group()  #   --exp-model and --twophase-model

    parser.add_argument("vcfname",type=str,help="vcf or vcf.gz input file")
    k_edges_group.add_argument("--alledges",dest="alledges",action="store_true",help="If set, will run invariant site estimator on all external edges.")
    k_edges_group.add_argument("--k1",dest="k1",action="store_true",help=("Estimate tc values for singletons"))  ## k1 for singletons is not really necessary as it is the default run mode, but use it for clarity
    k_edges_group.add_argument("--k-all",dest="k_all",action="store_true",help=("Estimate tc values for all variants"))
    k_edges_group.add_argument("--k-range",dest="k_range",type=int,nargs=2,help=("Estimate tc values for inclusive range of k"))
    k_edges_group.add_argument("--k-list",dest="k_list",type=int,nargs="+",help=("Estimate tc values for listed k values"))
    rec_group.add_argument("--map",dest="mapname",type=str,help="Genetic map to be used for calculating genetic distances (either --rec or --map must be used)")
    rec_group.add_argument("--rec",dest="rec_rate",type=float,help="Recombination rate per base per generation (either --rec or --map must be used)")
    parser.add_argument("--gzip",dest="gzip_check",action="store_true",help="Compress msh files")
    parser.add_argument("--msh-only",dest="msh_only",action="store_true",help="If set, will stop after generating files with msh values")
    parser.add_argument("--mut",dest="mut_rate",type=float,help="Mutation rate per base per generation (default = 1e-8)")
    parser.add_argument("--n0",dest="n0_model",type=int,default=10000,help="Effective population size of sampled population (default: 10000)")
    parser.add_argument("--outfn",dest="outfn",type=str,help="Output name for file of estimates")
    parser.add_argument("--outmsh",dest="outmsh",type=str,help="Output label for msh files")
    parser.add_argument("--output-all-est",dest="all_est",action="store_true",help=("Output all estimates instead of mean/max, use with --alledges or for non-singletons"))
    parser.add_argument("--positions",dest="posname",type=str,help="List of base positions for which to generate estimates. Use with --alledges")
    parser.add_argument("--rev-vcf",dest="revname",type=str,help="Path for reversed VCF file (by default this file will not be overwritten, even if --reuse is not invoked)")
    parser.add_argument("--reuse",dest="force_overwrite",action="store_false",help="Will reuse available msh files or available reversed VCF file named *_reversed.vcf")
    subgroup.add_argument("--randn",dest="random_n",type=int,default=-1,help="use random_n chromosomes selected at random")
    parser.add_argument("--seed",dest="random_n_seed",type=int,default=-1,help="seed value to use with --randn")
    parser.add_argument("--c-msh",dest="c_msh",action="store_true",help="Use C++ version of MSH code, must be compiled in script directory")  
    parser.add_argument("--keep-msh-files",action="store_true",help=argparse.SUPPRESS)
    parser.add_argument("--est-seed",dest="est_seed",type=int,help=argparse.SUPPRESS)
    subgroup.add_argument("--sub",dest="subname",type=str,help="name of file with list of chromosome numbers to which analyses are limited (values are 0-based)")
    # args with help suppressed as of 7/12/2019
 
    modelgroup.add_argument("--exp-model",dest="expmodel",type=float,help=argparse.SUPPRESS) ## suppressed 7/11/2019 JH  original help=("Use exponential population model with given growth rate"))
 
    modelgroup.add_argument("--twophase-model",dest="twophase",nargs=2,type=float,help=argparse.SUPPRESS) ## suppressed 7/11/2019 JH  original  help=("Use two-phase (constant pop size then exp growth) model, with growth rate as first argument and growth time in generations as second"))
    parser.add_argument("--dt-exp",dest="dt_exp",nargs=2,type=int,help=argparse.SUPPRESS) ## suppressed 7/11/2019 JH  original
    parser.add_argument("--nopypy",dest="no_pypy",action="store_true",default = False,help=argparse.SUPPRESS) ## suppressed and changed 7/11/2019 JH changed from  parser.add_argument("--pypy",dest="use_pypy",action="store_true",help="If set, use pypy to accelerate length generation")
    parser.add_argument("--exclude-singletons",dest="exc_sing",action="store_false",help=argparse.SUPPRESS) ## changed from --include-singletons and suppressed 7/11/2019 JH  original help="Includes singletons when determining MSH cutoffs")


    parser.add_argument("--round",dest="round",type=int,default=3,help=argparse.SUPPRESS) ## suppressed 7/11/2019 JH  original help="Round chi values to set number of significant digits for caching, -1 for no rounding")
    parser.add_argument("--snp-mode",dest="snp_mode",action="store_true",help=argparse.SUPPRESS) ## suppressed 7/11/2019 JH  original help="Report msh values from adjacent variants (instead of region)")
    parser.add_argument("--nosquish",dest="squish",action="store_false",help=argparse.SUPPRESS) ## suppressed 7/11/2019 JH  original help="Read all lines from genetic map even if regions end up with 0 cM/bp rate")
    parser.add_argument("--nocache",dest="cache",action="store_false",help=argparse.SUPPRESS) ## suppressed 7/11/2019 JH  original help="Turn off caching of estimates in estimator")
    parser.add_argument("--bin",dest="bin",action="store_true",help=argparse.SUPPRESS) ## suppressed 7/11/2019 JH  original help="Calculate estimates by linear interpolation of geometrically-distributed pre-calculated estimates")
    #parser.add_argument("--c-msh",dest="c_msh",action="store_true",help="Use C++ version of MSH code, must be compiled in script directory")


    return parser

def splitArgsForEstimator(args):
    arglist = []
    if args.random_n != -1:
        arglist.extend(['--n',str(args.random_n)])
    else:
        arglist.extend(['--n',str(getN(args))])
    arglist.extend(['--n0',str(args.n0_model)])
    arglist.extend(['--mut',str(args.mut_rate)])

    if args.alledges:
        arglist.extend(['--alledges'])
    if not args.cache:
        arglist.extend(['--nocache'])
    if args.bin:
        arglist.extend(['--bin'])
    if args.round != -1:
        arglist.extend(['--round',str(args.round)])
    if args.outfn is not None:
        arglist.extend(['--outfn',str(args.outfn)])
    if args.snp_mode:
        arglist.extend(['--snp-mode'])
    if args.mapname is not None:
        arglist.extend(['--gen',args.mapname])
    else:
        arglist.extend(['--rec',str(args.rec_rate)])
    if not args.squish:
        arglist.append('--nosquish')
    if args.posname is not None:
        arglist.append('--pos')
    if args.k_all or args.k_range is not None or args.k1 or args.k_list is not None:
        arglist.append("--kmode")
    if args.expmodel is not None:
        arglist.extend(['--exp-model',str(args.expmodel)])
    if args.twophase is not None:
        arglist.extend(['--twophase-model',str(args.twophase[0]),str(args.twophase[1])])
    if args.all_est:
        arglist.append("--output-all-est")
    if args.est_seed is not None:
        arglist.extend(['--seed',str(args.est_seed)])
    return arglist

def splitArgsForLengths(args,rvcfname):
    msh_left_args = ['--vcf',args.vcfname]
    msh_right_args = ['--vcf',rvcfname]
    vcftag = args.vcfname[0:args.vcfname.rfind(".vcf")]
    if args.outmsh is not None:
        mshtag = args.outmsh
    else:
        mshtag = vcftag
    leftmshfname = mshtag+"_left_msh.txt"
    rightmshfname = mshtag+"_right_msh.txt"
    if args.gzip_check:
        leftmshfname += ".gz"
        rightmshfname += ".gz"
    rightreversedmshfname = mshtag+"_reversed_right_msh.txt"
    msh_left_args.extend(['--out',leftmshfname])
    msh_right_args.extend(['--out',rightreversedmshfname])
    if args.mapname is not None:
        msh_left_args.extend(['--gen',args.mapname])
        msh_right_args.extend(['--gen',args.mapname])
    if args.subname is not None:
        msh_left_args.extend(['--sub',args.subname])
        msh_right_args.extend(['--sub',args.subname])
    if not args.squish:
        msh_left_args.append('--nosquish')
        msh_right_args.append('--nosquish')
    if args.round != -1:
        msh_left_args.extend(['--round',str(args.round)])
        msh_right_args.extend(['--round',str(args.round)])
    #if not args.alledges and not args.k_all and args.k_range is None:
    #    msh_left_args.append('--singleton')  ##Singleton mode: Generate two mshs, one for each haplotype of an individual with a singleton
    #    msh_right_args.append('--singleton') ##Singleton mode: Generate two mshs, one for each haplotype of an individual with a singleton
    if args.posname is not None:
        msh_left_args.extend(['--positions',str(args.posname)])
        msh_right_args.extend(['--positions',str(args.posname)])
        msh_right_args.extend(['--revpos'])
    if args.all_est:
        msh_left_args.append("--alledges-singleton-only")   ##Output mshs for alledges estimator at singleton sites only
        msh_right_args.append("--alledges-singleton-only")  ##Output mshs for alledges estimator at singleton sites only
    if args.exc_sing:
        msh_left_args.append('--exclude-singletons')   ##DO NOT Allow singletons to terminate MSH lengths
        msh_right_args.append('--exclude-singletons')  ##DO NOT Allow singletons to terminate MSH lengths
    if args.k_all:
        msh_left_args.append('--k-all')
        msh_right_args.append('--k-all')
    if args.k_range is not None:
        msh_left_args.extend(['--k-range',str(args.k_range[0]),str(args.k_range[1])])
        msh_right_args.extend(['--k-range',str(args.k_range[0]),str(args.k_range[1])])
    if args.k_list is not None:
        msh_left_args.extend(['--k-list']+[str(i) for i in args.k_list])
        msh_right_args.extend(['--k-list']+[str(i) for i in args.k_list])
    if args.k1:
        msh_left_args.append("--k1")
        msh_right_args.append("--k1")      
    if args.dt_exp is not None:
        msh_left_args.extend(['--dt-exp',str(args.dt_exp[0]),str(args.dt_exp[1])])
        msh_right_args.extend(['--dt-exp',str(args.dt_exp[0]),str(args.dt_exp[1])])
    if args.random_n != -1:
        totaln = getN(args,skipsub = True)
        if args.random_n > totaln:
            raise Exception("Random subsample of %d is larger than sample size %d"%(args.random_n,totaln))
        #msh_left_args.extend(['--n0',str(totaln)])
        #msh_right_args.extend(['--n0',str(totaln)])
        import random
        if args.random_n_seed == -1:
            seed = random.randint(1,1000000)
        else:
            seed = args.random_n_seed
        msh_left_args.extend(['--randn-seed',str(seed),str(args.random_n)])
        msh_right_args.extend(['--randn-seed',str(seed),str(args.random_n)])
    return msh_left_args,msh_right_args,leftmshfname,rightreversedmshfname,rightmshfname

def getN(args, skipsub = False):
    if args.subname is not None and skipsub != True:
        sf = open(args.subname,'r')
        sl = sf.readlines()
        return len(sl)
    else:
        if args.vcfname[-3:] == '.gz':
            tf = gzip.open(args.vcfname,'rt')
        else:
            tf = open(args.vcfname,'r')
        l = '#'
        while l[0] == '#':
            l = tf.readline()
        la = l.strip().split()
        n = 0
        for i in range(9,len(la)):
            laa = la[i].split(':')[0]
            sc = len(laa.replace('|','/').split('/'))
            n += sc
        return n
        
def decompressFile(fn,tempfn):
    f = gzip.open(fn)
    tempf = open(tempfn,'w')
    for line in f:
        tempf.write(line.decode())
        

def runWithPypy(pypy_version, script_name, args):
    exec_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),script_name)
    sub_string = pypy_version+' '+exec_path+' '+' '.join(map(str,args))
    subprocess_return = subprocess.run(sub_string,shell=True,stderr=subprocess.PIPE)
    if subprocess_return.returncode != 0:
        sys.stderr.write("Issue with pypy for %s"%(script_name)+'\n')
        if 'multiple chromosomes' in subprocess_return.stderr:
            sys.stderr.write(subprocess_return.stderr+'\n')
            exit()
    return subprocess_return.returncode

def runCMsh(args):
    exec_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'msh_vcf')
    sub_string = exec_path+' '+' '.join(map(str,args))
    if not isfile(exec_path):
        sys.stderr.write("C MSH code not compiled, run 'Make' in msh-python directory\n")
        return 1
    try:
        subprocess_return = subprocess.run(sub_string,shell=True,stderr=subprocess.PIPE)
    except:
        sys.stderr.write("Exception on subprocess C MSH call\n")
        return 1
    if subprocess_return.returncode != 0:
        sys.stderr.write("Issue with C msh_vcf\n")
        sys.stderr.write(str(subprocess_return.stderr)+'\n')
    return subprocess_return.returncode


def makemshfiles(args):
    #print ("Start : "+datetime.datetime.now())
    vcfname = args.vcfname
    vcftag = vcfname[0:vcfname.rfind(".vcf")]
    invcf_compressed = False
    if args.revname:
        rvcfname = args.revname
    else:
        rvcfname = vcftag+"_reversed.vcf"
        if vcfname[-6:] == 'vcf.gz':
            rvcfname += ".gz"
            invcf_compressed = True
    usepypy = args.no_pypy
    msh_left_args,msh_right_args,leftmshfname,rightreversedmshfname,rightmshfname = splitArgsForLengths(args,rvcfname)
    if (args.force_overwrite or not isfile(rvcfname)) and (args.revname is None or not isfile(args.revname) and not (not args.force_overwrite and isfile(rightmshfname))):
        sys.stderr.write("Reversing vcf %s into %s\n" % (vcfname,rvcfname))
        retcode = 1
        if invcf_compressed:
            trevvcf_name = vcftag+'_temp.vcf'
            decompressFile(vcfname,trevvcf_name)
        else:
            trevvcf_name = vcfname
        if usepypy:
            retcode = runWithPypy(PYPY_VERSION,'reversefile.py',[trevvcf_name,rvcfname])
        if retcode != 0:
            reverse_file(trevvcf_name,rvcfname)
        if invcf_compressed:
            os.remove(trevvcf_name)
            
    #exit()
    #print ("VCF reversed: "+datetime.datetime.now())
    #msh_left_args,msh_right_args,leftmshfname,rightreversedmshfname,rightmshfname = splitArgsForLengths(args,rvcfname)
    if args.force_overwrite or not isfile(leftmshfname):
        sys.stderr.write("Creating left msh values: %s\n" %(str(msh_left_args)))
        retcode = 1
        if args.c_msh:
            retcode = runCMsh(msh_left_args)
        elif usepypy:
            retcode = runWithPypy(PYPY_VERSION,'msh_from_vcf.py',msh_left_args)
        if retcode != 0:
            getmsh(msh_left_args)
    #print ("Left msh values: "+datetime.datetime.now())
    if args.force_overwrite or (not isfile(rightreversedmshfname) and not isfile(rightmshfname)):
        sys.stderr.write("Creating right msh values: %s\n" % (str(msh_right_args)))
        if args.c_msh:
            retcode = runCMsh(msh_right_args)
        elif usepypy:
            retcode = runWithPypy(PYPY_VERSION,'msh_from_vcf.py',msh_right_args)
        if retcode != 0:
            getmsh(msh_right_args)

    if args.force_overwrite or not isfile(rightmshfname):
        sys.stderr.write("Reversing right msh values\n")
        retcode = 1
        if usepypy:
            retcode = runWithPypy(PYPY_VERSION,'reversefile.py',[rightreversedmshfname,rightmshfname])
        if retcode != 0:
            reverse_file(rightreversedmshfname,rightmshfname)
    #print ("Right msh values: "+datetime.datetime.now())
    if isfile(rightreversedmshfname):
        os.remove(rightreversedmshfname)
    est_args = [leftmshfname,rightmshfname] + splitArgsForEstimator(args)
    return est_args

def main(argv):
    parser = createParser()
    if argv[-1] =='':
        argv = argv[0:-1]
    args = parser.parse_args(argv)
    #print (args.k_list)
    #exit()

##    k_edges_group
    if not (args.alledges or args.k1 or args.k_all or args.k_range or args.k_list):
        parser.error('one of [--alledges, --k1,  --k-range, --k-all, --k-list]  must be invoked')
    #if not (args.alledges or args.k_all or args.k_range):
    #    args.k1 = True
##  rec_group options
    if not (args.rec_rate or args.mapname):
        parser.error('one of [--rec, --map]  must be invoked')
    if args.msh_only and not args.outmsh:
        parser.error("--outmsh must be used with --msh-only")
    if args.k_range is not None and args.k_range[0] < 2:
        parser.error("--k-range requires a range start 2 or greater")
    if args.k_range is not None and args.k_range[0] > args.k_range[1]:
        parser.error("Start of --k-range must be not be greater than end")


    est_args = makemshfiles(args)
    if not args.msh_only:
        if args.k1:
            tempstr = " tc values for singletons only "
        elif args.k_all:
            tempstr = " tc values for all SNPs"
        elif args.k_range:
            tempstr = "tc values in SNPs with counts in range {}-{} ".format(args.k_range[0],args.k_range[1])
        else:
            assert args.alledges
            tempstr = "t0 values for all edges at positions"
            if args.posname:
                tempstr += " given in file {}".format(args.posname)
            else:
                tempstr += " adjacent to SNPs"
        sys.stderr.write("Generating estimates (%s): %s\n" % (tempstr,str(est_args)))
        estat = run_estimator(est_args)
        if not args.keep_msh_files:
            if os.path.exists(est_args[0]):
                os.remove(est_args[0])
            if os.path.exists(est_args[1]):
                os.remove(est_args[1])

if __name__ == "__main__":
    if len(sys.argv) < 2:
        #sys.stderr.write("No arguments. Use -h  or --help for help menu")
        main(['-h'])
    else:
        main(sys.argv[1:])
