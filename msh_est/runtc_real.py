import os
import sys
#import msh_from_vcf
#import reversefile
#import aae_work
import argparse
from os.path import isfile
import subprocess
import gzip
import datetime
#try:
#    from msh_est import reverse_file,getmsh,run_estimator
#except:
#    sys.stderr.write("msh_est not found, using local files\n")

from msh_from_vcf import getmsh
from reversefile import reverse_file
from aae_work import run_estimator

PYPY_VERSION="pypy3"



def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcfname",type=str,help="VCF input file")
    parser.add_argument("--map",dest="mapname",type=str,help="Genetic map to be used for calculating genetic distances")
    parser.add_argument("--n0",dest="n0_model",type=int,default=1000,help="Population size of sampled population")
    parser.add_argument("--mut",dest="mut_rate",type=float,default=1e-8,help="Mutation rate for estimator")
    parser.add_argument("--rec",dest="rec_rate",type=float,default=1e-8,help="Recombination rate for estimator")
    parser.add_argument("--alpha",dest="alpha",action="store_true",help="If set, will run invariant site estimator. If not, will run singleton mode estimator")
    parser.add_argument("--snp-mode",dest="snp_mode",action="store_true",help="Report lengths from adjacent variants (instead of region)")
    parser.add_argument("--reuse",dest="force_override",action="store_false",help="Will reuse any available length or reversed VCF file having run's tag")
    parser.add_argument("--nosquish",dest="squish",action="store_false",help="Read all lines from genetic map even if regions end up with 0 cM/bp rate")
    parser.add_argument("--outmsh",dest="outmsh",type=str,help="Output tag for length files")
    parser.add_argument("--outest",dest="outest",type=str,help="Output name for estimate file")
    parser.add_argument("--rev",dest="revname",type=str,help="Path to reversed VCF (by default this file will never be overwritten by --force)")
    parser.add_argument("--nocache",dest="cache",action="store_false",help="Turn off caching of estimates in estimator")
    parser.add_argument("--bin",dest="bin",action="store_true",help="Calculate estimates by linear interpolation of geometrically-distributed pre-calculated estimates")
    parser.add_argument("--round",dest="round",type=int,default=-1,help="Round chi values to set number of significant digits for caching")
    parser.add_argument("--gzip",dest="gzip_check",action="store_true",help="Will compress msh files and delete intermediate reversed length file")
    parser.add_argument("--positions",dest="posname",type=str,help="List of positions that should output regions for")
    parser.add_argument("--include-singletons",dest="inc_sing",action="store_true",help="Includes singletons when determining MSH cutoffs")
    parser.add_argument("--lengths-only",dest="lengths_only",action="store_true",help="If set, will stop after generating lengths")
    parser.add_argument("--pypy",dest="use_pypy",action="store_true",help="If set, use pypy to accelerate length generation")
    parser.add_argument("--k-all",dest="k_all",action="store_true",help=("Will provide outgroup MSH for all variants"))
    parser.add_argument("--k-range",dest="k_range",type=int,nargs=2,help=("Provide outgroup MSH for inclusive range of k"))
    parser.add_argument("--expmode",dest="expmode",action="store_true",help=("Use composite model in estimator"))
    parser.add_argument("--seed",dest="random_n_seed",type=int,default=-1,help="seed value to use with --randn")
    parser.add_argument("--output-all-est",dest="all_est",action="store_true",help=("Output all estimates instead of mean/max"))
    parser.add_argument("--dt-exp",dest="dt_exp",nargs=2,type=int)
    subgroup = parser.add_mutually_exclusive_group()
    subgroup.add_argument("--randn",dest="random_n",type=int,default=-1,help="use random_n chromosomes selected at random")
    subgroup.add_argument("--sub",dest="subname",type=str,help="file with list of chromosome numbers to use from vcf (0-based, e.g. individuals 0,3: 0,1,6,7")
    modelgroup = parser.add_mutually_exclusive_group()
    modelgroup.add_argument("--exp-model",dest="expmodel",type=float,help=("Use "
                            "exponential population model with given growth rate"))
    modelgroup.add_argument("--twophase-model",dest="twophase",nargs=2,type=float,
                            help=("Use two-phase (constant pop size then exp "
                            "growth) model, with growth rate as first argument "
                            "and growth time in generations as second"))
    return parser

def splitArgsForEstimator(args):
    arglist = []
    if args.random_n != -1:
        arglist.extend(['--n',str(args.random_n)])
    else:
        arglist.extend(['--n',str(getN(args))])
    arglist.extend(['--n0',str(args.n0_model)])
    arglist.extend(['--mut',str(args.mut_rate)])

    if args.alpha:
        arglist.extend(['--alpha'])
    if not args.cache:
        arglist.extend(['--nocache'])
    if args.bin:
        arglist.extend(['--bin'])
    if args.round != -1:
        arglist.extend(['--round',str(args.round)])
    if args.outest is not None:
        arglist.extend(['--outfn',str(args.outest)])
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
    if args.k_all or args.k_range is not None:
        arglist.append("--decmode")
    if args.expmode:
        arglist.append("--expmode")
    if args.expmodel is not None:
        arglist.extend(['--exp-model',str(args.expmodel)])
    if args.twophase is not None:
        arglist.extend(['--twophase-model',str(args.twophase[0]),str(args.twophase[1])])
    if args.all_est:
        arglist.append("--output-all-est")
    return arglist

def splitArgsForLengths(args,rvcfname):
    msh_left_args = ['--vcf',args.vcfname]
    msh_right_args = ['--vcf',rvcfname]
    vcftag = args.vcfname[0:args.vcfname.rfind(".vcf")]
    if args.outmsh is not None:
        mshtag = args.outmsh
    else:
        mshtag = vcftag
    leftmshfname = mshtag+"_left_lengths.txt"
    rightmshfname = mshtag+"_right_lengths.txt"
    if args.gzip_check:
        leftmshfname += ".gz"
        rightmshfname += ".gz"
    rightreversedmshfname = mshtag+"_reversed_right_lengths.txt"
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
    if not args.alpha and not args.k_all and args.k_range is None:
        msh_left_args.append('--singleton')
        msh_right_args.append('--singleton')
    if args.posname is not None:
        msh_left_args.extend(['--positions',str(args.posname)])
        msh_right_args.extend(['--positions',str(args.posname)])
        msh_right_args.extend(['--revpos'])
    if args.all_est:
        msh_left_args.append("--alpha-singleton-only")
        msh_right_args.append("--alpha-singleton-only")
    if args.inc_sing:
        msh_left_args.append('--include-singletons')
        msh_right_args.append('--include-singletons')
    if args.k_all:
        msh_left_args.append('--k-all')
        msh_right_args.append('--k-all')
    if args.k_range is not None:
        msh_left_args.extend(['--k-range',str(args.k_range[0]),str(args.k_range[1])])
        msh_right_args.extend(['--k-range',str(args.k_range[0]),str(args.k_range[1])])
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

def runWithPypy(pypy_version, script_name, args):
    exec_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),script_name)
    sub_string = pypy_version+' '+exec_path+' '+' '.join(map(str,args)) 
    subprocess_return = subprocess.run(sub_string,shell=True)
    if subprocess_return.returncode != 0:
        sys.stderr.write("Issue with pypy for %s"%(script_name)+'\n')
    return subprocess_return.returncode



def main(argv):
    parser = createParser()
    if argv[-1] =='':
        argv = argv[0:-1]
    args = parser.parse_args(argv)
    #print ("Start : "+datetime.datetime.now())
    vcfname = args.vcfname
    vcftag = vcfname[0:vcfname.rfind(".vcf")]
    if args.revname:
        rvcfname = args.revname
    else:
        rvcfname = vcftag+"_reversed.vcf"
        if vcfname[-6:] == 'vcf.gz':
            rvcfname += ".gz"
    usepypy = args.use_pypy

    if (args.force_override or not isfile(rvcfname)) and args.revname is None:
        sys.stderr.write("Reversing vcf %s into %s\n" % (vcfname,rvcfname))
        retcode = 1
        if usepypy:
            retcode = runWithPypy(PYPY_VERSION,'reversefile.py',[vcfname,rvcfname])
        if retcode != 0:
            reverse_file(vcfname,rvcfname)
    #print ("VCF reversed: "+datetime.datetime.now())
    msh_left_args,msh_right_args,leftmshfname,rightreversedmshfname,rightmshfname = splitArgsForLengths(args,rvcfname)
    if args.force_override or not isfile(leftmshfname):
        sys.stderr.write("Creating left lengths: %s\n" %(str(msh_left_args)))
        retcode = 1
        if usepypy:
            retcode = runWithPypy(PYPY_VERSION,'msh_from_vcf.py',msh_left_args)
        if retcode != 0:
            getmsh(msh_left_args)
    #print ("Left Lengths: "+datetime.datetime.now())
    if args.force_override or (not isfile(rightreversedmshfname) and not isfile(rightmshfname)):
        sys.stderr.write("Creating right lengths: %s\n" % (str(msh_right_args)))
        if usepypy:
            retcode = runWithPypy(PYPY_VERSION,'msh_from_vcf.py',msh_right_args)
        if retcode != 0:
            getmsh(msh_right_args)

    if args.force_override or not isfile(rightmshfname):
        sys.stderr.write("Reversing right lengths\n")
        retcode = 1
        if usepypy:
            retcode = runWithPypy(PYPY_VERSION,'reversefile.py',[rightreversedmshfname,rightmshfname])
        if retcode != 0:
            reverse_file(rightreversedmshfname,rightmshfname)
    #print ("Right lengths: "+datetime.datetime.now())
    if args.gzip_check and isfile(rightreversedmshfname):
        os.remove(rightreversedmshfname)
    if not args.lengths_only:
        est_args = [leftmshfname,rightmshfname] + splitArgsForEstimator(args)
        sys.stderr.write("Generating estimates: %s\n" % (str(est_args)))
        run_estimator(est_args)
    #print ("Done: "+datetime.datetime.now())

if __name__ == "__main__":
    main(sys.argv[1:])
