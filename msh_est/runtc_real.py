import os
import sys
#import msh_from_vcf
#import reversefile
#import aae_work
import argparse
from os.path import isfile
import subprocess
try:
    from msh_est import reverse_file,getmsh,run_estimator
except:
    from msh_from_vcf import getmsh
    from reversefile import reverse_file
    from aae_work import run_estimator

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcfname",type=str,help="VCF input file")
    parser.add_argument("--map",dest="mapname",type=str,help="Genetic map to be used for calculating genetic distances")
    parser.add_argument("--start",dest="start",type=int,default=-1,help="If set, start reporting estimates for SNP at this line")
    parser.add_argument("--end",dest="end",type=int,default=-1,help="If set, end reporting estimates for SNP at this line")
    parser.add_argument("--n",dest="n_model",type=int,default=100)
    parser.add_argument("--n0",dest="n0_model",type=int,default=1000,help="Population size of sampled population")
    parser.add_argument("--mut",dest="mut_rate",type=float,default=1e-8,help="Mutation rate for estimator")
    parser.add_argument("--rec",dest="rec_rate",type=float,default=1e-8,help="Recombination rate for estimator")
    parser.add_argument("--t0",dest="nosnp",action="store_true",help="If set, will run invariant site estimator")
    parser.add_argument("--mod-gen",dest="mod_gen",action="store_true",help="If set, will recalculate genetic distances in estimator using this formula: (1-e^2x)/2")
    parser.add_argument("--full-output",dest="full_out",action="store_true",help="Report lengths, chi, and estimate for all chromosomes in all sites")
    parser.add_argument("--region-mode",dest="region_mode",action="store_true",help="Report estimates for regions between SNPs instead of sites at and ignoring SNPs")
    parser.add_argument("--force",dest="force_override",action="store_true",help="Will override existing length files if they are present, default is to reuse if they exist")
    parser.add_argument("--nosquish",dest="squish",action="store_false",help="Read all lines from genetic map even if regions end up with 0 cM/bp rate")
    parser.add_argument("--outmsh",dest="outmsh",type=str,help="Output tag for length files")
    parser.add_argument("--outest",dest="outest",type=str,help="Output name for estimate file")
    parser.add_argument("--sub",dest="subname",type=str,help="file with list of chromosome numbers to use from vcf (0-based, e.g. individuals 0,3: 0,1,6,7")
    parser.add_argument("--rev",dest="revname",type=str,help="Path to reversed VCF (by default this file will never be overwritten by --force)")
    parser.add_argument("--nocache",dest="cache",action="store_false",help="Turn off caching of estimates in estimator")
    parser.add_argument("--bin",dest="bin",action="store_true",help="Calculate estimates by linear interpolation of geometrically-distributed pre-calculated estimates")
    parser.add_argument("--round",dest="round",type=int,default=-1,help="Round chi values to set number of significant digits for caching")
    parser.add_argument("--side-check",dest="side_check",action="store_true",help="Use only one-sided estimator if any chromosomes at a site are missing a side")
    return parser

def splitArgsForEstimator(args):
    arglist = []
    arglist.extend(['--start',str(args.start)])
    arglist.extend(['--end',str(args.end)])
    arglist.extend(['--n',str(args.n_model)])
    arglist.extend(['--n0',str(args.n0_model)])
    arglist.extend(['--mut',str(args.mut_rate)])
    arglist.extend(['--rec',str(args.rec_rate)])
    if args.nosnp:
        arglist.extend(['--t0'])
    if args.full_out:
        arglist.extend(['--full-output'])
    if args.mod_gen:
        arglist.extend(['--mod-gen'])
    if not args.cache:
        arglist.extend(['--nocache'])
    if args.bin:
        arglist.extend(['--bin'])
    if args.round != -1:
        arglist.extend(['--round',str(args.round)])
    if args.outest is not None:
        arglist.extend(['--outfn',str(args.outest)])
    if args.side_check:
        arglist.extend(['--side-check'])
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
    return msh_left_args,msh_right_args,leftmshfname,rightreversedmshfname,rightmshfname


def main(argv):
    parser = createParser()
    args = parser.parse_args(argv)
    vcfname = args.vcfname
    vcftag = vcfname[0:vcfname.rfind(".vcf")]
    if args.revname:
        rvcfname = args.revname
    else:
        rvcfname = vcftag+"_reversed.vcf"

    if (args.force_override or not isfile(rvcfname)) and args.revname is None:
        sys.stderr.write("Reversing vcf %s into %s\n" % (vcfname,rvcfname))
        temp_pypy_str = "\"from msh_est import reverse_file; reverse_file('%s','%s')\" "%(vcfname,rvcfname)
        temp_subprocess_return = subprocess.run("""pypy3 -c %s"""%(temp_pypy_str),shell=True)
        if temp_subprocess_return:
            reverse_file(vcfname,rvcfname)

    msh_left_args,msh_right_args,leftmshfname,rightreversedmshfname,rightmshfname = splitArgsForLengths(args,rvcfname)
    if args.force_override or not isfile(leftmshfname):
        sys.stderr.write("Creating left lengths: %s\n" %(str(msh_left_args)))
        temp_args_str = '['
        for a in msh_left_args:
            temp_args_str += "\'%s\',"%a
        temp_args_str = temp_args_str[:-1] + ']'  # replace ',' on end with ']'
        temp_pypy_str = "\"from msh_est import getmsh;  get_msh(%s)\""%temp_args_str
        temp_subprocess_return = subprocess.run("""pypy3 -c %s """%(temp_pypy_str),shell=True)
        if temp_subprocess_return:
            getmsh(msh_left_args)

    if args.force_override or not isfile(rightreversedmshfname):
        sys.stderr.write("Creating right lengths: %s\n" % (str(msh_right_args)))
        temp_args_str = "["
        for a in msh_right_args:
            temp_args_str += "\'%s\',"%a
        temp_args_str = temp_args_str[:-1] + ']' # replace ',' on end with ']'
        temp_pypy_str = "\"from msh_est import getmsh;  getmsh(%s)\""%temp_args_str
        temp_subprocess_return = subprocess.run("""pypy3 -c %s """%(temp_pypy_str),shell=True)
        if temp_subprocess_return:
            getmsh(msh_right_args)

    if args.force_override or not isfile(rightmshfname):
        sys.stderr.write("Reversing right lengths\n")
        temp_pypy_str = "\"from msh_est import reversefile;  reverse_file('%s','%s')\" "%(rightreversedmshfname,rightmshfname)
        temp_subprocess_return = subprocess.run("""pypy3 -c %s """%(temp_pypy_str),shell=True)
        if temp_subprocess_return:
            reverse_file(rightreversedmshfname,rightmshfname)

    est_args = [leftmshfname,rightmshfname] + splitArgsForEstimator(args)
    sys.stderr.write("Generating estimates: %s\n" % (str(est_args)))
    run_estimator(est_args)

if __name__ == "__main__":
    main(sys.argv[1:])
