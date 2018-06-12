import os
import sys
import msh_from_vcf
import reversefile
import aae_work
import argparse
from os.path import isfile
import subprocess

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcfname")
    parser.add_argument("--map",dest="mapname",type=str)
    parser.add_argument("--start",dest="start",type=int,default=-1)
    parser.add_argument("--end",dest="end",type=int,default=-1)
    parser.add_argument("--n",dest="n_model",type=int,default=100)
    parser.add_argument("--n0",dest="n0_model",type=int,default=1000)
    parser.add_argument("--mut",dest="mut_rate",type=float,default=1e-8)
    parser.add_argument("--rec",dest="rec_rate",type=float,default=1e-8)
    parser.add_argument("--t0",dest="nosnp",action="store_true")
    parser.add_argument("--no-cache",dest="cache_est",action="store_false")
    parser.add_argument("--mod-gen",dest="mod_gen",action="store_true")
    parser.add_argument("--full-output",dest="full_out",action="store_true")
    parser.add_argument("--region-mode",dest="region_mode",action="store_true")
    parser.add_argument("--force",dest="force_override",action="store_true")
    parser.add_argument("--nosquish",dest="squish",action="store_false")
    parser.add_argument("--outtag",dest="outtag",type=str)
    parser.add_argument("--sub",dest="subname",type=str,help="file with list of chromosome numbers to use from vcf (0-based, e.g. individuals 0,3: 0,1,6,7")
    parser.add_argument("--rev",dest="revname",type=str)
    parser.add_argument("--nocache",dest="cache",action="store_false")
    parser.add_argument("--bin",dest="bin",action="store_true")
    parser.add_argument("--round",dest="round",type=int,default=-1)
    parser.add_argument("--side-check",dest="side_check",action="store_true")
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
    if args.outtag is not None:
        arglist.extend(['--outfn',str(args.outtag)+"_estimates.txt"])
    if args.side_check:
        arglist.extend(['--side-check'])
    return arglist


def main(argv):
    parser = createParser()
    args = parser.parse_args(argv)
    vcfname = args.vcfname
    if args.revname is not None:
        rvcfname = args.revname
    else:
        rvcfname = vcfname[0:vcfname.rfind(".vcf")] + "_reversed.vcf"
    #rvcfname = '/home/jared/workspace/alex/bwt/chr22_reversed.vcf.gz'
    if (args.force_override or not isfile(rvcfname)) and args.revname is  None:
        try:
            temp_pypy_str = "\"import reversefile; reversefile.reverse('%s','%s')\" "%(vcfname,rvcfname)
            temp_subprocess_return = subprocess.run("""pypy3 -c %s"""%(temp_pypy_str),shell=True)
            if temp_subprocess_return:
                reversefile.reverse(vcfname,rvcfname)
        except Exception:
            reversefile.reverse(vcfname,rvcfname)
    if args.outtag is None:
        outtag = vcfname[0:vcfname.rfind(".vcf")]
    else:
        outtag = args.outtag
    msh_left_args = ['--vcf',vcfname]
    msh_right_args = ['--vcf',rvcfname]
    if args.mapname is not None:
        msh_left_args.extend(['--gen',args.mapname])
        msh_right_args.extend(['--gen',args.mapname])
    if args.subname is not None:
        msh_left_args.extend(['--sub',args.subname])
        msh_right_args.extend(['--sub',args.subname])
    if not args.squish:
        msh_left_args.append('--nosquish')
        msh_right_args.append('--nosquish')
    #msh_left_args.extend(['--dir','left'])
    #leftmshfname = vcfname[0:vcfname.rfind(".vcf")]+"_left_results.txt"
    leftmshfname = outtag+"_left_results.txt"
    msh_left_args.extend(['--out',leftmshfname])
    if args.force_override or not isfile(leftmshfname):
        try:
            temp_args_str = '['
            for a in msh_left_args:
                temp_args_str += "\'%s\',"%a
            temp_args_str = temp_args_str[:-1] + ']'  # replace ',' on end with ']'
            temp_pypy_str = "\"import msh_from_vcf;  msh_from_vcf.getmsh(%s)\""%temp_args_str
            temp_subprocess_return = subprocess.run("""pypy3 -c %s """%(temp_pypy_str),shell=True)
            if temp_subprocess_return:
                msh_from_vcf.getmsh(msh_left_args)
        except:
            msh_from_vcf.getmsh(msh_left_args)

    #msh_right_args.extend(['--dir','right'])
    #rightreversedmshfname = rvcfname[0:rvcfname.rfind(".vcf")]+"_right_results.txt"
    rightreversedmshfname = outtag+'_reversed_right_results.txt'
    msh_right_args.extend(['--out',rightreversedmshfname])
    if args.force_override or not isfile(rightreversedmshfname):
        try:
            temp_args_str = "["
            for a in msh_right_args:
                temp_args_str += "\'%s\',"%a
            temp_args_str = temp_args_str[:-1] + ']' # replace ',' on end with ']'
            temp_pypy_str = "\"import msh_from_vcf;  msh_from_vcf.getmsh(%s)\""%temp_args_str
            temp_subprocess_return = subprocess.run("""pypy3 -c %s """%(temp_pypy_str),shell=True)
            if temp_subprocess_return:
                msh_from_vcf.getmsh(msh_right_args)
        except:
            msh_from_vcf.getmsh(msh_right_args)

    #rightmshfname = rightreversedmshfname[0:rightreversedmshfname.find("_reversed")]  + rightreversedmshfname[rightreversedmshfname.find("_right_results"):]
    rightmshfname = outtag+"_right_results.txt"
    if args.force_override or not isfile(rightmshfname):
        try:
            temp_pypy_str = "\"import reversefile;  reversefile.reverse('%s','%s')\" "%(rightreversedmshfname,rightmshfname)
            temp_subprocess_return = subprocess.run("""pypy3 -c %s """%(temp_pypy_str),shell=True)
            if temp_subprocess_return:
                reversefile.reverse(rightreversedmshfname,rightmshfname)
        except:
            reversefile.reverse(rightreversedmshfname,rightmshfname)
    est_args = [leftmshfname,rightmshfname] + splitArgsForEstimator(args)
    aae_work.run(est_args)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        printf("Usage: runtc.py (VCF name) (optional name for genetic map file)")
