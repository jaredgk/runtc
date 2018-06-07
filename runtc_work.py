import sys
import msh_from_vcf
import reversefile
import aae_work
import argparse
from os.path import isfile

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
    parser.add_argument("--squish",dest="squish",action="store_true")
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
    return arglist


def main(argv):
    parser = createParser()
    args = parser.parse_args(argv)
    vcfname = args.vcfname
    rvcfname = vcfname[0:vcfname.rfind(".vcf")] + "_reversed.vcf"
    if args.force_override or not isfile(rvcfname):
        reversefile.reverse(argv[0],rvcfname)
    else:
        sys.stdout.write(str(args.force_override)+'\t'+str(isfile(rvcfname))+'\n')

    msh_left_args = ['--vcf',vcfname]
    msh_right_args = ['--vcf',rvcfname]
    if args.mapname is not None:
        msh_left_args.extend(['--gen',args.mapname])
        msh_right_args.extend(['--gen',args.mapname])
    msh_left_args.extend(['--dir','left'])
    leftmshfname = vcfname[0:vcfname.rfind(".vcf")]+"_left_results.txt"
    if args.force_override or not isfile(leftmshfname):
        sys.stdout.write(str(args.force_override)+'\t'+str(isfile(leftmshfname))+'\n')
        msh_from_vcf.getmsh(msh_left_args)

    msh_right_args.extend(['--dir','right'])
    rightreversedmshfname = rvcfname[0:rvcfname.rfind(".vcf")]+"_right_results.txt"
    if args.force_override or not isfile(rightreversedmshfname):
        msh_from_vcf.getmsh(msh_right_args)

    rightmshfname = rightreversedmshfname[0:rightreversedmshfname.find("_reversed")]  + rightreversedmshfname[rightreversedmshfname.find("_right_results"):]
    if args.force_override or not isfile(rightmshfname):
        reversefile.reverse(rightreversedmshfname,rightmshfname)
    est_args = [leftmshfname,rightmshfname] + splitArgsForEstimator(args)
    aae_work.run(est_args)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        printf("Usage: runtc.py (VCF name) (optional name for genetic map file)")
