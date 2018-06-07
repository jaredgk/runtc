import sys
import msh_from_vcf
import reversefile
import aae_work
import argparse

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
    parser.add_argument("--force")
    return parser



def main(argv):
    vcfname = argv[0]

    rvcfname = vcfname[0:vcfname.rfind(".vcf")] + "_reversed.vcf"

    reversefile.reverse(argv[0],rvcfname)


    if len(argv) > 1:
        leftmshfname = msh_from_vcf.getmsh(["--vcf",vcfname,"--dir","left","--gen",argv[1]])
    else:
        leftmshfname = msh_from_vcf.getmsh(["--vcf",vcfname,"--dir","left"])

    rightreversedmshfname = msh_from_vcf.getmsh(["--vcf",rvcfname,"--dir","right","--gen",argv[1]])

    rightmshfname = rightreversedmshfname[0:rightreversedmshfname.find("_reversed")]  + rightreversedmshfname[rightreversedmshfname.find("_right_results"):]

    reversefile.reverse(rightreversedmshfname,rightmshfname)

    aae_work.run([leftmshfname,rightmshfname,"--mut","1e-6","--n0","1000","--t0","--n","100"])

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        printf("Usage: runtc.py (VCF name) (optional name for genetic map file)")
