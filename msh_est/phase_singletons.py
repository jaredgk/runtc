#-------------------------------------------------------------------------------
# Name:        phase_singletons.py
# Authors:     Jared Knoblauch, Alex Platt, Jody Hey
# Created:     11/07/2019
# Copyright:   (c) Jody Hey 2019
#-------------------------------------------------------------------------------
import sys
import argparse
import gzip
from os.path import isfile
import os
import runtc


def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcfname",type=str,help="vcf or vcf.gz input file")
##    parser.add_argument("--mutrate",dest="mutrate",type=float,default=1e-8)  not needed
    parser.add_argument("--gzip",dest="gzip_check",action="store_true",help="gzip temporary msh files")
    parser.add_argument("--statfile",dest="statfile",type=str,help=argparse.SUPPRESS)
    parser.add_argument("--keep",dest='keepfiles',action="store_true",help=argparse.SUPPRESS)
    parser.add_argument("--c-msh",dest="c_msh",action="store_true")
    return parser

def findSingleton(la):
    alleles = []
    ac = 0
    an = 0
    for i in range(9,len(la)):
        laa = la[i].split(':')[0]
        gts = laa.replace('|','/').split('/')
        #sc = len(laa.replace('|','/').split('/'))
        for j in range(len(gts)):
            try:
                g1 = int(gts[j])
            except ValueError:
                return []
            alleles.append(g1)
            ac += g1
            an += 1
    if alleles.count(1) == 1:
        return (alleles.index(1)//2)+9,1
    elif alleles.count(0) == 1:
        return (alleles.index(0)//2)+9,0
    raise Exception("Allele count %d not singleton"%(alleles.count(1)))
    #return alleles

#def findSingleton(la):
#   for i in range(9,len(la)):
#        if '1' in la[i][0:3]:
#            return i
#    raise Exception("No singleton found")


def parseChi(f1,f2):
    if ':' in f1:
        fa1 = f1.split(':')
        fa2 = f2.split(':')
        phys1 = (int(fa1[0]) if '*' not in fa1[0] else int(fa1[0][:-1]))
        phys2 = (int(fa2[0]) if '*' not in fa2[0] else int(fa2[0][:-1]))
        gen1 = (float(fa1[1]) if '*' not in fa1[1] else float(fa1[1][:-1]))
        gen2 = (float(fa2[1]) if '*' not in fa2[1] else float(fa2[1][:-1]))
        phys = phys1+phys2
        gen = gen1 + gen2
        return phys, gen
        #chi = float(phys1)*mutrate+float(phys2)*mutrate+gen1+gen2
        #return chi
    phys1 = (int(f1) if '*' not in f1 else int(f1[:-1]))
    phys2 = (int(f2) if '*' not in f2 else int(f2[:-1]))
    return phys1+phys2, None

def computeChi(f1,f2,mutrate):
    phys, gen = parseChi(f1,f2)
    if gen is None:
        return float(phys)*mutrate
    return float(phys)*mutrate + gen

def getStatString(current_left_la,current_right_la,geno,singleton_allele):
    #if geno[0] = '0':
    if int(geno[0]) == singleton_allele:
        physl,genl = parseChi(current_left_la[-2],current_right_la[-2])
        physr,genr = parseChi(current_left_la[-1],current_right_la[-1])
        idx1 = -2
        idx2 = -1
    else:
        physr,genr = parseChi(current_left_la[-2],current_right_la[-2])
        physl,genl = parseChi(current_left_la[-1],current_right_la[-1])
        idx1 = -1
        idx2 = -2
    raw_vals = current_left_la[idx1]+'\t'+current_right_la[idx1]+'\t'+current_left_la[idx2]+'\t'+current_right_la[idx2]
    if genl is None:
        return str(physl)+'\t'+str(physr)+'\t'+raw_vals
    return str(physl)+'\t'+str(physr)+'\t'+str(genl)+'\t'+str(genr)+'\t'+raw_vals

def phaseSingleton(orig_geno,current_left_la,current_right_la,mutrate,singleton_allele):
    ga = '0|1' if singleton_allele == 1 else '1|0'
    gb = '1|0' if singleton_allele == 1 else '0|1'
    chi_1 = computeChi(current_left_la[-2],current_right_la[-2],mutrate)
    chi_2 = computeChi(current_left_la[-1],current_right_la[-1],mutrate)
    #new_geno = ('0|1' if chi_1 > chi_2 else '1|0')
    new_geno = (ga if chi_1 > chi_2 else gb)
    #new_geno = ('0|1' if chi_1 < chi_2 else '1|0')
    return new_geno, (new_geno[0] == orig_geno[0])


def phase_with_lengths(sysargs):
    parser = createParser()
    args = parser.parse_args(sysargs)
    vcfname = args.vcfname
    vcftag = vcfname[0:vcfname.rfind(".vcf")]
    rvcfname = vcftag+"_reversed.vcf"
    if vcfname[-6:] == 'vcf.gz':
        rvcfname += ".gz"

    ## set options for runtc.makemshfiles()
    args.outmsh = vcftag + "_singleton_phasing"
    #args.revname = rvcfname
    args.revname = None
    args.lengths_only = True
    args.exc_sing = True
##    args.revname = None
    args.no_pypy = False
    args.force_overwrite = True
    args.rec_rate = 1e-8
    args.mut_rate = 0.0  ## mutation does not matter for phasing singletons
    args.mapname = None
    args.subname = None
    args.squish = True
    args.round = 3
    args.all_est  = False
    args.posname = None
    args.k_all = None
    args.k_range = None
    args.dt_exp  = None
    args.random_n = -1
    args.alledges = False
    args.n0_model = 10000
    args.cache = None
    args.bin = None
    args.outfn = None
    args.snp_mode = None
    args.expmode = None
    args.expmodel = None
    args.twophase = None
    args.est_seed = None
    args.print_indiv = False
    args.k_list = None
    args.k1 = True
    args.singleton_phase = True
    args.genidx = None

    runtcestargs = runtc.makemshfiles(args)
    args.left_lengths = runtcestargs[0]
    args.right_lengths = runtcestargs[1]
    vcf_decode = False
    if args.vcfname[-3:] == '.gz':
        vcf = gzip.open(args.vcfname,'r')
        vcf_decode = True
    else:
        vcf = open(args.vcfname,'r')
    left_decode = False
    right_decode = False
    if args.left_lengths[-3:] == '.gz':
        left_f = gzip.open(args.left_lengths,'r')
        left_decode = True
    else:
        left_f = open(args.left_lengths,'r')
    if args.right_lengths[-3:] == '.gz':
        right_f = gzip.open(args.right_lengths,'r')
        right_decode = True
    else:
        right_f = open(args.right_lengths,'r')

    if args.statfile is not None:
        statf = open(args.statfile,'w')

    if left_decode:
        left_line = left_f.readline().decode()
    else:
        left_line = left_f.readline()
    if right_decode:
        right_line = right_f.readline().decode()
    else:
        right_line = right_f.readline()

    current_left_la = left_line.strip().split()
    current_right_la = right_line.strip().split()
    current_pos = int(current_left_la[0])
    for line in vcf:
        if vcf_decode:
            line = line.decode()
        if line[0] == '#':
            sys.stdout.write(line)
            continue
        la = line.strip().split()
        pos = int(la[1])
        if pos != current_pos:
            sys.stdout.write(line)
            continue
        singleton_idx,singleton_allele = findSingleton(la)
        for i in range(len(la)):
            if i != singleton_idx:
                sys.stdout.write(la[i])
            else:
                new_geno, changed = phaseSingleton(la[i],current_left_la,current_right_la,args.mut_rate,singleton_allele)
                sys.stdout.write(new_geno)
                if args.statfile is not None:
                    statf.write(la[1]+'\t'+str(singleton_idx)+'\t'+getStatString(current_left_la,current_right_la,la[i],singleton_allele)+'\n')
            if i != len(la)-1:
                sys.stdout.write('\t')
            else:
                sys.stdout.write('\n')
        if left_decode:
            left_line = left_f.readline().decode()
        else:
            left_line = left_f.readline()
        if right_decode:
            right_line = right_f.readline().decode()
        else:
            right_line = right_f.readline()
        current_left_la = left_line.strip().split()
        current_right_la = right_line.strip().split()

        try:
            current_pos = int(current_left_la[0])
        except IndexError:
            current_pos = 0
            continue
    if not args.keepfiles:
        left_f.close()
        right_f.close()
        if isfile(runtcestargs[0]):
            os.remove(runtcestargs[0])
        if isfile(runtcestargs[1]):
            os.remove(runtcestargs[1])
        if args.revname is not None and isfile(args.revname):
            os.remove(args.revname)
    return

if __name__ == '__main__':
    phase_with_lengths(sys.argv[1:])


