#-------------------------------------------------------------------------------
# Name:        aae_work.py
# Authors:     Jared Knoblauch, Alex Platt, Jody Hey
# Created:     11/07/2019
# Copyright:   (c) Jody Hey 2019
#-------------------------------------------------------------------------------
import sys
import numpy as np
import gzip
import argparse
import math
from msh_from_vcf import getGenMap
from tc_estimator import data,datalist,fitmodel

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("left_file")
    parser.add_argument("right_file")
    parser.add_argument("--start",dest="start",type=int,default=-1)
    parser.add_argument("--end",dest="end",type=int,default=-1)
    parser.add_argument("--n",dest="n_model",type=int,default=-1)
    parser.add_argument("--n0",dest="n0_model",type=int,default=1000)
    parser.add_argument("--mut",dest="mut_rate",type=float,default=1e-8)
    parser.add_argument("--rec",dest="rec_rate",type=float,default=1e-8)
    parser.add_argument("--alledges",dest="alledges",action="store_true",default = False)
    parser.add_argument("--output-all-est",dest="all_est",action="store_true",
                        help=("Output all estimates instead of mean/max"))
    parser.add_argument("--nocache",dest="cache",action="store_false")
    parser.add_argument("--bin",dest="bin",action="store_true")
    parser.add_argument("--round",dest="round",type=int,default=-1)
    parser.add_argument("--full-output",dest="full_out",action="store_true")
    parser.add_argument("--snp-mode",dest="snp_mode",action="store_true",default = False)
    parser.add_argument("--outfn",dest="outfilename",type=str)
    parser.add_argument("--positions",dest="posname",type=str)
    parser.add_argument("--gen",dest="genname",type=str)
    parser.add_argument("--nosquish",dest="squish",action="store_false")
    parser.add_argument("--pos",dest="pos",action="store_true")
    parser.add_argument("--decmode",dest="decmode",action="store_true")
    modelgroup = parser.add_mutually_exclusive_group()
    modelgroup.add_argument("--exp-model",dest="expmodel",type=float,help=("Use "
                            "exponential population model with given growth rate"))
    modelgroup.add_argument("--twophase-model",dest="twophase",nargs=2,type=float,
                            help=("Use two-phase (constant pop size then exp "
                            "growth) model, with growth rate as first argument "
                            "and growth time in generations as second"))
    return parser


def genFunction(val):
    return float(1-math.exp(-2*val))/2

def intw(v):
    if '*' in v:
        return None
    if v == '-2':
        return None
    return int(v)

def floatw(v):
    if v is None or '*' in v:
        return None
    if v == '-2':
        return None
    return float(v)

def intEndCheck(v):
    if v is None:
        return 0,True
    if v[-1] == '*':
        return int(v[:-1]),True
    return int(v),False

def floatEndCheck(v):
    if v is None:
        return 0,True
    if v[-1] == '*':
        return float(v[:-1]),True
    return float(v),False

def p(val,f):
    f.write(str(val)+'\t')

def getl(f):
    l = f.readline()
    try:
        return l.decode()
    except AttributeError:
        return l

def makeData(la1,la2,idx,gen_flag,rec,mu,mod,prp,prg,length_offset,forceleftnone=None,forcerightnone=None):
    d = None
    l1 = (la1[idx] if la1 is not None else ('0*:0*' if gen_flag else '0*'))
    l2 = (la2[idx] if la2 is not None else ('0*:0*' if gen_flag else '0*'))

    if gen_flag:
        d1,df1 = intEndCheck(l1.split(':')[0])
        d2,df2 = intEndCheck(l2.split(':')[0])
        m1,mf1 = floatEndCheck(l1.split(':')[1])
        m2,mf2 = floatEndCheck(l2.split(':')[1])
    else:
        d1,df1 = intEndCheck(l1)
        d2,df2 = intEndCheck(l2)
        m1 = None
        m2 = None
    if df1 and df2:
        return d
    if length_offset > 0:
        #if d1 is not None:
        if la1 is not None:
            pos1 = int(la1[0])
            inc1 = abs(prp-pos1)
            d1 += inc1
            if gen_flag:
                gen1 = float(la1[1])
                m1 += abs(prg-gen1)
        #if d2 is not None and (length_offset == 2 or d1 is None):
        if la2 is not None and (length_offset == 2 or la1 is None):
            pos2 = int(la2[0])
            inc2 = abs(prp-pos2)
            d2 += inc2
            if gen_flag:
                gen2 = float(la2[1])
                m2 += abs(prg-gen2)
    if forceleftnone:
        d1 = m1 = None
    if forcerightnone:
        d2 = m2 = None
    if d1 is None and d2 is None:
        return d
    d = data(dis1 = d1,dis2 = d2, morgans1 = m1, morgans2 = m2, rho1 = rec, rho2 = rec, mu=mu, model=mod, end1=df1, end2=df2)
    return d


def parseFilename(fn):
    fwd = fn.split('/')[-1]
    ac = int(fwd.split('.')[-1])
    return ac

def changeModel(dl,m):
    for d in dl:
        d.model = m

def parseGenLine(l,offset):
    la = l.strip().split()
    a = float(la[0])
    b = float(la[1+offset])
    return a,b

def getActiveIdx(firstla,start_inds):
    idx_list = []
    for i in range(start_inds,len(firstla)):
        if '*' in firstla[i]:
            idx_list.append(i)
    return idx_list

def hasmissing(la,idx_list):
    anymissing = False
    for i in idx_list:
        a = la[i]
        if '*' in a:
            anymissing = True
            break
    return anymissing

##def getGenMap(f):
##    l1 = []
##    l2 = []
##    for line in f:
##        a,b = parseGenLine(line,0)
##        if b != 0 and (len(l1) < 2 or b != l1[-1]):
##            l1.append(a)
##            l2.append(b)
##    return l1,l2

def getGenRate(pos,l1,l2):
    if pos < l1[0]:
        return (l2[1]-l2[0])*.01/(l1[1]-l1[0])
    i = 1
    while i < len(l1)-1 and pos > l1[i]:
        i += 1
    return (l2[i]-l2[i-1])*.01/(l1[i]-l1[i-1])

def geomean(l):
    ll = [float(1 if iiii == 0 else iiii) for iiii in l]
    b = np.array(ll).astype(np.float)
    a = np.log(b)
    return np.exp(a.sum()/len(a))

def fullStr(d):
    return str(d.dis1)+':'+str(d.dis2)+':'+str(d.chi)+':'+str(d.side)+':'+str(d.singlex)

def run_estimator(args):
    parser = createParser()
    args = parser.parse_args(args)
    if args.left_file[-3:] == '.gz':
        fl = gzip.open(args.left_file,'r')
    else:
        fl = open(args.left_file,'r')

    if args.right_file[-3:] == '.gz':
        fr = gzip.open(args.right_file,'r')
    else:
        fr = open(args.right_file,'r')

    position_mode = False
    if args.posname is not None:
        posf = open(args.posname,'r')
        pos_list = [int(l.strip()) for l in posf]
        pos_idx = 0
        position_mode = True
    map_for_rec = False
    if args.genname is not None:
##        print(args.genname)
        genf = open(args.genname,'r')
        g1,g2 = getGenMap(genf,squish=args.squish)
        map_for_rec = True


    snp_mode = args.snp_mode
    start_idx = args.start
    end_idx = args.end

    mu = args.mut_rate
    rec = args.rec_rate

    #n_model = args.n_model
    n0_model = args.n0_model
    #Means left and right lengths should line up
    lengths_set = (not args.alledges or args.pos)
    length_offset = 0
    if args.alledges and not args.pos:
        if args.snp_mode:
            length_offset = 2
        else:
            length_offset = 1
##    sys.stderr.write(str(length_offset)+'\n')  jh 7/11/2019  stopped writing this,  though can be useful for debugging
    ii = 0
    has_genetic_positions = False
    right_done = False
    dl1 = datalist()

    l2 = getl(fr)
    la2 = l2.strip().split()
    has_genetic_positions = (len(la2) > 2 and ':' in la2[2])
    input_length = len(la2)
    start_inds = 1
    if has_genetic_positions:
        start_inds = 2
    if args.n_model == -1:
        n_model = input_length - start_inds
    else:
        n_model = args.n_model
    prev_right_pos = None
    prev_right_gen = None
    prev_right_pos = int(la2[0])
    if has_genetic_positions:
        prev_right_gen = float(la2[1])


    if args.expmodel is not None:
        mc = fitmodel(n=n_model,N0=n0_model,popmodel="expgrowth",nosnp=args.alledges,g=args.expmodel)
    elif args.twophase is not None:
        mc = fitmodel(n=n_model,N0=n0_model,popmodel="twophase",nosnp=args.alledges,g=args.twophase[0],te=args.twophase[1])
    else:
         mc = fitmodel(n=n_model,N0=n0_model,popmodel="constant",nosnp = args.alledges)
    if args.bin:
        dl1.calc_tc_bins(1e-9,10,mu,mc)

    if args.outfilename:
        outf = open(args.outfilename,'w')
    else:
        outf = sys.stdout

    while True:
        #check_left = (not snp_mode or lengths_set or ii != 0)
        check_left = (length_offset <= 1 or ii != 0)
        if check_left:
            try:
                l1 = getl(fl)
                la1 = l1.strip().split()
                cur_left_pos = int(la1[0])
            except IndexError as si:
                break
        else:
            la1 = None
            cur_left_pos = prev_right_pos
        try:
            #if args.alledges or ii != 0:
            #if not lengths_set or ii != 0:
            if length_offset >= 1 or ii != 0:
                l2 = getl(fr)
                la2 = l2.strip().split()
            cur_right_pos = int(la2[0])
            if has_genetic_positions:
                cur_right_gen = float(la2[1])
        except IndexError as si:
            if length_offset < 2:
                break
            cur_right_pos = int(la1[0])
            if has_genetic_positions:
                cur_right_gen = float(la1[1])
            right_done = True
            la2 = None
        if start_idx != -1 and ii < start_idx:
            ii += 1
            continue
        elif end_idx != -1 and ii >= end_idx:
            break

        est_list_ml = []

        if map_for_rec:
            if length_offset > 0:
                rec_phys_pos = round((cur_right_pos+prev_right_pos)/2)
            else:
                rec_phys_pos = cur_right_pos
            recperbp = getGenRate(rec_phys_pos,g1,g2)
        else:
            recperbp = rec
        out_pos = 0

        hasleftmissing = (la1 == None)
        hasrightmissing = (la2 == None)
        est_str = ''
        #if args.nosnp:
        dl1.clear()
        if not args.decmode:
            for i in range(start_inds,input_length):
                d = makeData(la1,la2,i,has_genetic_positions,recperbp,mu,mc,
                            prev_right_pos,prev_right_gen,length_offset,forceleftnone=hasleftmissing,forcerightnone=hasrightmissing)
                if d is not None:
                    dl1.append(d)
            if len(dl1) != 0:
                chi_list = []
                for d in dl1:
                    if d is not None:
                        chi_list.append(d.chi)
                chi_geomean = geomean(chi_list)
                est_list_ml= dl1.estimate_tc_cache(cache=args.cache,round=args.round,bin=args.bin)
                if args.all_est:
                    est_str = '\t'.join(map(str,est_list_ml))+'\n'
                else:
                    if args.alledges:
                        est_all_ml = geomean(est_list_ml)
                    else:
                        est_all_ml = max(est_list_ml)
                    est_str =  "{:.3f}\n".format(est_all_ml)
##                    est_str = (str(est_all_ml)+'\n')

                if length_offset > 0:
                    outf.write(str(prev_right_pos)+'\t'+est_str)
                else:
                    outf.write(str(cur_right_pos)+'\t'+est_str)
            else:
                sys.stderr.write("pos %d: no valid data\n"%prev_right_pos)
                est_str = ''
        else:
            input_length = len(la2)
            for i in range(start_inds,input_length):
                d = makeData(la1,la2,i,has_genetic_positions,recperbp,mu,mc,
                            prev_right_pos,prev_right_gen,length_offset,forceleftnone=hasleftmissing,forcerightnone=hasrightmissing)
                if d is not None:
                    dl1.append(d)
            k_val = input_length - start_inds
            if len(dl1) != 0:
                dl1.estimate_tc_kgt1()
                est_str = (str(cur_right_pos)+'\t'+ "{:.3f}".format(dl1.tcest)+'\t'+str(k_val)+'\n')
##                est_str = (str(cur_right_pos)+'\t'+ str(dl1.tcest)+'\t'+str(k_val)+'\n')
                outf.write(est_str)
            else:
                outf.write(str(cur_right_pos)+'\t-1\t'+str(k_val)+'\n')
        ii += 1
        prev_right_pos = cur_right_pos
        if has_genetic_positions:
            prev_right_gen = cur_right_gen
        if right_done:
            break
##    jh 7/11/2019  stopped writing this,  though can be useful in development/debugging
##    if args.cache and args.alledges and not args.decmode:
##        sys.stderr.write("Cache: %d of %d hits (%f rate)\n" % (dl1.cache_hits,dl1.cache_total,float(dl1.cache_hits)/max(float(dl1.cache_total),1)))

if __name__ == "__main__":
    run_estimator(sys.argv[1:])
