import sys
from tc_estimator import data,datalist,fitmodel
import numpy as np
import gzip
import argparse
import math
from msh_from_vcf import getGenMap

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
    parser.add_argument("--alpha",dest="alpha",action="store_true",default = False)
    parser.add_argument("--nocache",dest="cache",action="store_false")
    parser.add_argument("--bin",dest="bin",action="store_true")
    parser.add_argument("--round",dest="round",type=int,default=-1)
    parser.add_argument("--mod-gen",dest="mod_gen",action="store_true")
    parser.add_argument("--full-output",dest="full_out",action="store_true")
    parser.add_argument("--snp-mode",dest="snp_mode",action="store_true",default = False)
    parser.add_argument("--side-check",dest="side_check",action="store_true")
    parser.add_argument("--outfn",dest="outfilename",type=str)
    parser.add_argument("--positions",dest="posname",type=str)
    parser.add_argument("--gen",dest="genname",type=str)
    parser.add_argument("--nosquish",dest="squish",action="store_false")
    parser.add_argument("--pos",dest="pos",action="store_true")
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
    #sys.stderr.write('d\t'+str(v)+'\t'+str(float(v))+'\n')
    if v is None or '*' in v:
        return None
    if v == '-2':
        return None
    return float(v)

def p(val,f):
    f.write(str(val)+'\t')

def getl(f):
    l = f.readline()
    try:
        return l.decode()
    except AttributeError:
        return l

def makeData(la1,la2,idx,gen_flag,rec,mu,mod,prp,prg,length_offset,mod_gen=False,forceleftnone=None,forcerightnone=None):
    d = None
    l1 = (la1[idx] if la1 is not None else ('-2:-2' if gen_flag else '-2'))
    l2 = (la2[idx] if la2 is not None else ('-2:-2' if gen_flag else '-2'))

    if gen_flag:
        d1 = intw(l1.split(':')[0])
        d2 = intw(l2.split(':')[0])
        m1 = floatw(l1.split(':')[1])
        m2 = floatw(l2.split(':')[1])
    else:
        d1 = intw(l1)
        d2 = intw(l2)
        m1 = None
        m2 = None
        #m1 = (d1*rec if d1 is not None else None)
        #m2 = (d2*rec if d2 is not None else None)
    if d1 is None and d2 is None:
        return d
    if length_offset > 0:
        if d1 is not None:
            pos1 = int(la1[0])
            inc1 = abs(prp-pos1)
            d1 += inc1
            if gen_flag:
                gen1 = float(la1[1])
                m1 += abs(prg-gen1)
        if d2 is not None and (length_offset == 2 or d1 is None):
            pos2 = int(la2[0])
            inc2 = abs(prp-pos2)
            d2 += inc2
            if gen_flag:
                gen2 = float(la2[1])
                m2 += abs(prg-gen2)
    if mod_gen and m1 is not None:
        m1 = genFunction(m1)
    if mod_gen and m2 is not None:
        m2 = genFunction(m2)
    if forceleftnone:
        d1 = m1 = None
    if forcerightnone:
        d2 = m2 = None
    if d1 is None and d2 is None:
        return d
    d = data(dis1 = d1,dis2 = d2, morgans1 = m1, morgans2 = m2, rho1 = rec, rho2 = rec, mu=mu, model=mod)
    return d


def parseFilename(fn):
    fwd = fn.split('/')[-1]
    nec = 0
    n = 0
    ac = int(fwd.split('.')[-1])
    return ac

def changeModel(dl,m):
    for d in dl:
        d.model = m

def parseGenLine(l,offset):
    la = line.strip().split()
    a = float(la[0])
    b = float(la[1+offset])
    return a,b

def getActiveIdx(firstla,start_inds):
    idx_list = []
    for i in range(start_inds,len(firstla)):
        if firstla[i][0:2] != '-2':
            idx_list.append(i)
    return idx_list

def hasmissing(la,idx_list):
    anymissing = False
    for i in idx_list:
        a = la[i]
        if a[0:2] == '-2':
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

def fullStr(dl):
    d = dl[0]
    left_side = ''
    if d.side != 2:
        left_side = str(d.dis1)+','+str(d.morgans1)
    else:
        left_side = '-1,-1'
    if d.side != 1:
        right_side = str(d.dis2)+','+str(d.morgans2)
    else:
        right_side = '-1,-1'
    return left_side+','+right_side+','+str(d.chi)+','+str(dl.tcest)

def run_estimator(args):
    parser = createParser()
    args = parser.parse_args(args)
    #fnl = str(sys.argv[1])
    if args.left_file[-3:] == '.gz':
        fl = gzip.open(args.left_file,'r')
    else:
        fl = open(args.left_file,'r')



    #fnr = str(sys.argv[2])
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
        print(args.genname)
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
    lengths_set = (not args.alpha or args.pos)
    length_offset = 0
    if args.alpha and not args.pos:
        if args.snp_mode:
            length_offset = 2
        else:
            length_offset = 1
    sys.stderr.write(str(length_offset)+'\n')
    ii = 0
    has_genetic_positions = False
    right_done = False
    dl1 = datalist()

    l2 = getl(fr)
    la2 = l2.strip().split()
    has_genetic_positions = (':' in la2[2])
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

    #if args.singleton:
    #    la1 = getl(fl).strip().split()

    if args.side_check:
        idx_list = getActiveIdx(la2,start_inds)

    mc = fitmodel(n=n_model,N0=n0_model,popmodel="constant",nosnp = args.alpha)

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
            #if args.alpha or ii != 0:
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
            #recperbp = (cur_right_gen-prev_right_gen)/float(cur_right_pos-prev_right_pos)
        else:
            recperbp = rec
##        if region_mode:
##            if position_mode:
##                outf.write(str(pos_list[pos_idx]))
##            else:
##                outf.write(str(prev_right_pos)+'-'+str(cur_right_pos))
##        else:
##            outf.write(str(prev_right_pos))
        out_pos = 0

        hasleftmissing = (la1 == None or (args.side_check and hasmissing(la1,idx_list)))
        hasrightmissing = (la2 == None or (args.side_check and hasmissing(la2,idx_list)))
        est_str = ''
        #if args.nosnp:
        dl1.clear()
        for i in range(start_inds,input_length):
            d = makeData(la1,la2,i,has_genetic_positions,recperbp,mu,mc,
                         prev_right_pos,prev_right_gen,length_offset,mod_gen=args.mod_gen,forceleftnone=hasleftmissing,forcerightnone=hasrightmissing)
            if d is not None:
                dl1.append(d)
        if len(dl1) != 0:
            chi_list = []
            for d in dl1:
                chi_list.append(d.chi)
            chi_geomean = geomean(chi_list)
            est_list_ml= dl1.estimate_tc_cache(cache=args.cache,round=args.round,bin=args.bin)
            if args.alpha:
                est_all_ml = geomean(est_list_ml)
            else:
                est_all_ml = max(est_list_ml)
            est_str += ("\t%.4g"%recperbp)
            est_str += ("\t%d\t%.4g"%(len(dl1),chi_geomean))
            est_str += ('\t'+str(est_all_ml)+'\n')
            if length_offset > 0:
                outf.write(str(prev_right_pos)+'\t'+est_str)
            else:
                outf.write(str(cur_right_pos)+'\t'+est_str)
        else:
            sys.stderr.write("pos %d: no valid data\n"%prev_right_pos)
            est_str = ''

        ii += 1
        prev_right_pos = cur_right_pos
        if has_genetic_positions:
            prev_right_gen = cur_right_gen
        if right_done:
            break
    if args.cache and args.alpha:
        sys.stderr.write("Cache: %d of %d hits (%f rate)\n" % (dl1.cache_hits,dl1.cache_total,float(dl1.cache_hits)/max(float(dl1.cache_total),1)))

if __name__ == "__main__":
    run_estimator(sys.argv[1:])
