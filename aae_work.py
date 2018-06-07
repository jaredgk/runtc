import sys
from tc_estimator import data,datalist,fitmodel
import numpy as np
import gzip
import argparse
import math

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("left_file")
    parser.add_argument("right_file")
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

def makeData(la1,la2,idx,gen_flag,rec,mu,mod,prp,prg,region_mode,mod_gen=False):
    d = None
    l1 = (la1[idx] if la1 is not None else '-2:-2')
    l2 = (la2[idx] if la2 is not None else '-2:-2')
    if gen_flag:
        d1 = intw(l1.split(':')[0])
        d2 = intw(l2.split(':')[0])
        m1 = floatw(l1.split(':')[1])
        m2 = floatw(l2.split(':')[1])
    else:
        d1 = intw(l1)
        d2 = intw(l2)
        m1 = (d1*rec if d1 is not None else None)
        m2 = (d2*rec if d2 is not None else None)
    if mod_gen and m1 is not None:
        m1 = genFunction(m1)
    if mod_gen and m2 is not None:
        m2 = genFunction(m2)
    if d1 is None and d2 is None:
        return d
    if d1 is not None:
        pos1 = int(la1[0])
        inc1 = abs(prp-pos1)
        d1 += inc1
        if gen_flag:
            gen1 = float(la1[1])
            m1 += abs(prg-gen1)
        else:
            gen1 = inc1*rec
            m1 += abs(prg-gen1)
    if d2 is not None and (not region_mode or d1 is None):
        pos2 = int(la2[0])
        inc2 = abs(prp-pos2)
        d2 += inc2
        if gen_flag:
            gen2 = float(la2[1])
            m2 += abs(prg-gen2)
        else:
            gen2 = inc2*rec
            m2 += abs(prg-gen2)
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

def getGenMap(f):
    l1 = []
    l2 = []
    for line in f:
        a,b = parseGenLine(line,0)
        if b != 0 and (len(l1) < 2 or b != l1[-1]):
            l1.append(a)
            l2.append(b)
    return l1,l2

def getGenRate(pos,l1,l2):
    if pos < l1[0]:
        return (l2[1]-l2[0])*.01/(l1[1]-l1[0])
    i = 1
    while i < len(l1)-1 and pos < l1[i]:
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
    return left_side+','+right_side+','+str(dl.tcest)

def run(args):
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

    region_mode = args.region_mode


    start_idx = args.start
    end_idx = args.end



    mu = args.mut_rate
    rec = args.rec_rate

    n_model = args.n_model
    n0_model = args.n0_model


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
    prev_right_pos = None
    prev_right_gen = None
    prev_right_pos = int(la2[0])
    if has_genetic_positions:
        prev_right_gen = float(la2[1])

    mc = fitmodel(n=n_model,N0=n0_model,popmodel="c",nosnp = args.nosnp)

    while True:
        check_left = (region_mode or ii != 0)
        if check_left:
            try:
                l1 = getl(fl)
                la1 = l1.strip().split()
                cur_left_pos = int(la1[0])
            except StopIteration as si:
                break
        else:
            la1 = None
            cur_left_pos = prev_right_pos
        try:
            l2 = getl(fr)
            la2 = l2.strip().split()
            cur_right_pos = int(la2[0])
            if has_genetic_positions:
                cur_right_gen = float(la2[1])
        except StopIteration as si:
            if region_mode:
                break
            right_done = True
            la2 = None
        if start_idx != -1 and ii < start_idx:
            ii += 1
            continue
        elif start_idx != -1 and ii >= end_idx:
            break

        est_list_ml = []
        #if not right_done:
        #    cur_right_pos = int(la2[0])
        #else: #only gets here in snp mode
        #    cur_right_pos = int(la1[0])

        if region_mode:
            sys.stdout.write(str(prev_right_pos)+'-'+str(cur_right_pos))
        else:
            sys.stdout.write(str(prev_right_pos))
        ppp,qqq = 0,0
        if la1 is not None:
            ppp = int(la1[0])
        if la2 is not None:
            qqq = int(la2[0])
        for i in range(start_inds,input_length):
            d = makeData(la1,la2,i,has_genetic_positions,rec,mu,mc,
                         prev_right_pos,prev_right_gen,region_mode,mod_gen=args.mod_gen)
            if d is None:
                if args.full_out:
                    sys.stdout.write('\t-1,-1,-1,-1,-1')
                continue
            if len(dl1) == 0:
                dl1.append(d)
            else:
                dl1[0] = d
            if len(dl1) != 0:
                dl1.estimate_tc_cache(cache_est=args.cache_est)
                est_list_ml.append(dl1.tcest)
                if args.full_out:
                    sys.stdout.write('\t'+fullStr(dl1))
        est_all_ml = geomean(est_list_ml)
        sys.stdout.write('\t'+str(est_all_ml))
        sys.stdout.write('\n')
        ii += 1
        prev_right_pos = cur_right_pos

if __name__ == "__main__":
    run(sys.argv[1:])
