#-------------------------------------------------------------------------------
# Name:        msh_from_vcf.py
# Authors:     Jared Knoblauch, Alex Platt, Jody Hey
# Created:     11/07/2019
# Copyright:   (c) Jody Hey 2019
#-------------------------------------------------------------------------------

import sys
import argparse
import gzip
import math
import gc
#from memory_profiler import profile


def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf",dest="vcfname",help="Input VCF filename, either uncompressed or gzipped")
    parser.add_argument("--gen",dest="genname",help="Name of map file")
    parser.add_argument("--sub",dest="subname",help="Name of file with subsample indices")
    parser.add_argument("--out",dest="outname",help="Name of output msh file")
    parser.add_argument("--map-col-idx",dest="genidx",type=int,nargs=2,default=[1,2],help="Two numbers, first is column of physical position in map file and second is column of genetic position (default: 1 2)")
    parser.add_argument("--nosquish",dest="squish",action="store_false",help="Read all rows in map file, not just rows with differing genetic distances")
    parser.add_argument("--round",dest="round",type=int,default=-1,help="Round floats to this many significant figures")
    runtype = parser.add_mutually_exclusive_group()
    runtype.add_argument("--k1",dest="k1",action="store_true",help=("Return estimates for singletons"))
    runtype.add_argument("--k-all",dest="k_all",action="store_true",help="Return estimates for every variant")
    runtype.add_argument("--k-range",dest="k_range",nargs=2,type=int,help="Return estimates for every variant with allele count in given range")
    runtype.add_argument("--k-list",dest="k_list",type=int,nargs="+",help=("Return estimates for every variant with allele count in provided list"))
    runtype.add_argument("--singleton",dest="singleton",action="store_true",help="Singleton mode: Generate two msh values, one for each haplotype of an individual with a singleton")
    parser.add_argument("--singleton-phase",action="store_true",help=("Used to return two lengths for singleton phasing on individual with singleton"))    
    parser.add_argument("--alledges-singleton-only",dest="alledges_sing_only",action="store_true",help=("Output lengths for alledges estimator at singleton sites only"))
    parser.add_argument("--positions",dest="posname",type=str,help="File with positions for alledges estimates")
    parser.add_argument("--exclude-singletons",dest="exc_sing",action="store_true",help="DO NOT Allow singletons to terminate MSH lengths") ## changed from --include-singletons JH 7/11/2019
    parser.add_argument("--revpos",dest="revpos",action="store_true",help="If using --positions, indicate this is for right length generation and input VCF is reversed")
    #parser.add_argument("--k",dest="k_val",type=int,help=("Return outgroup "
    #                    "lengths for every k-ton"))
    parser.add_argument("--dt-exp",dest="dt_exp",nargs=2,type=int)
    return parser

def splitAllelesAll(la):
    alleles = []
    ac = 0
    an = 0
    for i in range(9,len(la)):
        laa = la[i].split(':')[0]
        gts = laa.replace('|','/').split('/')
        for j in range(len(gts)):
            try:
                g1 = int(gts[j])
            except ValueError:
                return []
            alleles.append(g1)
            ac += g1
            an += 1
    return alleles

def writeToFile(outf,write_line,compress_out):
    if compress_out:
        write_line = write_line.encode()
    outf.write(write_line)

def subsampToIdx(la,sub_list):
    sub_list.sort()
    current_idx = 0
    sub_idx = 0
    idx_list = []
    for i in range(9,len(la)):
        laa = la[i].split(':')[0]
        gts = laa.replace('|','/').split('/')
        for j in range(len(gts)):
            if current_idx == sub_list[sub_idx]:
                idx_list.append([i,j])
                sub_idx += 1
            current_idx += 1
            if sub_idx == len(sub_list):
                break
        if sub_idx == len(sub_list):
            break
    if sub_idx != len(sub_list):
        raise Exception("Value %d in subsample list is too large for data" %(sub_list[sub_idx]))
    return idx_list

def splitAllelesSub(la,idx_list):
    alleles = []
    for i in range(len(idx_list)):
        ip = idx_list[i]
        laa = la[ip[0]].split(':')[0]
        gts = laa.replace('|','/').split('/')
        try:
            geno = int(gts[ip[1]])
        except ValueError:
            return []
        alleles.append(geno)
    return alleles

def splitAlleles(la,idx_list=None):
    if idx_list is None:
        return splitAllelesAll(la)
    else:
        return splitAllelesSub(la,idx_list)

def getVectors(a_prev,d_prev,cur_row,k):
    p = k+1
    q = k+1
    a = []
    d = []
    b = []
    e = []
    for i in range(len(a_prev)):
        p = max(d_prev[i],p)
        q = max(d_prev[i],q)
        if cur_row[a_prev[i]] == 0:
            a.append(a_prev[i])
            d.append(p)
            p = 0
        else:
            b.append(a_prev[i])
            e.append(q)
            q = 0
    a.extend(b)
    d.extend(e)
    return a,d

def msh(a,d,pos_list,pos,sample_count):
    if a is None:
        return ['0*' for i in range(sample_count)]
    l = len(a)
    y_msh = [0 for i in range(l)]
    site_msh = [0 for i in range(l)]

    for i in range(l):
        if i == l-1:
            c_idx = d[l-1]
        elif i == 0:
            c_idx = d[1]
        else:
            c_idx = min(d[i],d[i+1])
        if c_idx == 0:
            y_msh[i] = str(abs(pos-pos_list[0]))+'*'
        else:
            y_msh[i] = abs(pos - pos_list[c_idx-1])
    for a_i,a_v in enumerate(a):
        site_msh[a_v] = y_msh[a_i]
    return site_msh

def getDSingle(a,d,idx,sample_count):
    if a is None:
        return -1
    aidx = a.index(idx)

    if aidx == sample_count - 1:
        return d[sample_count-1]
    elif aidx == 0:
        return d[1]
    return min(d[aidx],d[aidx+1])

def reconstructDVector(d,idx):
    d_prime = [0 for i in d]
    for i in range(idx-1,-1,-1):
        d_prime[i] = max(d[i+1:idx+1])
    for i in range(idx+1,len(d)):
        d_prime[i] = max(d[idx+1:i+1])
    return d_prime

def reconstructDMatrix(a,d):
    d_prime = []
    for i in range(len(d)):
        d_prime.append(reconstructDVector(d,i))
    real_d_prime = [[0 for i in range(len(a))] for j in range(len(a))]
    for i in range(len(a)):
        for j in range(len(a)):
            real_d_prime[a[i]][a[j]] = d_prime[i][j]
    return real_d_prime




#def printAD(a,d,idx,idx_list):
def printAD(a,d):
    #d_prime = reconstructDVector(d,idx)
    #sys.stdout.write('\t'.join(map(str,idx_list))+'\n')
    #sys.stdout.write('\t'.join(map(str,d_prime))+'\n')
    sys.stdout.write('\t'.join(map(str,d))+'\n')
    sys.stdout.write('\t'.join(map(str,a))+'\n')
    #a_prime = [(aa in idx_list) for aa in a]
    #sys.stdout.write('\t'.join(map(str,a_prime))+'\n')
    #for i,rec in enumerate(zip(map(str,d),map(str,d_prime),map(str,a),map(str,a_prime))):
    #    sys.stdout.write(str(i)+'\t' +'\t'.join(rec)+'\n')

def getDKton(a,d,idx,sample_count,idx_list):
    if a is None:
        return -1
    if idx not in idx_list:
        raise Exception("Index %d not in index list %s" %(a[idx],str(idx_list)))
    d_low = 0
    d_hi = 0
    a_idx = a.index(idx)
    hi_valid = False
    low_valid = False
    if a_idx != 0:
        d_low = d[a_idx]
        low_valid = True
    if a_idx != len(a)-1:
        d_hi = d[a_idx+1]
        hi_valid = True
    i = a_idx-1
    if a_idx > 0 and a[a_idx-1] in idx_list:
        low_valid = False
        while i >= 1:
            d_low = max(d_low,d[i])
            if a[i-1] not in idx_list:
                low_valid = True
                break
            i -= 1
    i = a_idx+2
    if a_idx < sample_count-1 and a[a_idx+1] in idx_list:
        hi_valid = False
        while i < sample_count:
            d_hi = max(d_hi,d[i])
            if a[i] not in idx_list:
                hi_valid = True
                break
            i += 1
    if not hi_valid:
        return d_low
    if not low_valid:
        return d_hi
    return min(d_hi,d_low)



def roundSig(f,n):
    """
    Rounds chi to n significant digits

    """
    if n == -1:
        return f
    if n < 1:
        raise Exception("N value %d is not valid" % (n))
    if f == 0.0:
        return 0.0
    return round(f,-int(math.floor(math.log10(abs(f))))+(n-1))

def roundStar(f,n):
    sf = str(f)
    if sf[-1] == '*':
        return str(roundSig(int(sf[:-1]),n))+'*'
    return str(roundSig(f,n))

def parseGenLine(la,lb,offset):
    a = float(la[0])
    b = float(la[1+offset])
    return a,b

def getGenMap(f,pos_idx=1,gen_idx=2,squish=False):
    l1 = []
    l2 = []
    for line in f:
        if line[0] == '#' or sum(c.isdigit() for c in line) < sum(c.isalpha() for c in line): # skip a line with more letters than numbers # first line may be column headers
            continue
        la = line.strip().split()
        a = float(la[pos_idx-1])
        b = float(la[gen_idx-1])
        if not squish or len(l2) == 0 or (len(l2) > 0 and l2[-1] != b):
            l1.append(a)
            l2.append(b)
##    sys.stderr.write('\t'.join(map(str,[l1[0],l2[0],l1[-1],l2[-1]]))+'\n')
##    sys.stderr.write(str(len(l1))+'\n')
    return l1,l2

def getGenPos(pos,l1,l2,cm=True):
    rate = 1.0
    if cm:
        rate = .01
    if float(pos) < l1[0]:
        gr = (l2[1]-l2[0])*rate/(l1[1]-l1[0])
        return gr*(pos-l1[0])
    i = 1
    while i < len(l1)-1 and float(pos) > float(l1[i]):
        i += 1
    gr = (l2[i]-l2[i-1])*rate/float(l1[i]-l1[i-1])
    return rate*l2[i-1]+gr*float(pos-l1[i-1])

def getSingletonIdx(alleles):
    t = alleles.count(1)
    if t == 1:
        return (alleles.index(1)//2)*2
    else:
        return (alleles.index(0)//2)*2


def readsubfile(subname):
    """
        reads a file that contains a series of integers
        returns a list of integers
        should work regardless of what separates the integers
    """
    sf = open(subname,'r')
    idx_list = []
    inint = False
    cc = ''
    while True:
        c = sf.read(1)
        if c:
            if c.isdigit():
                if inint:
                    cc += c
                else:
                    cc = c
                inint = True
            else:
                if inint:
                    idx_list.append(int(cc))
                    cc = ''
                    inint = False
        else:
            if inint:
                idx_list.append(int(cc))
            break
    sf.close()
    return idx_list

def getKIdxList(alleles,k):
    ac = sum(alleles)
    target_allele = 0
    out_idx = []
    if ac == k:
        target_allele = 1
    for i in range(len(alleles)):
        if alleles[i] == target_allele:
            out_idx.append(i)
    return out_idx

def getMshString(args,a,d,out_range,pos_list,gen_list,pos,gpos,sample_count,idx_list=None):
    out_string = ''
    gen_flag = (gen_list is not None)
    out_string += str(pos)
    if gen_flag:
        out_string += '\t'+str(roundSig(gpos,args.round))
    for i in out_range:
        if idx_list is None or len(idx_list) <= 1:
            d_idx = getDSingle(a,d,i,sample_count)
        else:
            d_idx = getDKton(a,d,i,sample_count,idx_list)
            #sys.stderr.write(str(a)+'\t'+str(d_idx)+'\n')
        #out_string += '\t'+str(d_idx)+'\t'+str(aidx)
        #sys.stderr.write(str(d_idx)+'\n')
        if d_idx == -1:
            pos_len = '0*'
        elif d_idx == 0:
            pos_len = str(abs(pos-pos_list[d_idx]))+'*'
        else:
            pos_len = str(abs(pos-pos_list[d_idx-1]))
        out_string += '\t'+pos_len
        aaa = (0 if a is None else a.index(i))
        #out_string += '\t'+str(d_idx)+'\t'+str(aaa)
        if gen_flag:
            if d_idx == -1:
                gen_len = '0.0*'
            elif d_idx == 0:
                gen_len = str(abs(roundSig(gpos-gen_list[0],args.round)))+'*'
            else:
                gen_len = str(abs(roundSig(gpos-gen_list[d_idx-1],args.round)))
            out_string += ':'+gen_len
    #exit()
    out_string += '\n'
    return out_string

def positionCondition(outpos_list,outpos_idx,pos,revpos_flag):
    if outpos_idx == len(outpos_list):
        return False
    if revpos_flag and outpos_list[outpos_idx] >= pos:
        return True
    elif not revpos_flag and outpos_list[outpos_idx] < pos:
        return True
    return False

def getmsh(args):
    parser = createParser()
    args = parser.parse_args(args)
    if args.singleton_phase and not args.k1:
        raise Exception("Singleton phasing requires --k1")

    gen_flag = False
    sub_flag = False
    pos_flag = False
    if args.genname is not None:
        gf = open(args.genname,'r')
        l1,l2 = getGenMap(gf,args.genidx[0],args.genidx[1],squish=args.squish)
        gen_flag = True

    idx_list = None
    if args.subname is not None:
        #sf = open(args.subname,'r')
        #idx_list = [int(i.strip()) for i in sf]
        sub_list = readsubfile(args.subname)
        sub_flag = True
    if args.posname is not None:
        posf = open(args.posname,'r')
        outpos_list = [int(l.strip()) for l in posf]
        if args.revpos:
            outpos_list = outpos_list[::-1]
        outpos_idx = 0
        pos_flag = True
    compressed_mode = False
    if args.vcfname[-3:] == '.gz':
        compressed_mode = True
        vcf = gzip.open(args.vcfname,'r')
    else:
        vcf = open(args.vcfname,'r')
    #k = args.k_val
    k_vals = None
    #if args.k_all:
    #    k_vals = 
    #k_mode = (args.k_all or args.k_range is not None or args.k_list is not None)
    k_mode = (args.k_all or args.k_range is not None or args.k_list is not None or args.k1)
    snp_count = 0
    pos_list = []
    gen_list = None
    a = None
    d = None
    current_chrom = None
    if gen_flag:
        gen_list = []
    if args.outname is not None:
        outfn = args.outname
    else:
        raise Exception("--out option required for msh_to_vcf")
    if pos_flag and k_mode:
        raise Exception("K-mode and position modes are incompatible")
    if outfn[-3:] == '.gz':
        compress_out = True
        outf = gzip.open(outfn,"w")
    else:
        compress_out = False
        outf = open(outfn,"w")
    for line in vcf:
        if compressed_mode:
            line = line.decode()
        if line[0] == '#' or line[0] == '\n':
            continue
        la = line.strip().split()
        if 'CPG_TAG' in la[7]:
            continue
        if len(la[3]) != 1 or len(la[4]) != 1:
            continue
        if current_chrom is None:
            current_chrom = la[0]
        if la[0] != current_chrom:
            raise Exception("Source VCF cannot contain multiple chromosomes")
        noninf_pos = None
        noninf_gen = None
        singleton_idx = None
        k_idxlist = None
        k_pos = None
        k_gen = None
        is_singleton = False
        if snp_count == 0 and sub_flag:
            idx_list = subsampToIdx(la,sub_list)
        alleles = splitAlleles(la,idx_list)
        if snp_count == 0 and k_mode:
            if args.k_all:
                k_vals = [i for i in range(len(alleles))]
            elif args.k_range is not None:
                k_vals = [i for i in range(args.k_range[0],args.k_range[1]+1)]
            elif args.k_list is not None:
                k_vals = args.k_list
            elif args.k1:
                k_vals = [1]
        if args.dt_exp is not None and args.dt_exp[0] == int(la[1]):
            alleles[args.dt_exp[1]] = 1
        ac = sum(alleles)
        if len(alleles) == 0 or ac == 0 or ac == len(alleles):
            continue
        if ac == 1 or ac == len(alleles)-1:
            #if not args.singleton and not k_mode:
            is_singleton = True
            if not k_mode:
##                if not args.inc_sing:  ## not clear what this does
                if args.exc_sing:  ## changed to this JH 7/11/2019
                    continue
            #elif not k_mode:
            #    noninf_pos = int(la[1])
            #    singleton_idx = getSingletonIdx(alleles)
            #    if gen_flag:
            #        noninf_gen = float(getGenPos(noninf_pos,l1,l2))
            #if args.k_all and ac == 1 and args.exc_sing:
            #    continue
        #do_k = args.k_all and ac != 1
        #if args.k_all:
        #if do_k:
        #if args.k_all:
        #    k_idxlist = getKIdxList(alleles,ac)
        #    k_pos = int(la[1])
        #    if gen_flag:
        #        k_gen = float(getGenPos(k_pos,l1,l2))
        #if args.k_all or (args.k_range is not None and ac >= args.k_range[0] and ac <= args.k_range[1]):
        if k_vals is not None and ac in k_vals:
            k_idxlist = getKIdxList(alleles,ac)
            k_pos = int(la[1])
            if gen_flag:
                k_gen = float(getGenPos(k_pos,l1,l2))

        if snp_count == 0:
            sample_count = len(alleles)
            a_prev = [i for i in range(sample_count)]
            d_prev = [0 for i in range(sample_count)]
        if pos_flag:
            while positionCondition(outpos_list,outpos_idx,int(la[1]),args.revpos):
                opos = outpos_list[outpos_idx]
                out_range = range(sample_count)
                gpos = None
                if gen_flag:
                    gpos = float(getGenPos(opos,l1,l2))
                out_string = getMshString(args,a,d,out_range,pos_list,gen_list,opos,gpos,sample_count,k_idxlist)
                writeToFile(outf,out_string,compress_out)
                outpos_idx += 1
            if outpos_idx == len(outpos_list):
                break
        #if (args.singleton and noninf_pos is not None):
            #Singleton mode, singletons included in cutoffs
            #Current snp is singleton
            #Outputs lengths starting at current position
        #    out_range = [singleton_idx,singleton_idx+1]
        #    out_string = getMshString(args,a,d,out_range,pos_list,gen_list,noninf_pos,noninf_gen,sample_count,k_idxlist)
        #    writeToFile(outf,out_string,compress_out)
        if k_idxlist is not None:
            if args.dt_exp is not None and args.dt_exp[0] != int(la[1]):
                continue
            if args.k1 and args.singleton_phase:
                s_i = k_idxlist[0]
                out_range = ([s_i-1,s_i] if s_i%2==1 else [s_i,s_i+1])
            else:
                out_range = k_idxlist
            out_string = getMshString(args,a,d,out_range,pos_list,gen_list,k_pos,k_gen,sample_count,k_idxlist)
            writeToFile(outf,out_string,compress_out)
            if args.dt_exp is not None and args.dt_exp[0] == int(la[1]):
                break
        #if not args.singleton or noninf_pos is None or not args.exc_sing:  ## changed from or args.inc_sing  JH 7/11/2019
        if not k_mode or not is_singleton or not args.exc_sing:
            #Only skips when in singleton mode when singleton is hit
            #and shouldn't affect a/d
            pos_list.append(int(la[1]))
            if gen_flag:
                gen_list.append(float(getGenPos(int(la[1]),l1,l2)))
            a,d = getVectors(a_prev,d_prev,alleles,snp_count)
            a_prev = a
            d_prev = d
            snp_count += 1
        if not args.singleton and not pos_flag and not k_mode:
            #Standard alledges use
            if args.alledges_sing_only and ac != 1:
                continue
            out_range = range(sample_count)
            gpos = (None if gen_list is None else gen_list[-1])
            out_string = getMshString(args,a,d,out_range,pos_list,gen_list,pos_list[-1],gpos,sample_count,k_idxlist)
            writeToFile(outf,out_string,compress_out)

    outf.close()
    vcf.close()
    return outfn

if __name__ == "__main__":
    getmsh(sys.argv[1:])
