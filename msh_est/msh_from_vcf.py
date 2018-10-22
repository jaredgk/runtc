import sys
import argparse
import gzip
import math
import gc
#from memory_profiler import profile


def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen",dest="genname")
    parser.add_argument("--vcf",dest="vcfname")
    parser.add_argument("--sub",dest="subname")
    parser.add_argument("--out",dest="outname")
    parser.add_argument("--gen-idx",dest="genidx",type=int,default=0)
    parser.add_argument("--nosquish",dest="squish",action="store_false")
    parser.add_argument("--round",dest="round",type=int,default=-1)
    parser.add_argument("--singleton",dest="singleton",action="store_true")
    parser.add_argument("--positions",dest="posname",type=str)
    parser.add_argument("--include-singletons",dest="inc_sing",action="store_true")
    parser.add_argument("--revpos",dest="revpos",action="store_true")
    return parser

def splitAllelesAll(la):
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
    ac = 0
    an = 0
    #for i in range(len(idx_list)):
    #    reg = (idx_list[i]%2)*2
    #    f_idx = (idx_list[i]//2)+9
    #    try:
    #        geno = int(la[f_idx][reg])
    #    except ValueError:
    #        return []
    #    alleles.append(geno)
    #    ac += geno
    #    an += 1
    for i in range(len(idx_list)):
        ip = idx_list[i]
        laa = la[ip[0]].split(':')[0]
        gts = laa.replace('|','/').split('/')
        try:
            geno = int(gts[ip[1]])
        except ValueError:
            return []
        alleles.append(geno)

    #if ac == 1 or ac == an - 1:
    #    return []
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
        return [-2 for i in range(sample_count)]
    l = len(a)
    y_msh = [0 for i in range(l)]
    site_msh = [0 for i in range(l)]
    #if noninf_pos is not None:
    #    pos = noninf_pos
    #else:
    #    pos = pos_list[-1]
    for i in range(l):
        if i == l-1:
            c_idx = d[l-1]
        elif i == 0:
            c_idx = d[1]
        else:
            c_idx = min(d[i],d[i+1])
        if c_idx == 0:
            y_msh[i] = -2
        else:
            y_msh[i] = abs(pos - pos_list[c_idx-1])
    for a_i,a_v in enumerate(a):
        site_msh[a_v] = y_msh[a_i]
    return site_msh


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

def parseGenLine(la,offset):
    a = float(la[0])
    b = float(la[1+offset])
    return a,b

def getGenMap(f,idx=0,squish=False):
    l1 = []
    l2 = []
    for line in f:
        if sum(c.isdigit() for c in line) < sum(c.isalpha() for c in line): # skip a line with more letters than numbers # first line may be column headers
            continue
        la = line.strip().split()
        if la[0][0].isalpha():
            if len(la)==4: ## only work with columns 1 and 3
                a,b = parseGenLine([la[1],la[3]],idx)
            else:
                a,b = parseGenLine(la[1:],idx)
        else:
            if len(la)==3: ## only work with columns 0 and 2
                a,b = parseGenLine([la[0],la[2]],idx)
            else:
                a,b = parseGenLine(la,idx)

        #if b != 0:
        #For future: Use last position with 0.0 cm instead of first
        if not squish or len(l2) == 0 or (len(l2) > 0 and l2[-1] != b):
            l1.append(a)
            l2.append(b)
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
    s = sum(alleles)
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

def getMshString(args,a,d,out_range,pos_list,gen_list,pos,gpos,sample_count):
    out_string = ''
    msh_vec = msh(a,d,pos_list,pos,sample_count)
    gen_flag = (gen_list is not None)
    if gen_flag:
        g_vec = msh(a,d,gen_list,gpos,sample_count)
    out_string += str(pos)
    if gen_flag:
        out_string += '\t'+str(roundSig(gpos,args.round))
    for i in out_range:
        out_string += '\t'+str(msh_vec[i])
        if gen_flag:
            out_string += ':'+str(roundSig(g_vec[i],args.round))
    out_string += '\n'
    return out_string

def positionCondition(outpos_list,outpos_idx,pos,revpos_flag):
    if outpos_idx == len(outpos_list):
        return False
    if revpos_flag and outpos_list[outpos_idx] >= pos:
        return True
    elif not revpos_flag and outpos_list[outpos_idx] <= pos:
        return True
    return False

def getmsh(args):
    parser = createParser()
    args = parser.parse_args(args)

    gen_flag = False
    sub_flag = False
    pos_flag = False
    if args.genname is not None:
        gf = open(args.genname,'r')
        l1,l2 = getGenMap(gf,idx=args.genidx,squish=args.squish)
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
    a_mat = []
    d_mat = []
    m_mat = []

    k = 0
    pos_list = []
    gen_list = None
    a = None
    d = None
    if gen_flag:
        gen_list = []
    if args.outname is not None:
        outfn = args.outname
    else:
        raise Exception("--out option required for msh_to_vcf")
    if pos_flag and args.singleton:
        raise Exception("Singleton and position modes are incompatible")
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
        try:
            if 'CPG_TAG' in la[7]:
                continue
            lastline = la
        except Exception:
            pass
        if len(la[3]) != 1 or len(la[4]) != 1:
            continue
        noninf_pos = None
        noninf_gen = None
        singleton_idx = None
        if k == 0 and sub_flag:
            idx_list = subsampToIdx(la,sub_list)
        alleles = splitAlleles(la,idx_list)
        if len(alleles) == 0:
            continue
        if sum(alleles) == 1 or sum(alleles) == len(alleles)-1:
            if not args.singleton:
                if not args.inc_sing:
                    continue
                #else:
                #    singleton_idx = getSingletonIdx(alleles)
            else:
                noninf_pos = int(la[1])
                singleton_idx = getSingletonIdx(alleles)
                if gen_flag:
                    noninf_gen = float(getGenPos(noninf_pos,l1,l2))
        if k == 0:
            sample_count = len(alleles)
            a_prev = [i for i in range(sample_count)]
            d_prev = [0 for i in range(sample_count)]
            #msh_vec = [-2 for ii in range(sample_count)]
            #if gen_flag:
            #    g_vec = [-2.0 for ii in range(sample_count)]
        if pos_flag:
            #if outpos_list[outpos_idx] <= int(la[1]):
            #while outpos_idx < len(outpos_list) and outpos_list[outpos_idx] <= int(la[1]):
            while positionCondition(outpos_list,outpos_idx,int(la[1]),args.revpos):
                opos = outpos_list[outpos_idx]
                out_range = range(sample_count)
                gpos = None
                if gen_flag:
                    gpos = float(getGenPos(opos,l1,l2))
                out_string = getMshString(args,a,d,out_range,pos_list,gen_list,opos,gpos,sample_count)
                writeToFile(outf,out_string,compress_out)
                outpos_idx += 1
            if outpos_idx == len(outpos_list):
                break
        if args.singleton and args.inc_sing and noninf_pos is not None:
            out_range = [singleton_idx,singleton_idx+1]
            out_string = getMshString(args,a,d,out_range,pos_list,gen_list,noninf_pos,noninf_gen,sample_count)
            writeToFile(outf,out_string,compress_out)
        if not args.singleton or noninf_pos is None or args.inc_sing:
            pos_list.append(int(la[1]))
            if gen_flag:
                gen_list.append(float(getGenPos(int(la[1]),l1,l2)))
            a,d = getVectors(a_prev,d_prev,alleles,k)
            a_prev = a
            d_prev = d
            k += 1
        if args.singleton and noninf_pos is not None and not args.inc_sing:
            out_range = [singleton_idx,singleton_idx+1]
            out_string = getMshString(args,a,d,out_range,pos_list,gen_list,noninf_pos,noninf_gen,sample_count)
            writeToFile(outf,out_string,compress_out)
            #pos = noninf_pos
        if not args.singleton and not pos_flag:
            out_range = range(sample_count)
            gpos = (None if gen_list is None else gen_list[-1])
            out_string = getMshString(args,a,d,out_range,pos_list,gen_list,pos_list[-1],gpos,sample_count)
            writeToFile(outf,out_string,compress_out)
        #if args.singleton and args.inc_sing:
        #    out_range = []
        #if not args.singleton or noninf_pos is not None:
        #    if args.singleton:
        #        out_range = [singleton_idx,singleton_idx+1]
        #    else:
        #        out_range = range(sample_count)
        #    if noninf_pos is not None:
        #        pos = noninf_pos
        #        if gen_flag:
        #            gpos = noninf_gen
        #    else:
        #        pos = pos_list[-1]
        #        if gen_flag:
        #            gpos = gen_list[-1]
        #    msh_vec = msh(a,d,pos_list,pos)
        #    if gen_flag:
        #        g_vec = msh(a,d,gen_list,gpos)
        #    writeToFile(outf,str(pos),compress_out)
        #    if gen_flag:
        #        genpos = float(getGenPos(pos,l1,l2))
        #     writeToFile(outf,'\t'+str(roundSig(genpos,args.round)),compress_out)
        #    for i in out_range:
        #        writeToFile(outf,'\t'+str(msh_vec[i]),compress_out)
        #        if gen_flag:
        #            writeToFile(outf,':'+str(roundSig(g_vec[i],args.round)),compress_out)
        #    writeToFile(outf,'\n',compress_out)


    outf.close()
    return outfn

if __name__ == "__main__":
    getmsh(sys.argv[1:])
