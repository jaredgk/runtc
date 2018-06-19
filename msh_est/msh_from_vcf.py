import sys
import argparse
import gzip
import math


def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gen",dest="genname")
    parser.add_argument("--vcf",dest="vcfname")
    parser.add_argument("--sub",dest="subname")
    parser.add_argument("--out",dest="outname")
    parser.add_argument("--gen-idx",dest="genidx",type=int,default=0)
    parser.add_argument("--nosquish",dest="squish",action="store_false")
    parser.add_argument("--round",dest="round",type=int,default=-1)
    return parser

def splitAllelesAll(la):
    alleles = []
    ac = 0
    an = 0
    for i in range(9,len(la)):
        try:
            g1 = int(la[i][0])
            try:
                g2 = int(la[i][2])
            except Exception:
                print (la)
        except ValueError:
            return []
        alleles.append(g1)
        alleles.append(g2)
        ac += g1
        ac += g2
        an += 2
    if ac == 1 or ac == an - 1:
        return []
    return alleles

def writeToFile(outf,write_line,compress_out):
    if compress_out:
        write_line = write_line.encode()
    outf.write(write_line)

def splitAllelesSub(la,idx_list):
    alleles = []
    ac = 0
    an = 0
    for i in range(len(idx_list)):
        reg = (idx_list[i]%2)*2
        f_idx = (idx_list[i]//2)+9
        try:
            geno = int(la[f_idx][reg])
        except ValueError:
            return []
        alleles.append(geno)
        ac += geno
        an += 1
    if ac <= 1 or ac >= an - 1:
        return []
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

def msh(a,d,k,pos_list):
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
            y_msh[i] = -2
        else:
            y_msh[i] = abs(pos_list[-1] - pos_list[c_idx-1])
    for a_i,a_v in enumerate(a):
        site_msh[a_v] = y_msh[a_i]
    return site_msh

def msh_gen(a,d,k,gen_list):
    l = len(a)
    g_msh = [0.0 for i in range(l)]
    site_msh = [0.0 for i in range(l)]
    for i in range(l):
        if i == l-1:
            c_idx = d[l-1]
        elif i == 0:
            c_idx = d[1]
        else:
            c_idx = min(d[i],d[i+1])
        if c_idx == 0:
            g_msh[i] = -2
        else:
            g_msh[i] = abs(gen_list[-1] - gen_list[c_idx-1])
    for a_i,a_v in enumerate(a):
        site_msh[a_v] = g_msh[a_i]
    return site_msh

def roundSig(f,n):
    """
    Rounds chi to n significant digits

    """
    if n == -1:
        return f
    if n < 1:
        raise Exception("N value %d is not valid" % (n))
    sign = (-1 if f < 0 else 1)
    return sign*round(f,-int(math.floor(math.log10(abs(f))))+(n-1))

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


def getmsh(args):
    parser = createParser()
    args = parser.parse_args(args)

    gen_flag = False
    sub_flag = False
    if args.genname is not None:
        gf = open(args.genname,'r')
        l1,l2 = getGenMap(gf,idx=args.genidx,squish=args.squish)
        gen_flag = True

    idx_list = None
    if args.subname is not None:
        #sf = open(args.subname,'r')
        #idx_list = [int(i.strip()) for i in sf]
        idx_list = readsubfile(args.subname)
        sub_flag = True
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
    if gen_flag:
        gen_list = []
    if args.outname is not None:
        outfn = args.outname
    else:
        raise Exception("--out option required for msh_to_vcf")
    if outfn[-3:] == '.gz':
        compress_out = True
        outf = gzip.open(outfn,"w")
    else:
        compress_out = False
        outf = open(outfn,"w")
    for line in vcf:
        if compressed_mode:
            line = line.decode()
        if line[0] == '#':
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
        alleles = splitAlleles(la,idx_list)
        if len(alleles) == 0:
            continue
        if k == 0:
            sample_count = len(alleles)
            a_prev = [i for i in range(sample_count)]
            d_prev = [0 for i in range(sample_count)]
        #sys.stderr.write(str(k)+'\n')
        pos_list.append(int(la[1]))
        if gen_flag:
            gen_list.append(float(getGenPos(int(la[1]),l1,l2)))
        a,d = getVectors(a_prev,d_prev,alleles,k)
        msh_vec = msh(a,d,k+1,pos_list)
        if gen_flag:
            g_vec = msh(a,d,k+1,gen_list)
        #tv = [str(ll) for ll in sorted(msh_vec)]
        #outf.write(str(pos_list[-1]))
        writeToFile(outf,str(pos_list[-1]),compress_out)
        if gen_flag:
            writeToFile(outf,'\t'+str(roundSig(gen_list[-1],args.round)),compress_out)
            #outf.write('\t'+str(gen_list[-1]))
        for i in range(len(msh_vec)):
            writeToFile(outf,'\t'+str(msh_vec[i]),compress_out)
            #outf.write('\t'+str(msh_vec[i]))
            if gen_flag:
                writeToFile(outf,':'+str(roundSig(g_vec[i],args.round)),compress_out)
                #outf.write(':'+str(g_vec[i]))
        writeToFile(outf,'\n',compress_out)
        #outf.write('\n')

        a_prev = a
        d_prev = d
        k += 1
    outf.close()
    return outfn

if __name__ == "__main__":
    getmsh(sys.argv[1:])
